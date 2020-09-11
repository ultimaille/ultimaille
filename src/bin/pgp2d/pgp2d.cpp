#include <iostream>
#include <cstdlib>
#include "ultimaille/mesh_io.h"
#include "ultimaille/surface.h"
#include "ultimaille/attributes.h"
#include "ultimaille/range.h"

#include <OpenNL_psm/OpenNL_psm.h>

vec3 rotate_vector_around_axis(const vec3 axis, const double angle, const vec3 v) {
    vec3 u = vec3(axis).normalize();
    double c = cos(angle);
    double s = sin(angle);
    return vec3((c+u.x*u.x*(1-c))*v.x     + (u.x*u.y*(1-c)-u.z*s)*v.y + (u.x*u.z*(1-c)+u.y*s)*v.z,
                (u.x*u.y*(1-c)+u.z*s)*v.x + (c+u.y*u.y*(1-c))*v.y     + (u.y*u.z*(1-c)-u.x*s)*v.z,
                (u.x*u.z*(1-c)-u.y*s)*v.x + (u.z*u.y*(1-c)+u.x*s)*v.y + (c+u.z*u.z*(1-c))*v.z     );
}

mat<3,3> Bi(PointAttribute<vec3> &f0, PointAttribute<vec3> &f1, int i) {
    mat<3,3> m = mat<3,3>::identity();
    m[0] = f0[i];
    m[1] = f1[i];
    return m;
}

mat<3,3> Rij(PointAttribute<vec3> &f0, PointAttribute<vec3> &f1, int i, int j) {
    mat<3,3> m = mat<3,3>::identity();
    mat<3,3> best_m = m;
    double best_norm = (Bi(f0, f1, i) - Bi(f0, f1, j)).norm();

    mat<3,3> R = mat<3,3>::identity();
    R[0][0] = R[1][1] = 0;
    R[0][1] = -1;
    R[1][0] = 1;

    for (int k=0; k<3; k++) {
        m = m*R;
        double norm = (Bi(f0, f1, i) - Bi(f0, f1, j)*m).norm();
        if (norm<best_norm) {
            best_norm = norm;
            best_m = m;
        }
    }
    return best_m;
}

vec3 gij(PointAttribute<vec3> &f0, PointAttribute<vec3> &f1, Triangles &m, int i, int j) {
    return (Bi(f0, f1, i) + Rij(f0, f1, i, j)*Bi(f0, f1, j))*(m.points[j] - m.points[i])*.5;
}


int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }
    nlInitialize(argc, argv);


    Triangles m;
    {
        Polygons pm;
        read_wavefront_obj(argv[1], pm);
        pm.extract_triangles(m);
        assert(m.nfacets() == pm.nfacets());
        for (vec3 &p : *m.points.data) p.z = 0; // make sure it is 2D
    }

    PointAttribute<vec3> f0(m.points), f1(m.points), a(m.points);
    PointAttribute<double> theta(m.points), constraints(m.points);
    FacetAttribute<int> sing(m);

    { // init the attributes
        for (int vi : range(m.nverts())) {
            theta[vi] = 0;
            a[vi]  = {1, 0, 0};
            f0[vi] = rotate_vector_around_axis(vec3(0,0,1), theta[vi], vec3(1,0,0));
            f1[vi] = rotate_vector_around_axis(vec3(0,0,1), theta[vi], vec3(0,1,0));
        }
    }

    SurfaceConnectivity fec(m);
    { // compute the constraints
        for (int vi : range(m.nverts())) {
            if (!fec.is_border_vert(vi)) continue;
            vec3 n = {0,0,0};

            int ci = fec.v2c[vi];
            do {
                if (fec.opposite(ci)<0) {
                    vec3 hv = m.points[fec.to(ci)] - m.points[fec.from(ci)];
                    if (n.norm2()<1e-10) {
                        n = hv;
                    } else {
                        double dist = -1;
                        vec3 best = hv;
                        for (int i=0; i<4; i++) {
                            vec3 tv = rotate_vector_around_axis(vec3(0,0,1), M_PI/2.*i, hv);
                            double tdist = (tv-n).norm2();
                            if (dist<0 || tdist<dist) {
                                dist = tdist;
                                best = tv;
                            }
                        }
                        n += best;
                    }
                }
                ci = fec.c2c[ci];
            } while(ci != fec.v2c[vi]);

            n.normalize();
            double tmp = atan2(n.y, n.x);
            while (std::abs(tmp)>M_PI/4.) tmp += (tmp>0 ? -1. : 1.)*M_PI/2.;
            constraints[vi] = tmp;
        }

        for (int vi : range(m.nverts())) {
            if (!fec.is_border_vert(vi)) continue;
            theta[vi] = constraints[vi];
            a[vi]  = {cos(4.*theta[vi]), sin(4.*theta[vi]), 0};
            f0[vi] = rotate_vector_around_axis({0,0,1}, theta[vi], {1,0,0});
            f1[vi] = rotate_vector_around_axis({0,0,1}, theta[vi], {0,1,0});
        }
    }

    { // generate the frame field
        nlNewContext();
        nlSolverParameteri(NL_NB_VARIABLES, 2*m.nverts());
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlBegin(NL_SYSTEM);
        nlBegin(NL_MATRIX);

        for (int ci : range(m.ncorners())) {
            if (fec.opposite(ci)>=0 && fec.from(ci)>fec.to(ci)) continue;
            for (int d : range(2)) {
                nlBegin(NL_ROW);
                nlCoefficient(fec.to(ci)  + d*m.nverts(),  1);
                nlCoefficient(fec.from(ci)+ d*m.nverts(), -1);
                nlEnd(NL_ROW);
            }
        }

        for (int vi : range(m.nverts())) {
            if (!fec.is_border_vert(vi)) continue;
            vec3 goal = rotate_vector_around_axis({0,0,1}, 4.*constraints[vi], {1,0,0});
            for (int d=0; d<2; d++) {
                nlRowScaling(100);
                nlBegin(NL_ROW);
                nlCoefficient(vi + d*m.nverts(),  1);
                nlRightHandSide(goal[d]);
                nlEnd(NL_ROW);
            }
        }

        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);

        nlSolve();

        for (int vi : range(m.nverts())) {
            vec2 tmp = vec2(nlGetVariable(vi), nlGetVariable(vi+m.nverts())).normalize();
            a[vi] = {tmp.x, tmp.y, 0};
            theta[vi] = atan2(a[vi].y, a[vi].x)/4.;// + (rand()%4)*M_PI/2.;
            f0[vi] = rotate_vector_around_axis({0,0,1}, theta[vi], {1,0,0});
            f1[vi] = rotate_vector_around_axis({0,0,1}, theta[vi], {0,1,0});
        }
        for (int fi : range(m.nfacets())) {
            mat<3,3> R = mat<3,3>::identity();
            for (int lvi : range(3)) {
                R = R*Rij(f0, f1, m.vert(fi, lvi), m.vert(fi, (lvi+1)%3));
            }
            sing[fi] = ((R-mat<3,3>::identity()).norm()>1e-10);
        }

        nlDeleteContext(nlGetCurrent());
    }


    write_geogram("pgp.geogram", m, { {{"theta", theta.ptr}, {"constraints", constraints.ptr}}, {{"sing", sing.ptr}}, {} });

    return 0;
}

