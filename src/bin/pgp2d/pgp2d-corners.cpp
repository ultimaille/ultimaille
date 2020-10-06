#include <iostream>
#include <cstdlib>
#include "ultimaille/mesh_io.h"
#include "ultimaille/surface.h"
#include "ultimaille/attributes.h"
#include "ultimaille/range.h"

#include <OpenNL_psm/OpenNL_psm.h>


double average_edge_length(const Surface &m) {
    double sum = 0;
    int nb = 0;
    for (int f : range(m.nfacets())) {
        for (int lv : range(m.facet_size(f))) {
            vec3 a = m.points[m.vert(f, lv)];
            vec3 b = m.points[m.vert(f, (lv+1)%m.facet_size(f))];
            sum += (a-b).norm();
            nb++;
        }
    }
    return sum/nb;
}

void compute_frame_field(const Triangles &m, const SurfaceConnectivity &fec, FacetAttribute<double> &theta, FacetAttribute<mat<2,2>> &Bi) {
    { // compute the frame field at the boundary
        for (int c : range(m.ncorners())) {
            if (fec.opposite(c)>=0) continue;
            vec3 edge = m.points[fec.to(c)] - m.points[fec.from(c)];
            vec2 n = {edge.y, -edge.x};
            double angle = atan2(n.y, n.x);
            while (std::abs(angle)>M_PI/4) angle += (angle>0 ? -1 : 1)*M_PI/2;
            theta[fec.c2f[c]] = angle;
        }
    }

    { // interpolate the frame field
        for (int iter : range(10)) {
            nlNewContext();
            nlSolverParameteri(NL_NB_VARIABLES, 2*m.nfacets());
            nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
            nlBegin(NL_SYSTEM);
            nlBegin(NL_MATRIX);

            for (int c : range(m.ncorners())) { // smoothing term
                int opp = fec.opposite(c);
                if (opp<0 || opp<c) continue;
                for (int d : range(2)) {
                    nlBegin(NL_ROW);
                    nlCoefficient(2*fec.c2f[c]  + d,  1);
                    nlCoefficient(2*fec.c2f[opp]+ d, -1);
                    nlEnd(NL_ROW);
                }
            }

            for (int c : range(m.ncorners())) { // boundary constraint
                if (fec.opposite(c)>=0) continue;
                int f = fec.c2f[c];
                vec2 goal = { cos(4*theta[f]), sin(4*theta[f]) };
                for (int d : range(2)) {
                    nlRowScaling(100);
                    nlBegin(NL_ROW);
                    nlCoefficient(2*f+d,  1);
                    nlRightHandSide(goal[d]);
                    nlEnd(NL_ROW);
                }
            }

            if (iter) for (int f : range(m.nfacets())) { // enforce non-zero norm
                vec2 goal = { cos(4*theta[f]), sin(4*theta[f]) };
                for (int d : range(2)) {
                    nlRowScaling(.05*iter);
                    nlBegin(NL_ROW);
                    nlCoefficient(2*f+d, 1);
                    nlRightHandSide(goal[d]);
                    nlEnd(NL_ROW);
                }
            }

            nlEnd(NL_MATRIX);
            nlEnd(NL_SYSTEM);
            nlSolve();

            for (int f : range(m.nfacets()))
                theta[f] = atan2(nlGetVariable(2*f+1), nlGetVariable(2*f))/4;//  + M_PI/2 * (m.points[vi].x>5);// + (rand()%4)*M_PI/2.;

            nlDeleteContext(nlGetCurrent());
        }
    }

/*
    { // compute the frame field at the boundary
        for (int vi : range(m.nverts())) {
            if (!fec.is_border_vert(vi)) continue;
            vec2 n = {0,0};

            int ci = fec.v2c[vi];
            do {
                if (fec.opposite(ci)<0) {
                    vec3 edge = m.points[fec.to(ci)] - m.points[fec.from(ci)];
                    vec2 rote = {edge.y, -edge.x};
                    if (n.norm2()<1e-10) {
                        n = rote;
                    } else {
                        double dist = -1;
                        vec2 best = rote;
                        for (int i=0; i<4; i++) {
                            double tdist = (rote-n).norm2();
                            if (dist<0 || tdist<dist) {
                                dist = tdist;
                                best = rote;
                            }
                            rote = {rote.y, -rote.x};
                        }
                        n += best;
                    }
                }
                ci = fec.c2c[ci];
            } while (ci != fec.v2c[vi]);

            double angle = atan2(n.y, n.x);
            while (std::abs(angle)>M_PI/4) angle += (angle>0 ? -1 : 1)*M_PI/2;
            theta[vi] = angle;
        }
    }


/*
    std::vector<mat<2,2>> Rij(m.ncorners());

    for (int vi : range(m.nverts()))
        Bi[vi] = {{ {cos(theta[vi]), cos(M_PI/2+theta[vi])}, {sin(theta[vi]), sin(M_PI/2+theta[vi])} }};

    const mat<2,2> R90 = {{ {0,-1}, {1,0} }};
    for (int c : range(m.ncorners())) {
        int i=fec.from(c), j=fec.to(c), opp=fec.opposite(c);
        if (opp>=0 && i>j) continue;
        mat<2,2> curr = mat<2,2>::identity();
        Rij[c] = curr;
        double best_norm = (Bi[i] - Bi[j]).norm();
        for (int k=0; k<3; k++) {
            curr = R90*curr;
            double norm = (Bi[i] - curr*Bi[j]).norm();
            if (norm<best_norm) {
                best_norm = norm;
                Rij[c] = curr;
            }
        }
        if (opp>=0) {
            Rij[opp] = Rij[c].transpose();
        }
    }
*/
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
        for (vec3 &p : m.points) p.z = 0; // make sure it is 2D
        // TODO assert manifoldity
    }


    FacetAttribute<double> theta(m);
    FacetAttribute<mat<2,2>> Bi(m);

//  FacetAttribute<int> ffsing(m);

    std::cerr << "Computing fec...";
    SurfaceConnectivity fec(m);
    std::cerr << "ok\n";

    compute_frame_field(m, fec, theta, Bi);

    write_geogram("pgp.geogram", m, { {}, {{"theta", theta.ptr}}, {} });

    return 0;
}

