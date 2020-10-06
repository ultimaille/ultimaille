#include <iostream>
#include <cstdlib>
#include "ultimaille/mesh_io.h"
#include "ultimaille/surface.h"
#include "ultimaille/attributes.h"
#include "ultimaille/range.h"

#include <OpenNL_psm/OpenNL_psm.h>

mat<2,2> Bi(const double theta) {
    return {{ {cos(theta), cos(M_PI/2+theta)}, {sin(theta), sin(M_PI/2+theta)} }};
}

mat<2,2> Rij(const double thi, const double thj) {
    mat<2,2> R90 = {{ {0,-1}, {1,0} }};
    mat<2,2> bi = Bi(thi);
    mat<2,2> bj = Bi(thj);

    mat<2,2> curr = mat<2,2>::identity();
    mat<2,2> rij = curr;
    double best_norm = (bi - bj).norm();

    for (int k=0; k<3; k++) {
        curr = R90*curr;
        double norm = (bj - curr*bi).norm();  // TODO why?? it DOES correspond to the definiton in section 2.2 from the PGP3D paper, but i do not understand it
        if (norm<best_norm) {
            best_norm = norm;
            rij = curr;
        }
    }
    return rij;
}

vec2 gij(const double thi, const double thj, const vec2 edge) {
    return (Bi(thi).transpose() + Rij(thi, thj)*Bi(thj).transpose())*edge*.5;
}

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
    }

    PointAttribute<double> theta(m.points);
    FacetAttribute<int> sing(m);

    SurfaceConnectivity fec(m);
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

    { // interpolate the frame field
        for (int iter : range(5)) {
            nlNewContext();
            nlSolverParameteri(NL_NB_VARIABLES, 2*m.nverts());
            nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
            nlBegin(NL_SYSTEM);
            nlBegin(NL_MATRIX);

            for (int ci : range(m.ncorners())) {
                if (fec.opposite(ci)>=0 && fec.from(ci)>fec.to(ci)) continue;
                for (int d : range(2)) {
                    nlBegin(NL_ROW);
                    nlCoefficient(2*fec.to(ci)  + d,  1);
                    nlCoefficient(2*fec.from(ci)+ d, -1);
                    nlEnd(NL_ROW);
                }
            }
            if (iter) for (int vi : range(m.nverts())) {
                vec2 goal = { cos(4*theta[vi]), sin(4*theta[vi]) };
                for (int d : range(2)) {
                    nlRowScaling(.1);
                    nlBegin(NL_ROW);
                    nlCoefficient(2*vi+d, 1);
                    nlRightHandSide(goal[d]);
                    nlEnd(NL_ROW);
                }
            }

            for (int vi : range(m.nverts())) {
                if (!fec.is_border_vert(vi)) continue;
                vec2 goal = { cos(4*theta[vi]), sin(4*theta[vi]) };
                for (int d : range(2)) {
                    nlRowScaling(100);
                    nlBegin(NL_ROW);
                    nlCoefficient(2*vi+d,  1);
                    nlRightHandSide(goal[d]);
                    nlEnd(NL_ROW);
                }
            }

            nlEnd(NL_MATRIX);
            nlEnd(NL_SYSTEM);
            nlSolve();

            for (int vi : range(m.nverts()))
                theta[vi] = atan2(nlGetVariable(2*vi+1), nlGetVariable(2*vi))/4 ;// + M_PI/2 * (m.points[vi].x>5) + (rand()%4)*M_PI/2.;

            for (int fi : range(m.nfacets())) {
                mat<2,2> R = mat<2,2>::identity();
                for (int lvi : range(3))
                    R = Rij(theta[m.vert(fi, lvi)], theta[m.vert(fi, (lvi+1)%3)])*R;
                sing[fi] = ((R-mat<2,2>::identity()).norm()>1e-10);
            }
            nlDeleteContext(nlGetCurrent());
        }
    }

    CornerAttribute<vec2> cij(m);
    { // curl correction step
        nlNewContext();
        nlSolverParameteri(NL_NB_VARIABLES, 2*m.ncorners());
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlBegin(NL_SYSTEM);
        nlBegin(NL_MATRIX);

        for (int hi : range(m.ncorners())) {
            for (int d : range(2)) { // \sum_{ij} ||c_{ij}||^2
                nlRowScaling(.1);
                nlBegin(NL_ROW);
                nlCoefficient(2*hi + d,  1);
                nlEnd(NL_ROW);
            }

            int opp = fec.opposite(hi);
            if (opp>=0) { // symmetry of the curl correction
                int i = fec.from(hi);
                int j = fec.to(hi);
                mat<2,2> Rij_ = Rij(theta[i], theta[j]);
                for (int d : range(2)) {
                    nlRowScaling(100);
                    nlBegin(NL_ROW);
                    nlCoefficient(2*hi+d, 1);
                    nlCoefficient(2*opp+0, Rij_[d][0]);
                    nlCoefficient(2*opp+1, Rij_[d][1]);
                    nlEnd(NL_ROW);
                }
            } else {
                vec3 edge = m.points[fec.to(hi)] - m.points[fec.from(hi)];
                vec2 n = {edge.y, -edge.x};

                int vi = fec.from(hi);
                mat<2,2> BiT = Bi(theta[vi]).transpose();
                int tmp = (std::abs(BiT[0]*n) > std::abs(BiT[1]*n));

                nlRowScaling(100);
                nlBegin(NL_ROW);
                nlCoefficient(2*hi+0,    tmp);
                nlCoefficient(2*hi+1,  1-tmp);
                nlEnd(NL_ROW);
            }
        }

        for (int fi : range(m.nfacets())) {
            if (sing[fi]) continue;
            int hij = m.facet_corner(fi, 0);
            int hjk = fec.next(hij);
            int hki = fec.prev(hij);
            int i = fec.from(hij);
            int j = fec.from(hjk);
            int k = fec.from(hki);
            mat<2,2> Rij_ = Rij(theta[i], theta[j]);
            mat<2,2> Rik_ = Rij(theta[i], theta[k]);
            vec3 eij = m.points[j] - m.points[i];
            vec3 ejk = m.points[k] - m.points[j];
            vec3 eki = m.points[i] - m.points[k];

            vec2 gij_ = gij(theta[i], theta[j], {eij.x, eij.y});
            vec2 gjk_ = gij(theta[j], theta[k], {ejk.x, ejk.y});
            vec2 gki_ = gij(theta[k], theta[i], {eki.x, eki.y});
            for (int d : range(2)) {
                nlRowScaling(100);
                nlBegin(NL_ROW);
                nlCoefficient(2*hij+d,  1);
                nlCoefficient(2*hjk+0,  Rij_[d][0]);
                nlCoefficient(2*hjk+1,  Rij_[d][1]);
                nlCoefficient(2*hki+0,  Rik_[d][0]);
                nlCoefficient(2*hki+1,  Rik_[d][1]);
                double rhs = (-gij_[d] - (Rij_[d][0]*gjk_[0]+Rij_[d][1]*gjk_[1]) - (Rik_[d][0]*gki_[0]+Rik_[d][1]*gki_[1]));
                nlRightHandSide(rhs);
                nlEnd(NL_ROW);
            }
        }

        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);
        nlSolve();

        for (int hi : range(m.ncorners())) {
            int i = fec.from(hi);
            int j = fec.to(hi);
            vec3 eij = m.points[j] - m.points[i];
            vec2 gij_ = gij(theta[i], theta[j], {eij.x, eij.y});

            cij[hi] = vec2{nlGetVariable(2*hi+0), nlGetVariable(2*hi+1)};
            cij[hi] = std::min(1., .35*gij_.norm()/cij[hi].norm())*cij[hi];
        }
        nlDeleteContext(nlGetCurrent());
    }


    CornerAttribute<vec2> uti(m);
    PointAttribute<vec2> ui(m.points);

    { // PGP
        double av_length = average_edge_length(m);
        std::cerr << "average edge length is " << av_length << std::endl;
        const double scale = 0.4/av_length;
        for (int iter : range(5)) {
            nlNewContext();
            nlSolverParameteri(NL_NB_VARIABLES, 4*m.nverts());
            nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
            nlBegin(NL_SYSTEM);

            for (int ci : range(m.ncorners())) {
                if (fec.opposite(ci)>=0) continue;

                vec3 edge = m.points[fec.to(ci)] - m.points[fec.from(ci)];
                vec2 n = {edge.y, -edge.x};

                int vi = fec.from(ci);
                mat<2,2> BiT = Bi(theta[vi]).transpose();

                if (std::abs(BiT[0]*n) > std::abs(BiT[1]*n)) {
                    nlSetVariable(4*vi+0, 1); // u is integer
                    nlSetVariable(4*vi+1, 0);
                    nlLockVariable(4*vi+0);
                    nlLockVariable(4*vi+1);
                } else {
                    nlSetVariable(4*vi+2, 1); // v is integer
                    nlSetVariable(4*vi+3, 0);
                    nlLockVariable(4*vi+2);
                    nlLockVariable(4*vi+3);
                }
            }

            nlBegin(NL_MATRIX);

            /*
            if (iter) for (int vi : range(m.nverts())) {
                double old[4] = { cos(2*M_PI*ui[vi].x), sin(2*M_PI*ui[vi].x), cos(2*M_PI*ui[vi].y), sin(2*M_PI*ui[vi].y) };
                for (int z : range(4)) {
                    nlRowScaling(.1);
                    nlBegin(NL_ROW);
                    nlCoefficient(2*vi+z, 1);
                    nlRightHandSide(old[z]);
                    nlEnd(NL_ROW);
                }
            }
            */

            for (int ci : range(m.ncorners())) {
                int i = fec.from(ci);
                int j = fec.to(ci);
                vec3 edge = m.points[fec.to(ci)] - m.points[fec.from(ci)];

                mat<2,2> Rij_ = Rij(theta[i], theta[j]);
                vec2 gij_ = scale*(gij(theta[i], theta[j], {edge.x, edge.y}) + cij[ci]);

                for (int d : range(2)) {
                    assert(Rij_[d][0] == 0 || Rij_[d][1]==0);
                    double rdij = (Rij_[d][0] ? 0 : 1);
                    double sdij = Rij_[d][0] + Rij_[d][1];

                    nlBegin(NL_ROW);
                    nlCoefficient(4*i+d*2+0,  cos(2.*M_PI*gij_[d]));
                    nlCoefficient(4*i+d*2+1, -sin(2.*M_PI*gij_[d]));
                    nlCoefficient(4*j+2*rdij, -1);
                    nlEnd(NL_ROW);

                    nlBegin(NL_ROW);
                    nlCoefficient(4*i+d*2+1, cos(2.*M_PI*gij_[d]));
                    nlCoefficient(4*i+d*2+0, sin(2.*M_PI*gij_[d]));
                    nlCoefficient(4*j+2*rdij+1, -sdij);
                    nlEnd(NL_ROW);
                }
            }
            nlEnd(NL_MATRIX);
            nlEnd(NL_SYSTEM);
            nlSolve();

            for (int i : range(m.nverts())) {
                ui[i][0] = atan2(nlGetVariable(4*i+1), nlGetVariable(4*i+0))/(2.*M_PI);
                ui[i][1] = atan2(nlGetVariable(4*i+3), nlGetVariable(4*i+2))/(2.*M_PI);
            }
        }

        for (int hi : range(m.ncorners())) // transfer tex coords from verts to corners
            uti[hi] = ui[fec.from(hi)];

        for (int fi : range(m.nfacets())) {
            int hr = m.facet_corner(fi, 2);
            int i = fec.from(hr);

            for (int e : range(2)) {
                int he = m.facet_corner(fi, e);
                int j = fec.from(he);
                uti[he] = Rij(theta[i], theta[j])*uti[he];

                vec3 edge = fec.m.points[j] - m.points[i];
                vec2 gij_ = scale*(gij(theta[i], theta[j], {edge.x, edge.y}) + cij[he]);

                for (int d : range(2)) {
                    while (uti[hr][d] + gij_[d] - uti[he][d] >  .5) uti[he][d] += 1;
                    while (uti[hr][d] + gij_[d] - uti[he][d] < -.5) uti[he][d] -= 1;
                }
            }
        }
        nlDeleteContext(nlGetCurrent());
    }
    write_geogram("pgp.geogram", m, { {{"theta", theta.ptr}, {"ui", ui.ptr}}, {{"sing", sing.ptr}}, {{"uti", uti.ptr}} });

    return 0;
}

