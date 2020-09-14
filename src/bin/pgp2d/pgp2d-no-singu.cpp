#include <iostream>
#include <cstdlib>
#include "ultimaille/mesh_io.h"
#include "ultimaille/surface.h"
#include "ultimaille/attributes.h"
#include "ultimaille/range.h"

#include <OpenNL_psm/OpenNL_psm.h>

#define STAR_LOOP(c, v, fec) for (int c=fec.v2c[v], __FLAG__=1; c!=fec.v2c[v] || __FLAG__; c=fec.c2c[c], __FLAG__=0)

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

void compute_frame_field(const Triangles &m, const SurfaceConnectivity &fec,
        PointAttribute<double> &theta,
        FacetAttribute<int> &ffsing,
        PointAttribute<mat<2,2>> &Bi,
        CornerAttribute<mat<2,2>> &Rij,
        CornerAttribute<int> &debugRij) {
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
        for (int iter : range(10)) {
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
                    nlRowScaling(.05*iter);
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
                theta[vi] = atan2(nlGetVariable(2*vi+1), nlGetVariable(2*vi))/4  + M_PI/2 * (m.points[vi].x>5);// + (rand()%4)*M_PI/2.;

            nlDeleteContext(nlGetCurrent());
        }
    }

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
                debugRij[c] = k+1;
            }
        }
        if (opp>=0) {
            Rij[opp] = Rij[c].transpose();
            debugRij[opp] = (4-debugRij[c])%4;
        }
    }

    for (int fi : range(m.nfacets())) {
        mat<2,2> R = mat<2,2>::identity();
        for (int lci : range(3))
            R = Rij[m.facet_corner(fi, lci)]*R;
        ffsing[fi] = ((R-mat<2,2>::identity()).norm()>1e-10);
    }
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

    PointAttribute<double> theta(m.points);
    FacetAttribute<int> ffsing(m);
    PointAttribute<mat<2,2>> Bi(m.points);
    CornerAttribute<mat<2,2>> Rij(m);
    CornerAttribute<int> debugRij(m);

    SurfaceConnectivity fec(m);

    compute_frame_field(m, fec, theta, ffsing, Bi, Rij, debugRij);

    CornerAttribute<vec2> gij(m);
    for (int c : range(m.ncorners())) {
        int i=fec.from(c), j=fec.to(c);//, opp=fec.opposite(c);
//        if (opp>=0 && i>j) continue;
        vec3 edge = m.points[j] - m.points[i];
        gij[c] = ((Bi[i] + Rij[c]*Bi[j])/2).transpose()*vec2(edge.x, edge.y);
//        if (opp>=0) gij[opp] = -(Bi[j].transpose()*Bi[i]*gij[c]);

//        if (opp>=0) gij[opp] = Rij[c]*(-gij[c]); // TODO WTF?
    }



    PointAttribute<int> debugBi(m.points);
    CornerAttribute<int> lci(m);
    for (int f : range(m.nfacets()))
        for (int i : range(3))
            lci[m.facet_corner(f, i)] = i;


    const mat<2,2> R90 = {{ {0,-1}, {1,0} }};
    for (int v : range(m.nverts())) {
        mat<2,2> curr = mat<2,2>::identity();
        double best_norm = (Bi[v] - curr).norm();
        for (int k=0; k<3; k++) {
            curr = R90*curr;
            double norm = (Bi[v] - curr).norm();  // TODO why?? it DOES correspond to the definiton in section 2.2 from the PGP3D paper, but i do not understand it
            if (norm<best_norm) {
                best_norm = norm;
                debugBi[v] = k+1;
            }
        }
    }


    CornerAttribute<vec2> cij(m);

/*
    { // curl correction step
        nlNewContext();
        nlSolverParameteri(NL_NB_VARIABLES, 2*m.ncorners());
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlBegin(NL_SYSTEM);
        nlBegin(NL_MATRIX);

        for (int c : range(m.ncorners())) {
            int i = fec.from(c), j = fec.to(c);
            for (int d : range(2)) { // min \sum_{ij} ||c_{ij}||^2
                nlRowScaling(.1);
                nlBegin(NL_ROW);
                nlCoefficient(2*c + d,  1);
                nlEnd(NL_ROW);
            }

            int opp = fec.opposite(c);
            if (opp>=0) { // symmetry of the curl correction
                for (int d : range(2)) {
                    nlRowScaling(100);
                    nlBegin(NL_ROW);
                    nlCoefficient(2*c+d, 1);
                    nlCoefficient(2*opp+0, Rij[c][d][0]);
                    nlCoefficient(2*opp+1, Rij[c][d][1]);
                    nlEnd(NL_ROW);
                }
            } else {
                vec3 edge = m.points[j] - m.points[i];
                vec2 n = {edge.y, -edge.x};

                mat<2,2> BiT = Bi[i].transpose();
                int tmp = (std::abs(BiT[0]*n) > std::abs(BiT[1]*n));

                nlRowScaling(100);
                nlBegin(NL_ROW);
                nlCoefficient(2*c+0,    tmp);
                nlCoefficient(2*c+1,  1-tmp);
                nlEnd(NL_ROW);
            }
        }

        for (int f : range(m.nfacets())) {
            if (ffsing[f]) continue;
            int hij = m.facet_corner(f, 0);
            int hjk = m.facet_corner(f, 1);
            int hki = m.facet_corner(f, 2);
            mat<2,2> Rij_ = Rij[hij];
            mat<2,2> Rik_ = Rij[hki].transpose();

            vec2 gij_ = gij[hij];
            vec2 gjk_ = gij[hjk];
            vec2 gki_ = gij[hki];
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

        for (int c : range(m.ncorners())) {
            cij[c] = vec2{nlGetVariable(2*c+0), nlGetVariable(2*c+1)};
            cij[c] = std::min(1., .35*gij[c].norm()/cij[c].norm())*cij[c];
        }
        nlDeleteContext(nlGetCurrent());
    }
*/

    CornerAttribute<vec2> uti(m);
    PointAttribute<vec2> ui(m.points);

    { // PGP
        double av_length = average_edge_length(m);
        std::cerr << "average edge length is " << av_length << std::endl;
        const double scale = 0.4/av_length;
        nlNewContext();
        nlSolverParameteri(NL_NB_VARIABLES, 4*m.nverts());
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlBegin(NL_SYSTEM);

        for (int c : range(m.ncorners())) {
            if (fec.opposite(c)>=0) continue;

            int i = fec.from(c), j = fec.to(c);
            vec3 edge = m.points[j] - m.points[i];
            vec2 n = {edge.y, -edge.x};

            if (std::abs(Bi[i].col(0)*n) > std::abs(Bi[i].col(1)*n)) { // TODO both constraints on corners
                nlSetVariable(4*i+0, 1); // u is integer
                nlSetVariable(4*i+1, 0);
                nlLockVariable(4*i+0);
                nlLockVariable(4*i+1);
            } else {
                nlSetVariable(4*i+2, 1); // v is integer
                nlSetVariable(4*i+3, 0);
                nlLockVariable(4*i+2);
                nlLockVariable(4*i+3);
            }
        }

        nlBegin(NL_MATRIX);

        for (int c : range(m.ncorners())) {
            int i = fec.from(c), j = fec.to(c);

            vec2 gij_ = scale*(gij[c] + cij[c]);
            mat<2,2> &Rij_ = Rij[c];

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

        for (int hi : range(m.ncorners())) // transfer tex coords from verts to corners
            uti[hi] = ui[fec.from(hi)];

        for (int i : range(3)) {
            std::cerr << "gij: " << scale*gij[m.facet_corner(1,i)] << std::endl << std::endl;
        }

        for (int i : range(3)) {
            std::cerr << "Rij: " << Rij[m.facet_corner(1,i)] << std::endl << std::endl;
        }

        for (int i : range(3)) {
            std::cerr << "Bi: " << Bi[m.vert(1,i)] << " ui: " << ui[m.vert(1,i)] << std::endl << std::endl;
        }

        for (int f : range(m.nfacets())) {
        if (f!=1) continue;
            int hij = m.facet_corner(f, 0);
            int hjk = m.facet_corner(f, 1);
            int hki = m.facet_corner(f, 2);
            mat<2,2> Rij_ = Rij[hij];
            mat<2,2> Rik_ = Rij[hij]*Rij[hjk];

            uti[hjk] = Rij_*uti[hjk];
            uti[hki] = Rik_*uti[hki];

              std::cerr << "uti: " <<uti[hij] <<  uti[hjk] << " " << uti[hki] << std::endl;

            vec2 gij_ = gij[hij];
            vec2 gik_ = gij[hij] + (Rij_*gij[hjk]);

              std::cerr << gij_ << " " << gik_ << std::endl;
              continue;
                for (int d : range(2)) {
                    while (uti[hij][d] + gij_[d] - uti[hjk][d] >  .5) uti[hjk][d] += 1;
                    while (uti[hij][d] + gij_[d] - uti[hjk][d] < -.5) uti[hjk][d] -= 1;
                    while (uti[hij][d] + gik_[d] - uti[hki][d] >  .5) uti[hki][d] += 1;
                    while (uti[hij][d] + gik_[d] - uti[hki][d] < -.5) uti[hki][d] -= 1;
                }

         /*
            mat<2,2> R = mat<2,2>::identity();
            for (int e : range(2)) {
                int ci = m.facet_corner(f, e);
                int cj = m.facet_corner(f, e+1);
                vec2 gij_ =  scale*(R*(gij[ci] + cij[ci]));
              R = Rij[ci]*R;
              uti[cj] = R*uti[cj];
//              std::cerr << R << "\n\n";


                for (int d : range(2)) {
                    while (uti[ci][d] + gij_[d] - uti[cj][d] >  .5) uti[cj][d] += 1;
                    while (uti[ci][d] + gij_[d] - uti[cj][d] < -.5) uti[cj][d] -= 1;
                }
            }
            */
        }

        nlDeleteContext(nlGetCurrent());
    }
    write_geogram("pgp.geogram", m, { {{"theta", theta.ptr}, {"debugBi", debugBi.ptr}, {"ui", ui.ptr}}, {{"ffsing", ffsing.ptr}}, {{"uti", uti.ptr},{"debugRij", debugRij.ptr},{"lci", lci.ptr}} });

    return 0;
}

