#include <iostream>
#include <queue>
#include <cstdlib>
#include <ultimaille/mesh_io.h>
#include <ultimaille/surface.h>
#include <ultimaille/attributes.h>
#include <ultimaille/range.h>

#include <OpenNL_psm/OpenNL_psm.h>
# define M_PI          3.141592653589793238462643383279502884L

using namespace UM;

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

bool compute_frame_field_wsing(const Triangles &m, const SurfaceConnectivity &fec, PointAttribute<double> &theta, PointAttribute<mat<2,2>> &Bi) {
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
                theta[vi] = atan2(nlGetVariable(2*vi+1), nlGetVariable(2*vi))/4;//  + M_PI/2 * (m.points[vi].x>5);// + (rand()%4)*M_PI/2.;

            nlDeleteContext(nlGetCurrent());
        }
    }

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

    for (int fi : range(m.nfacets())) {
        mat<2,2> R = mat<2,2>::identity();
        for (int lci : range(3))
            R = Rij[m.facet_corner(fi, lci)]*R;
        if (((R-mat<2,2>::identity()).norm()>1e-10)) return false;
    }

    std::vector<bool> visited(m.nverts());
    std::queue<int> q;
    q.push(0);
    visited[0] = true;

    while (!q.empty()) {
        int i = q.front();
        q.pop();
//        if (visited[i]) continue;

        int c = fec.v2c[i];
        do {
            int j = fec.to(c);
            if (!visited[j]) {
                visited[j] = true;
                q.push(j);
                while (theta[j] - theta[i]> M_PI/4) theta[j] -= M_PI/2;
                while (theta[j] - theta[i]<-M_PI/4) theta[j] += M_PI/2;
            }
            c = fec.c2c[c];
        } while (c != fec.v2c[i]);

    }
    for (int vi : range(m.nverts()))
        Bi[vi] = {{ {cos(theta[vi]), cos(M_PI/2+theta[vi])}, {sin(theta[vi]), sin(M_PI/2+theta[vi])} }};
    return true;
}


void compute_frame_field(const Triangles &m, const SurfaceConnectivity &fec, PointAttribute<double> &theta, PointAttribute<mat<2,2>> &Bi) {
    std::cerr << "Computing the ff constraints...";
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

            if      (angle> (M_PI*3)/4) angle -= M_PI;
            else if (angle<-(M_PI*3)/4) angle += M_PI;
            else if (angle>M_PI/4)      angle -= M_PI/2;
            else if (angle<-M_PI/4)     angle += M_PI/2;

            theta[vi] = angle;
        }
    }
    std::cerr << "ok\n";

    std::cerr << "Interpolating the ff...";
    { // interpolate the frame field
        nlNewContext();
        nlSolverParameteri(NL_NB_VARIABLES, m.nverts());
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlBegin(NL_SYSTEM);
        nlBegin(NL_MATRIX);

        for (int c : range(m.ncorners())) {
            if (fec.opposite(c)>=0 && fec.from(c)>fec.to(c)) continue;
            nlBegin(NL_ROW);
            nlCoefficient(fec.to(c)  ,  1);
            nlCoefficient(fec.from(c), -1);
            nlEnd(NL_ROW);
        }

        for (int v : range(m.nverts())) {
            if (!fec.is_border_vert(v)) continue;
            nlRowScaling(100);
            nlBegin(NL_ROW);
            nlCoefficient(v,  1);
            nlRightHandSide(theta[v]);
            nlEnd(NL_ROW);
        }

        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);
        nlSolve();

        for (int v : range(m.nverts()))
            theta[v] = nlGetVariable(v);

        nlDeleteContext(nlGetCurrent());
    }
    std::cerr << "ok\n";

    for (int vi : range(m.nverts()))
        Bi[vi] = {{ {cos(theta[vi]), cos(M_PI/2+theta[vi])}, {sin(theta[vi]), sin(M_PI/2+theta[vi])} }};
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
    PointAttribute<mat<2,2>> Bi(m.points);
    FacetAttribute<int> ffsing(m);

    std::cerr << "Computing fec...";
    SurfaceConnectivity fec(m);
    std::cerr << "ok\n";
//  PointAttribute<vec2> A(m.points);
//  for (int v : range(m.nverts())) {
//      A[v].x = 1 + (10-m.points[v].y);
//      A[v].y = 1 + (10-m.points[v].x);
//  }

    bool ret = compute_frame_field_wsing(m, fec, theta, Bi);
    if (!ret) 
    compute_frame_field(m, fec, theta, Bi);

    write_geogram("pgp.geogram", m, { {{"theta", theta.ptr}}, {}, {} });


    CornerAttribute<vec2> gij(m);
    for (int c : range(m.ncorners())) {
        int i=fec.from(c), j=fec.to(c);
        vec3 edge = m.points[j] - m.points[i];
        gij[c] = ((Bi[i].transpose() + Bi[j].transpose())/2)*vec2(edge.x, edge.y);
//          gij[c].x *= A[fec.from(c)].x;
//          gij[c].y *= A[fec.from(c)].y;
    }

    CornerAttribute<vec2> cij(m);

    if (1) { // compute curl correction
        std::cerr << "Computing the curl correction...";
        nlNewContext();
        nlSolverParameteri(NL_NB_VARIABLES, 2*m.ncorners());
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlBegin(NL_SYSTEM);
        nlBegin(NL_MATRIX);

        for (int c : range(m.ncorners())) {
            int i = fec.from(c), j = fec.to(c);
            mat<2,2> Rij_ = mat<2,2>::identity();

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
                    nlCoefficient(2*opp+0, Rij_[d][0]);
                    nlCoefficient(2*opp+1, Rij_[d][1]);
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
            mat<2,2> Rij_ = mat<2,2>::identity();
            mat<2,2> Rik_ = mat<2,2>::identity();

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
        std::cerr << "solving...";
        nlSolve();

        for (int c : range(m.ncorners())) {
            cij[c] = vec2{nlGetVariable(2*c+0), nlGetVariable(2*c+1)};
//            cij[c] = std::min(1., .35*gij[c].norm()/cij[c].norm())*cij[c];
        }
        nlDeleteContext(nlGetCurrent());
        std::cerr << "ok\n";
    }

    { // apply curl correction
        double av_length = average_edge_length(m);
        std::cerr << "average edge length is " << av_length << std::endl;
        const double scale = 0.3/av_length;

        for (int c : range(m.ncorners())) {
            gij[c] = scale*(gij[c]+cij[c]);
//          gij[c].x *= A[fec.from(c)].x;
//          gij[c].y *= A[fec.from(c)].y;
        }
    }

    CornerAttribute<vec2> uti(m);
    PointAttribute<vec2> ui(m.points);

    { // PGP
        for (int iter : range(3)) {

            nlNewContext();
            nlSolverParameteri(NL_NB_VARIABLES, 4*m.nverts());
            nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
            nlBegin(NL_SYSTEM);

            for (int c : range(m.ncorners())) {
                if (fec.opposite(c)>=0) continue;

                int i = fec.from(c), j = fec.to(c);
                vec3 edge = m.points[j] - m.points[i];
                vec2 n = {edge.y, -edge.x};

                int offset_i = std::abs(Bi[i].col(0)*n) > std::abs(Bi[i].col(1)*n) ? 4*i : 4*i+2;
                int offset_j = std::abs(Bi[j].col(0)*n) > std::abs(Bi[j].col(1)*n) ? 4*j : 4*j+2;

                nlSetVariable(offset_i+0, 1); // integer constraint
                nlSetVariable(offset_i+1, 0);
                nlLockVariable(offset_i+0);
                nlLockVariable(offset_i+1);

                nlSetVariable(offset_j+0, 1); // integer constraint
                nlSetVariable(offset_j+1, 0);
                nlLockVariable(offset_j+0);
                nlLockVariable(offset_j+1);
            }

            nlBegin(NL_MATRIX);

            if (iter) for (int vi : range(m.nverts())) {
                for (int d : range(2)) {
                    nlRowScaling(.05*iter);
                    nlBegin(NL_ROW);
                    nlCoefficient(4*vi+2*d + 0,  1);
                    nlRightHandSide(cos(2.*M_PI*ui[vi][d]));
                    nlEnd(NL_ROW);
                    nlRowScaling(.05*iter);
                    nlBegin(NL_ROW);
                    nlCoefficient(4*vi+2*d + 1,  1);
                    nlRightHandSide(sin(2.*M_PI*ui[vi][d]));
                    nlEnd(NL_ROW);
                }

            }


            for (int c : range(m.ncorners())) {
                int i = fec.from(c), j = fec.to(c);

                mat<2,2> Rij_ = mat<2,2>::identity();

                for (int d : range(2)) {
                    assert(Rij_[d][0] == 0 || Rij_[d][1]==0);
                    double rdij = (Rij_[d][0] ? 0 : 1);
                    double sdij = Rij_[d][0] + Rij_[d][1];

                    nlBegin(NL_ROW);
                    nlCoefficient(4*i+d*2+0,  cos(2.*M_PI*gij[c][d]));
                    nlCoefficient(4*i+d*2+1, -sin(2.*M_PI*gij[c][d]));
                    nlCoefficient(4*j+2*rdij, -1);
                    nlEnd(NL_ROW);

                    nlBegin(NL_ROW);
                    nlCoefficient(4*i+d*2+1, cos(2.*M_PI*gij[c][d]));
                    nlCoefficient(4*i+d*2+0, sin(2.*M_PI*gij[c][d]));
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

            for (int f : range(m.nfacets())) {
                for (int e : range(2)) {
                    int ci = m.facet_corner(f, e);
                    int cj = m.facet_corner(f, e+1);
                    for (int d : range(2)) {
                        while (uti[ci][d] + gij[ci][d] - uti[cj][d] >  .5) uti[cj][d] += 1;
                        while (uti[ci][d] + gij[ci][d] - uti[cj][d] < -.5) uti[cj][d] -= 1;
                    }
                }
            }

            nlDeleteContext(nlGetCurrent());
        }
    }
    write_geogram("pgp.geogram", m, { {{"theta", theta.ptr}, {"ui", ui.ptr}}, {}, {{"uti", uti.ptr}} });

    return 0;
}

