#include <iostream>
#include <limits>
#include <cstdlib>
#include "ultimaille/constraints.h"
#include "ultimaille/mesh_io.h"
#include "ultimaille/polyline.h"
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

#if 0
vec3 facet_centroid(const Surface &m, const int f) {
    vec3 ave(0, 0, 0);
    for (int i : range(m.facet_size(f)))
        ave += m.points[m.vert(f, i)];
    return ave / double(m.facet_size(f));
}
#endif

mat<2,2> R90_k(const int k) {
    assert(k>=0 && k<4);
    const mat<2,2> R90 = {{ {0,-1}, {1,0} }};
    mat<2,2> curr = mat<2,2>::identity();
    for (int i=0; i<k; i++)
        curr = R90*curr;
    return curr;
}

void compute_cross_field(const Triangles &m, const SurfaceConnectivity &fec, FacetAttribute<double> &theta, FacetAttribute<mat<2,2>> &Bi, CornerAttribute<int> &Rij) {
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

            if (iter) for (int f : range(m.nfacets())) { // enforce non-zero field norm
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
                theta[f] = atan2(nlGetVariable(2*f+1), nlGetVariable(2*f))/4 + M_PI/2 * (m.points[m.vert(f, 0)].x>5) + (rand()%4)*M_PI/2.;

            nlDeleteContext(nlGetCurrent());
        }
    }

    for (int f : range(m.nfacets()))
        Bi[f] = {{ {cos(theta[f]), cos(M_PI/2+theta[f])}, {sin(theta[f]), sin(M_PI/2+theta[f])} }};

    const mat<2,2> R90 = {{ {0,-1}, {1,0} }};
    for (int c : range(m.ncorners())) {
        int opp = fec.opposite(c);
        if (opp<0 || opp<c) continue;
        int f = fec.c2f[c], fopp = fec.c2f[opp];

        double best_norm = std::numeric_limits<double>::max();
        mat<2,2> curr = mat<2,2>::identity();
        for (int k=0; k<4; k++) {
            double norm = (Bi[f] - curr*Bi[fopp]).norm();
            if (norm<best_norm) {
                best_norm = norm;
                Rij[c] = k;
            }
            curr = R90*curr;
        }
    }

    for (int c : range(m.ncorners())) {
        int opp = fec.opposite(c);
        if (opp<0 || opp>c) continue;
        Rij[c] = (4-Rij[opp])%4;
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

    CornerAttribute<int> Rij(m);
    FacetAttribute<double> theta(m);
    FacetAttribute<mat<2,2>> Bi(m);

    std::cerr << "Computing fec...";
    SurfaceConnectivity fec(m);
    std::cerr << "ok\n";

    std::cerr << "Computing cross field...";
    compute_cross_field(m, fec, theta, Bi, Rij);
    std::cerr << "ok\n";

    int nvars = m.ncorners()*2;
    SignedPairwiseEquality speq(nvars);

    for (int c : range(m.ncorners())) {
        int rij = Rij[c];
        assert(rij>=0 && rij<4);
        int opp = fec.opposite(c);
        if (opp<0) continue;
        int ngh = fec.next(opp);
        // rij=0 -> ui=+uj, vi=+vj
        // rij=1 -> ui=+vj, vi=-uj
        // rij=2 -> ui=-uj, vi=-vj
        // rij=3 -> ui=-vj, vi=+uj
        const bool sign[8] = { true, true, true, false, false, false, false, true };
        speq.merge(c*2+0, ngh*2+  rij%2, sign[rij*2  ]);
        speq.merge(c*2+1, ngh*2+1-rij%2, sign[rij*2+1]);
    }

    std::vector<int> redvar, redsgn;
    int nsets = speq.reduce(redvar, redsgn);

#if 0
    { // drop var reduction debug
        std::vector<int> resp(nsets, -1);
        for (int i : range(redvar.size()))
            resp[redvar[i]] = speq.root(i);

        PolyLine pl;
        SegmentAttribute<int> pls(pl);
        pl.create_segments(m.ncorners()*2);
        for (int c : range(m.ncorners())) {
            int i = fec.from(c);
            vec3 p = m.points[i];
            vec3 b = facet_centroid(m, fec.c2f[c]);
            for (int d : range(2)) {
                pl.points.push_back(p + (b-p)/10*(d+1));
                pl.vert(c*2+d, 0) = c*2+d;
                int r = resp[redvar[c*2+d]];
                pl.vert(c*2+d, 1) = r;
                pls[c*2+d] = redsgn[c*2+d];
            }
        }
        write_geogram("reduction_test.geogram", pl, {{}, {{"sign", pls.ptr}}});
    }
#endif


    CornerAttribute<vec2> tex_coord(m);
    {
        CornerAttribute<vec2> frac_coord(m);
        double av_length = average_edge_length(m);
        const double scale = 0.3/av_length;

        CornerAttribute<vec2> gij(m);
        { // N.B. this is an edge attribute, so only the resp. half-edges are filled in the 1st pass, and the non-resp are assigned in the 2nd pass
            for (int c : range(m.ncorners())) {
                int opp = fec.opposite(c);
                int i = fec.from(c), j = fec.to(c);
                mat<2,2> &B = Bi[fec.c2f[c]];
                vec3 edge = m.points[j] - m.points[i];
                vec2 gij_ = B.transpose()*vec2(edge.x, edge.y)*scale;
                if (opp<0 || i<j)
                    gij[c] += gij_/(2-int(opp<0));
                else
                    gij[opp] += -(R90_k(Rij[c])*gij_)/2; // note Rij[c] and not Rij[opp]!
            }
            for (int c : range(m.ncorners())) { // 2nd pass: fill non-resp. half-edges
                int opp = fec.opposite(c);
                int i = fec.from(c), j = fec.to(c);
                if (opp<0 || i<j) continue;
                gij[c] = -(R90_k(Rij[opp])*gij[opp]);
            }
        }

        for (int iter : range(3)) { // few passes to enforce non-zero field norm
            nlNewContext();
            nlSolverParameteri(NL_NB_VARIABLES, 2*nsets);
            nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
            nlBegin(NL_SYSTEM);

            for (int v : range(nvars)) { // integer constraints for the singular vertices
                if (redsgn[v]) continue;
                int ind = redvar[v];
                nlSetVariable(ind*2+0, 1);
                nlSetVariable(ind*2+1, 0);
                nlLockVariable(ind*2+0);
                nlLockVariable(ind*2+1);
            }

            for (int c : range(m.ncorners())) { // integer constraints for the boundaries
                if (fec.opposite(c)>=0) continue;
                vec3 edge = m.points[fec.to(c)] - m.points[fec.from(c)];
                vec2 n = {edge.y, -edge.x};
                mat<2,2> &B = Bi[fec.c2f[c]];
                bool v = (std::abs(B.col(1)*n) > std::abs(B.col(0)*n));
                int vari = redvar[v+2*c];
                int varj = redvar[v+2*fec.next(c)];
                nlSetVariable(vari*2+0, 1); // no need to use the sign for integer constraints
                nlSetVariable(vari*2+1, 0);
                nlSetVariable(varj*2+0, 1);
                nlSetVariable(varj*2+1, 0);
                nlLockVariable(vari*2+0);
                nlLockVariable(vari*2+1);
                nlLockVariable(varj*2+0);
                nlLockVariable(varj*2+1);
            }

            nlBegin(NL_MATRIX);

            if (iter) for (int c : range(m.ncorners())) { // enforce non-zero field norm
                for (int d : range(2)) {
                    int var = redvar[c*2+d];
                    int sgn = redsgn[c*2+d];

                    nlRowScaling(.05*iter);
                    nlBegin(NL_ROW);
                    nlCoefficient(var*2+0, 1);
                    nlRightHandSide(cos(2*M_PI*frac_coord[c][d]));
                    nlEnd(NL_ROW);

                    nlRowScaling(.05*iter);
                    nlBegin(NL_ROW);
                    nlCoefficient(var*2+1, sgn);
                    nlRightHandSide(sin(2*M_PI*frac_coord[c][d]));
                    nlEnd(NL_ROW);
                }
            }

            for (int c : range(m.ncorners())) { // actual PGP energy
                int i = fec.from(c), j = fec.to(c);
                if (fec.opposite(c)>=0 && i>j) continue;

                for (int d : range(2)) {
                    int vari = redvar[d+2*c];
                    int sgni = redsgn[d+2*c];
                    int varj = redvar[d+2*fec.next(c)];
                    int sgnj = redsgn[d+2*fec.next(c)];

                    nlRowScaling(1);
                    nlBegin(NL_ROW);
                    nlCoefficient(vari*2+0,  cos(2*M_PI*gij[c][d]));
                    nlCoefficient(vari*2+1, -sin(2*M_PI*gij[c][d])*sgni);
                    nlCoefficient(varj*2+0, -1);
                    nlEnd(NL_ROW);

                    nlRowScaling(1);
                    nlBegin(NL_ROW);
                    nlCoefficient(vari*2+1, cos(2*M_PI*gij[c][d])*sgni);
                    nlCoefficient(vari*2+0, sin(2*M_PI*gij[c][d]));
                    nlCoefficient(varj*2+1, -sgnj);
                    nlEnd(NL_ROW);
                }
            }

            nlEnd(NL_MATRIX);
            nlEnd(NL_SYSTEM);
            nlSolve();

            for (int c : range(m.ncorners())) { // project the solution to fractional coords
                for (int d : range(2)) {
                    int var = redvar[c*2+d];
                    int sgn = redsgn[c*2+d];
                    frac_coord[c][d] = atan2(sgn*nlGetVariable(var*2+1), nlGetVariable(var*2+0))/(2.*M_PI);
                    tex_coord[c][d] = frac_coord[c][d];
                }
            }

            for (int f : range(m.nfacets())) { // recover the integer part
                for (int lv : range(2)) { // N.B. 2 and not 3: tex_coord[m.facet_corner(f,0)] is fixed
                    int ci = m.facet_corner(f, lv);
                    int cj = fec.next(ci);
                    for (int d : range(2)) {
                        while (tex_coord[cj][d] - tex_coord[ci][d] - gij[ci][d] >  .5) tex_coord[cj][d] -= 1;
                        while (tex_coord[cj][d] - tex_coord[ci][d] - gij[ci][d] < -.5) tex_coord[cj][d] += 1;
                    }
                }
            }

            nlDeleteContext(nlGetCurrent());
        }
    }

    write_geogram("pgp.geogram", m, { {}, {{"theta", theta.ptr}}, {{"tex_coord", tex_coord.ptr},{"Rij", Rij.ptr}} });

    return 0;
}

