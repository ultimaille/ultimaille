#include <iostream>
#include <limits>
#include <cstdlib>
#include "ultimaille/disjointset.h"
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
                theta[f] = atan2(nlGetVariable(2*f+1), nlGetVariable(2*f))/4;// + M_PI/2 * (m.points[m.vert(f, 0)].x>5);// + (rand()%4)*M_PI/2.;

            for (int i : range(100)) {
                theta[rand()%m.nfacets()] += M_PI/2;
            }

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
    SignedPairwiseEquality dset(nvars);

    for (int c : range(m.ncorners())) {
        assert(Rij[c]>=0 && Rij[c]<4);
        int opp = fec.opposite(c);
        if (opp<0) continue;
        int ngh = fec.next(opp);
        const bool sign[8] = { true, true, true, false, false, false, false, true };
        int rij = Rij[c];
        dset.merge(c*2+0, ngh*2+0+rij%2, sign[rij*2  ]);
        dset.merge(c*2+1, ngh*2+1-rij%2, sign[rij*2+1]);
    }

    std::vector<int> redvar, redsign;
    int nsets = dset.reduce(redvar, redsign);
    std::cerr << nsets << " " << nvars << std::endl;

    CornerAttribute<vec2> ui(m);
    {
        double av_length = average_edge_length(m);
        const double scale = 0.3/av_length;

        nlNewContext();
        nlSolverParameteri(NL_NB_VARIABLES, 2*nsets);
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlBegin(NL_SYSTEM);

        for (int v : range(nvars)) { // integer constraints for the singular vertices
            if (redsign[v]) continue;
            int ind = redvar[v];
            nlSetVariable(ind*2+0, 1);
            nlSetVariable(ind*2+1, 0);
            nlLockVariable(ind*2+0);
            nlLockVariable(ind*2+1);
        }

        for (int c : range(m.ncorners())) { // integer constraints for the boundaries
            if (fec.opposite(c)>=0) continue;
            int i = fec.from(c), j = fec.to(c);
            vec3 edge = m.points[j] - m.points[i];
            vec2 n = {edge.y, -edge.x};
            mat<2,2> &B = Bi[fec.c2f[c]];
            bool v = (std::abs(B.col(1)*n) > std::abs(B.col(0)*n));
            int vari = v+2*c;
            int varj = v+2*fec.next(c);
            nlSetVariable(redvar[vari]*2+0, 1);
            nlSetVariable(redvar[vari]*2+1, 0);
            nlSetVariable(redvar[varj]*2+0, 1);
            nlSetVariable(redvar[varj]*2+1, 0);
            nlLockVariable(redvar[vari]*2+0);
            nlLockVariable(redvar[vari]*2+1);
            nlLockVariable(redvar[varj]*2+0);
            nlLockVariable(redvar[varj]*2+1);
        }

        nlBegin(NL_MATRIX);

#if 1
//      CornerAttribute<vec2> gij(m);
//      for (int c : range(m.ncorners())) {
//          mat<2,2> &B = Bi[fec.c2f[c]];
//          vec3 edge = m.points[fec.to(c)] - m.points[fec.from(c)];
//          gij[c] = B.transpose()*vec2(edge.x, edge.y)*scale;
//      }

        for (int c : range(m.ncorners())) { // TODO
//            if (fec.opposite(c)>=0 && i>j) continue;
            mat<2,2> &B = Bi[fec.c2f[c]];
            int i = fec.from(c), j = fec.to(c);
            vec3 edge = m.points[j] - m.points[i];
            vec2 gij = B.transpose()*vec2(edge.x, edge.y)*scale;
            for (int d : range(2)) {
                int indi = redvar[c*2+d];
                int indj = redvar[fec.next(c)*2+d];
                nlBegin(NL_ROW);
                nlCoefficient(indi*2+0,  cos(2*M_PI*gij[d]));
                nlCoefficient(indi*2+1, -sin(2*M_PI*gij[d])*redsign[c*2+d]);
                nlCoefficient(indj*2+0, -1);
                nlEnd(NL_ROW);

                nlBegin(NL_ROW);
                nlCoefficient(indi*2+1, cos(2*M_PI*gij[d])*redsign[c*2+d]);
                nlCoefficient(indi*2+0, sin(2*M_PI*gij[d]));
                nlCoefficient(indj*2+1, -1*redsign[fec.next(c)*2+d]);
                nlEnd(NL_ROW);
            }
        }
#else
        for (int v : range(nvars)) {
             int ind = redvar[v];
             nlRowScaling(.01);
             nlBegin(NL_ROW);
             nlCoefficient(ind*2+0, 1);
             nlRightHandSide(.134234);
             nlEnd(NL_ROW);
             nlBegin(NL_ROW);
             nlCoefficient(ind*2+1, 1);
             nlRightHandSide(.134234);
             nlEnd(NL_ROW);
         }
#endif
        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);
        nlSolve();

        for (int c : range(m.ncorners())) {
            for (int d : range(2)) {
                int v = redvar[c*2+d];
                ui[c][d] = atan2(redsign[c*2+d]*nlGetVariable(v*2+1), nlGetVariable(v*2+0))/(2.*M_PI);
//                ui[c][d] = atan2(nlGetVariable(v*2+1), nlGetVariable(v*2+0))/(2.*M_PI);
            }
        }

        nlDeleteContext(nlGetCurrent());
    }


    CornerAttribute<int> uvarid(m), vvarid(m), uzero(m),vzero(m);

    for (int c : range(m.ncorners())) {
        uzero[c] = dset.is_zero(c*2);
        vzero[c] = dset.is_zero(c*2+1);
        uvarid[c] = redvar[c*2];
        vvarid[c] = redvar[c*2+1];
    }


    write_geogram("pgp.geogram", m, { {}, {{"theta", theta.ptr}}, {{"tex_coord", ui.ptr},{"Rij", Rij.ptr},{"uvarid", uvarid.ptr},{"vvarid", vvarid.ptr}, {"uzero", uzero.ptr}, {"vzero", vzero.ptr}} });
    return 0;
}

