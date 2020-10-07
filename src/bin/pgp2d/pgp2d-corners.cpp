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
    DisjointSetWithSign dset(nvars);

    for (int c : range(m.ncorners())) {
        assert(Rij[c]>=0 && Rij[c]<4);
        int opp = fec.opposite(c);
        if (opp<0) continue;
        int ngh = fec.next(opp);
        const bool sign[8] = { true, true, false, true, false, false, true, false };
        int rij = Rij[c];
//      if (rij!=0) continue;
        dset.merge(c*2,   ngh*2+  rij%2, sign[rij*2  ]);
        dset.merge(c*2+1, ngh*2+1-rij%2, sign[rij*2+1]);
    }

    std::vector<int> indices;
    int nsets = dset.get_sets_id(indices);
    std::cerr << nsets << " " << nvars << std::endl;
    CornerAttribute<int> uvarid(m), vvarid(m);
//  for (int z : range(nvars)) {
//      std::cerr << indices[z] << std::endl;
//  }
    for (int c : range(m.ncorners())) {
        uvarid[c] = indices[c*2];
        vvarid[c] = indices[c*2+1];
    }


    write_geogram("pgp.geogram", m, { {}, {{"theta", theta.ptr}}, {{"Rij", Rij.ptr},{"uvarid", uvarid.ptr},{"vvarid", vvarid.ptr}} });
    return 0;
}

