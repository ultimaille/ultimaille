#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>

#include <ultimaille/all.h>
#include <OpenNL_psm/OpenNL_psm.h>

#define UNUSED(x) (void)(x)

using namespace UM;

void compute_cross_field(Triangles &m, CornerAttribute<bool> &feature_edge, CornerAttribute<vec2> &uv, CornerAttribute<int> &layer_shift, int nb_FF_iters=1);
vec3 facet_normal(const Surface &m, const int f);


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

vec3 triangle_gradient(const Triangles &m, const int f, const double a, const double b, const double c) {
    vec3 ea = m.points[m.vert(f,2)] - m.points[m.vert(f,1)];
    vec3 eb = m.points[m.vert(f,0)] - m.points[m.vert(f,2)];
    vec3 ec = m.points[m.vert(f,1)] - m.points[m.vert(f,0)];
    vec3 n = cross(ec, -eb);
    vec3 g90 = -(a*ea + b*eb + c*ec)/n.norm();
    return cross(g90, n.normalize());
}

double triangle_area(const vec2 A, const vec2 B, const vec2 C) {
    const vec2 pts[3] = {A, B, C};
    double area = 0;
    for (int v : range(3)) {
        vec2 a = pts[v];
        vec2 b = pts[(v+1)%3];
        area += (b.y-a.y)*(b.x+a.x)/2;
    }
    return area;
}


int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }
    nlInitialize(argc, argv);

    // TODO assert manifoldity
    Triangles m;
#if 0
    read_wavefront_obj(argv[1], m);
//  for (vec3 &p : m.points) p.z = 0; // make sure it is 2D

    CornerAttribute<vec2> uv(m);
    CornerAttribute<bool> feature_edge(m);
    CornerAttribute<int> layer_shift(m);
    compute_cross_field(m, feature_edge, uv, layer_shift, 5);

    write_geogram("ff.geogram", m, {{}, {}, {{"uv", uv.ptr}, {"feature_edge", feature_edge.ptr}, {"layer_shift", layer_shift.ptr}}});

#else
    SurfaceAttributes attr = read_geogram(argv[1], m);

    CornerAttribute<vec2> uv("uv", attr, m);
    CornerAttribute<bool> feature_edge("feature_edge", attr, m);
    CornerAttribute<int> layer_shift("layer_shift", attr, m);
#endif

    SurfaceConnectivity fec(m);

    SignedPairwiseEquality param_vars(m.ncorners()*2);
    SignedPairwiseEquality curlc_vars(m.ncorners()*2);

    const bool perform_curl_correction = true;

    for (int h : corner_iter(m)) {
        int opp = fec.opposite(h);
        if (feature_edge[h]) { // zero feautre edge condition
            assert(uv[h].x==uv[fec.next(h)].x || uv[h].y==uv[fec.next(h)].y); // uv is supposed to be snapped
            bool v = (uv[h].y==uv[fec.next(h)].y);
            if (perform_curl_correction)
                curlc_vars.merge(v+2*h, v+2*h, false);
            param_vars.merge(v+2*h, v+2*h, false);
            param_vars.merge(v+2*fec.next(h), v+2*fec.next(h), false);
            if (opp<0) continue;
        }
        assert(opp>=0);
        int ngh = fec.next(opp);
        int rij = layer_shift[h];
        assert(rij>=0 && rij<4);
        // rij=0 -> ui=+uj, vi=+vj
        // rij=1 -> ui=+vj, vi=-uj
        // rij=2 -> ui=-uj, vi=-vj
        // rij=3 -> ui=-vj, vi=+uj
        const bool sign[8] = { true, true, true, false, false, false, false, true };
        param_vars.merge(h*2+0, ngh*2+  rij%2, sign[rij*2  ]);
        param_vars.merge(h*2+1, ngh*2+1-rij%2, sign[rij*2+1]);
        if (perform_curl_correction) {
            curlc_vars.merge(h*2+0, opp*2+  rij%2, !sign[rij*2  ]);
            curlc_vars.merge(h*2+1, opp*2+1-rij%2, !sign[rij*2+1]);
        }
    }

    double av_length = average_edge_length(m);
    const double scale = .3/av_length;

    CornerAttribute<vec2> gkl(m);
    { // N.B. this is an edge attribute, so only the resp. half-edges are filled in the 1st pass, and the non-resp are assigned in the 2nd pass
        for (int h : corner_iter(m)) {
            int f = fec.facet(h);
            vec3 gu = triangle_gradient(m, f, uv[m.corner(f, 0)].x,uv[m.corner(f, 1)].x,uv[m.corner(f, 2)].x).normalize();
            vec3 gv = triangle_gradient(m, f, uv[m.corner(f, 0)].y,uv[m.corner(f, 1)].y,uv[m.corner(f, 2)].y).normalize();
            vec3 n = facet_normal(m, f);
            mat<3,3> B = { { gu, gv, n } };

            int opp = fec.opposite(h);
            int i = fec.from(h), j = fec.to(h);
            vec3 gkl_ = B*fec.geom(h)*scale;
            vec2 tmp = {gkl_.x, gkl_.y};
            if (opp<0 || i<j)
                gkl[h] += tmp/(2-int(opp<0));
            else {
                for (int r : range(layer_shift[h])) {
                    UNUSED(r);
                    tmp = vec2{-tmp.y, tmp.x};
                }
                gkl[opp] += -tmp/2;
            }
        }
        for (int h : corner_iter(m)) { // 2nd pass: fill non-resp. half-edges
            int opp = fec.opposite(h);
            int i = fec.from(h), j = fec.to(h);
            if (opp<0 || i<j) continue;
            gkl[h] = -gkl[opp];
            for (int r : range(layer_shift[opp])) {
                UNUSED(r);
                gkl[h] = vec2{-gkl[h].y, gkl[h].x};
            }
        }
    }

    if (perform_curl_correction) {
        std::vector<int> ccvar, ccsgn;
        int nsets = curlc_vars.reduce(ccvar, ccsgn);
        nlNewContext();
        nlSolverParameteri(NL_NB_VARIABLES, nsets);
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlBegin(NL_SYSTEM);

        for (int v : range(curlc_vars.size())) { // zero boundary condition
            if (ccsgn[v]) continue;
            nlSetVariable(ccvar[v], 0);
            nlLockVariable(ccvar[v]);
        }

        nlBegin(NL_MATRIX);

        for (int v : range(curlc_vars.size())) { // min \sum_{ij} ||c_{ij}||^2
            nlRowScaling(.01);
            nlBegin(NL_ROW);
            nlCoefficient(ccvar[v],  1);
            nlEnd(NL_ROW);
        }

        for (int f : facet_iter(m)) { // actual correction
            int c0 = m.corner(f, 0);
            int c1 = m.corner(f, 1);
            int c2 = m.corner(f, 2);
            for (int d : range(2)) {
                nlRowScaling(100);
                nlBegin(NL_ROW);
                nlCoefficient(ccvar[d+2*c0], ccsgn[d+2*c0]);
                nlCoefficient(ccvar[d+2*c1], ccsgn[d+2*c1]);
                nlCoefficient(ccvar[d+2*c2], ccsgn[d+2*c2]);
                nlRightHandSide(-gkl[c0][d] - gkl[c1][d] - gkl[c2][d]);
                nlEnd(NL_ROW);
            }
        }

        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);
        nlSolve();

        CornerAttribute<vec2> ckl(m);
        for (int c : corner_iter(m)) {
            for (int d : range(2))
                ckl[c][d] = ccsgn[c*2+d]*nlGetVariable(ccvar[c*2+d]);
            ckl[c] = std::min(1., .27*gkl[c].norm()/(1e-4+ckl[c].norm()))*ckl[c]; // clamp the curl correction
            gkl[c] += ckl[c]; // apply the curl correction
        }
        nlDeleteContext(nlGetCurrent());
    }

    std::vector<int> redvar, redsgn;
    int nsets = param_vars.reduce(redvar, redsgn);

    FacetAttribute<bool> invalid(m);
    CornerAttribute<vec2> tex_coord(m);
    {
        CornerAttribute<vec2> frac_coord(m);
        for (int iter : range(5)) { // few passes to enforce non-zero field norm
            nlNewContext();
            nlSolverParameteri(NL_NB_VARIABLES, 2*nsets);
            nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
            nlBegin(NL_SYSTEM);

            for (int v : range(param_vars.size())) { // integer constraints for the singular vertices
                if (redsgn[v]) continue;
                int ind = redvar[v];
                nlSetVariable(ind*2+0, 1);
                nlSetVariable(ind*2+1, 0);
                nlLockVariable(ind*2+0);
                nlLockVariable(ind*2+1);
            }


            nlBegin(NL_MATRIX);
            if (iter) for (int c : corner_iter(m)) { // enforce non-zero field norm
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

            for (int c : corner_iter(m)) { // actual PGP energy
                int i = fec.from(c), j = fec.to(c);
                if (fec.opposite(c)>=0 && i>j) continue;

                for (int d : range(2)) {
                    int vari = redvar[d+2*c];
                    int sgni = redsgn[d+2*c];
                    int varj = redvar[d+2*fec.next(c)];
                    int sgnj = redsgn[d+2*fec.next(c)];

                    nlRowScaling(1);
                    nlBegin(NL_ROW);
                    nlCoefficient(vari*2+0,  cos(2*M_PI*gkl[c][d]));
                    nlCoefficient(vari*2+1, -sin(2*M_PI*gkl[c][d])*sgni);
                    nlCoefficient(varj*2+0, -1);
                    nlEnd(NL_ROW);

                    nlRowScaling(1);
                    nlBegin(NL_ROW);
                    nlCoefficient(vari*2+1, cos(2*M_PI*gkl[c][d])*sgni);
                    nlCoefficient(vari*2+0, sin(2*M_PI*gkl[c][d]));
                    nlCoefficient(varj*2+1, -sgnj);
                    nlEnd(NL_ROW);
                }
            }

            nlEnd(NL_MATRIX);
            nlEnd(NL_SYSTEM);
            nlSolve();

            for (int c : corner_iter(m)) { // project the solution to fractional coords
                for (int d : range(2)) {
                    int var = redvar[c*2+d];
                    int sgn = redsgn[c*2+d];
                    frac_coord[c][d] = atan2(sgn*nlGetVariable(var*2+1), nlGetVariable(var*2+0))/(2.*M_PI);
                }
            }

            for (int f : facet_iter(m)) { // recover the integer part
                invalid[f] = false;
                vec2 t[3] = {{0,0},{0,0},{0,0}}; // integer translation
                for (int lv : range(3)) {
                    int ci = m.corner(f, lv);
                    int cj = fec.next(ci);
                    for (int d : range(2)) {
                        while (frac_coord[cj][d] - frac_coord[ci][d] - gkl[ci][d] + t[lv][d] >  .5) t[lv][d] -= 1;
                        while (frac_coord[cj][d] - frac_coord[ci][d] - gkl[ci][d] + t[lv][d] < -.5) t[lv][d] += 1;
                    }
                }
                if ((t[0]+t[1]+t[2]).norm2()>1e-2) { // PGP-singular triangle
                    invalid[f] = true;
                    for (int lv : range(3))
                        tex_coord[m.corner(f, lv)] = {0,0};
                } else {
                    int c = m.corner(f, 0);
                    tex_coord[c] = frac_coord[c];
                    tex_coord[fec.next(c)] = frac_coord[fec.next(c)] + t[0];
                    tex_coord[fec.prev(c)] = frac_coord[fec.prev(c)] + t[0] + t[1];
                }
                invalid[f] = invalid[f] || (triangle_area(tex_coord[m.corner(f, 0)], tex_coord[m.corner(f, 1)], tex_coord[m.corner(f, 2)])<0); // TODO find a good positive constaint instead of 0
            }

            nlDeleteContext(nlGetCurrent());
        }
    }

    write_geogram("pgp.geogram", m, { {}, {{"invalid", invalid.ptr}}, {{"tex_coord", tex_coord.ptr},{"layer_shift", layer_shift.ptr}} });

    return 0;
}

