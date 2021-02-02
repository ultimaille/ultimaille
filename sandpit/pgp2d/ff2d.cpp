#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>

#include <ultimaille/all.h>
#include <OpenNL_psm/OpenNL_psm.h>

using namespace UM;

#define DEBUG_DROP 1

// assume signed i, positive n. The result is in [0..n-1]
int modulo(const int i, const int n) {
    assert(n>0);
    return (i%n + n)%n; // Based on the C99 specification: i == (i / n) * n + i % n
}

double average_edge_length(const Surface &m);
vec3 triangle_gradient(const Triangles &m, const int f, const double a, const double b, const double c);

vec3 rotate_vector_around_axis(const vec3 axis, const double angle, const vec3 v) {
    vec3 u = vec3(axis).normalize();
    double c = cos(angle);
    double s = sin(angle);
    return vec3((c+u.x*u.x*(1-c))*v.x     + (u.x*u.y*(1-c)-u.z*s)*v.y + (u.x*u.z*(1-c)+u.y*s)*v.z,
                (u.x*u.y*(1-c)+u.z*s)*v.x + (c+u.y*u.y*(1-c))*v.y     + (u.y*u.z*(1-c)-u.x*s)*v.z,
                (u.x*u.z*(1-c)-u.y*s)*v.x + (u.z*u.y*(1-c)+u.x*s)*v.y + (c+u.z*u.z*(1-c))*v.z     );
}


vec3 facet_bary(const Surface &m, const int f) {
    vec3 ave(0, 0, 0);
    for (int v : facet_vert_iter(m, f))
        ave += m.points[v];
    return ave / double(m.facet_size(f));
}

vec3 facet_normal(const Surface &m, const int f) {
    vec3 bary = facet_bary(m, f);
    vec3 res(0, 0, 0);
    for (int fv : range(m.facet_size(f)))
        res = res + cross(m.points[m.vert(f, fv)] - bary, m.points[m.vert(f, (fv+1)%m.facet_size(f))] - bary);
    return res.normalize();
}

double vector_angle(const vec3 &v0, const vec3 &v1) {
    return atan2(cross(v0, v1).norm(), v0*v1);
}

double edge_angle_in_ref(const Triangles &m, const SurfaceConnectivity &fec, const int h) {
    int h_ref = m.corner(fec.facet(h), 0);
    if (h == h_ref) return 0;
    double angle = vector_angle(fec.geom(h), fec.geom(h_ref));
    if (h == m.corner(fec.facet(h), 1)) return angle;
    if (h == m.corner(fec.facet(h), 2)) return 2. * M_PI - angle;
    assert(false);
}

double c_ij(const Triangles &m, const SurfaceConnectivity &fec, int h) {
    assert(fec.opposite(h)>=0);
    if (fec.opposite(h) > h) return -c_ij(m, fec, fec.opposite(h));
    return M_PI - edge_angle_in_ref(m, fec, fec.opposite(h)) + edge_angle_in_ref(m, fec, h);
}

void compute_cross_field(const Triangles &m, const SurfaceConnectivity &fec, const CornerAttribute<bool> &feature_edge, FacetAttribute<double> &alpha, CornerAttribute<int> &layer_shift, int nb_FF_iters=1) {
    for (int h : corner_iter(m)) // compute the frame field at the feature edges
        if (feature_edge[h])
            alpha[fec.facet(h)] = edge_angle_in_ref(m, fec, h);

    for (int iter : range(nb_FF_iters)) { // interpolate the frame field
        auto context = nlNewContext();
        nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
        nlSolverParameteri(NL_NB_VARIABLES, 2*m.nfacets());
        nlBegin(NL_SYSTEM);

        for (int h : corner_iter(m)) { // lock the facet incident to feature edges
            if (!feature_edge[h]) continue;
            int f = fec.facet(h);
            vec2 goal = { cos(4*alpha[f]), sin(4*alpha[f]) };
            for (int d : range(2)) {
                nlSetVariable (2*f+d, goal[d]);
                nlLockVariable(2*f+d);
            }
        }

        nlBegin(NL_MATRIX);
        for (int h : corner_iter(m)) { // smoothing term
            if (feature_edge[h]) continue;
            int opp = fec.opposite(h); assert(opp>=0);

            double angle = 4*c_ij(m, fec, h);
            double rot[2][2] = {{cos(angle),-sin(angle)},{sin(angle),cos(angle)}};
            for (int d : range(2)) {
                nlBegin(NL_ROW);
                nlCoefficient(2*fec.facet(h) + d, -1);
                for (int dd : range(2))
                    nlCoefficient(fec.facet(opp)*2 + dd, rot[d][dd]);
                nlEnd(NL_ROW);
            }
        }

        if (iter) for (int f : facet_iter(m)) { // enforce non-zero field norm
            vec2 goal = { cos(4*alpha[f]), sin(4*alpha[f]) };
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

        for (int f : facet_iter(m))
            alpha[f] = atan2(nlGetVariable(2*f+1), nlGetVariable(2*f))/4;// + M_PI/2 * (m.points[m.vert(f, 0)].x>5) + (rand()%4)*M_PI/2.;
        nlDeleteContext(context);
    }

    for (int h : corner_iter(m)) { // compute the layer shift
        layer_shift[h] = 0;
        if (fec.opposite(h)<0) continue;
        int f   = fec.facet(h);
        int opp = fec.facet(fec.opposite(h));

        double noLSrot = alpha[opp] - alpha[f] /*- omega[h]*/ + c_ij(m, fec, h);
        double rot = noLSrot;
        while (rot >  M_PI_2/2) rot -= M_PI_2;
        while (rot < -M_PI_2/2) rot += M_PI_2;

        layer_shift[h] = -int(floor(0.5 + (noLSrot - rot) / M_PI_2));
        layer_shift[fec.opposite(h)] = -layer_shift[h];
    }

    for (int h : corner_iter(m)) { // normalize the layer shift
        layer_shift[h] = modulo(layer_shift[h], 4);
        assert(layer_shift[h]>=0 && layer_shift[h]<4);
    }
}

void local_basis(const SurfaceConnectivity &fec, const int h, vec3& x, vec3& y, vec3& z) {
    vec3 n = cross(fec.geom(fec.prev(h)), fec.geom(h));
    assert(n.norm2() > 1e-20);
    z = n.normalize();
    x = fec.geom(h).normalize();
    y = cross(z, x).normalize();
}

void export_local_uv(const Triangles &m, const SurfaceConnectivity &fec, const FacetAttribute<double> &alpha, const FacetAttribute<double> &beta, CornerAttribute<vec2> &uv) {
    for (int f : facet_iter(m)) {
        vec3 x, y, z;
        local_basis(fec, m.corner(f, 0), x, y, z);
        double a = alpha[f], b = beta[f];
        vec3 grads[2] = {cos(a)*x + sin(a)*y, cos(b)*x + sin(b)*y};
        for (int i : range(3)) {
            vec3 P = m.points[m.vert(f, i)];
            for (int d : range(2))
                uv[m.corner(f, i)][d] = grads[d]*P;
        }

        for (int i : range(3)) // snap tex coords
            for (int j : range(3))
                for (int d : range(2))
                    if (std::abs(uv[m.corner(f, i)][d] - uv[m.corner(f, j)][d]) < 1e-10)
                        uv[m.corner(f, i)][d] = uv[m.corner(f, j)][d];
    }
}

void init_feature_edge(const Triangles &m, const SurfaceConnectivity &fec, CornerAttribute<bool> &feature_edge) {
    for (int h : corner_iter(m)) {
        feature_edge[h] = false;
        if (fec.opposite(h)<0) {
            feature_edge[h] = true;
            continue;
        }
        vec3 n1 = facet_normal(m, fec.facet(h));
        vec3 n2 = facet_normal(m, fec.facet(fec.opposite(h)));
        if (n1*n2 < cos(M_PI / 5.))
            feature_edge[h] = true;
    }
}

void split_overconstrained_facets(Triangles &m, CornerAttribute<bool> &feature_edge) {
    std::vector<bool> to_kill(m.nfacets(), false);

    for (int f : facet_iter(m)) {
        int nb_features = 0;
        for (int i : range(3))
            nb_features += feature_edge[m.corner(f, i)];
        if (nb_features < 2) continue;

        int nv = m.points.push_back(facet_bary(m, f));
        int offset = m.create_facets(3);
        for (int i : range(3)) {
            to_kill.push_back(false);
            m.vert(offset + i, 0) = m.vert(f, i);
            m.vert(offset + i, 1) = m.vert(f, (i + 1) % 3);
            m.vert(offset + i, 2) = nv;
            feature_edge[m.corner(offset + i, 0)] = feature_edge[m.corner(f, i)];
            for (int c : range(2))
                feature_edge[m.corner(offset+i, 1+c)] = false;
        }
        to_kill[f] = true;
    }

    m.delete_facets(to_kill);
}


void compute_cross_field(Triangles &m, CornerAttribute<bool> &feature_edge, CornerAttribute<vec2> &uv, CornerAttribute<int> &layer_shift, int nb_FF_iters=1) {
    SurfaceConnectivity fec(m);
    init_feature_edge(m, fec, feature_edge);
    split_overconstrained_facets(m, feature_edge);
    fec.reset();


    std::cerr << "Computing cross field...";
    FacetAttribute<double> alpha(m), beta(m);
    compute_cross_field(m, fec, feature_edge, alpha, layer_shift, nb_FF_iters);
    for (int f : facet_iter(m)) beta[f] = alpha[f] + M_PI_2;
    std::cerr << "ok\n";

    export_local_uv(m, fec, alpha, beta, uv);

#if DEBUG_DROP
    { // drop var reduction debug
        PolyLine pl;
        pl.create_segments(m.nfacets()*2);
        double avlen = average_edge_length(m);
        for (int f : facet_iter(m)) {
            vec3 g = facet_bary(m, f);
            vec3 e1 = fec.geom(m.corner(f, 0));
            vec3 e2 = fec.geom(m.corner(f, 2));

            vec3 a = rotate_vector_around_axis(cross(e1, e2), -alpha[f], e1);
            vec3 b = rotate_vector_around_axis(cross(e1, e2), -beta[f], e1);

            int off = pl.nverts();
            pl.points.push_back(g);
            pl.points.push_back(g+a.normalize()*avlen*.2);
            pl.points.push_back(g+b.normalize()*avlen*.2);


            pl.vert(f*2+0, 0) = off;
            pl.vert(f*2+0, 1) = off+1;
            pl.vert(f*2+1, 0) = off;
            pl.vert(f*2+1, 1) = off+2;
        }
        write_geogram("field.geogram", pl, {{}, {}});
    }
    { // drop var reduction debug
        PolyLine pl;
        pl.create_segments(m.nfacets()*2);
        double avlen = average_edge_length(m);
        for (int f : facet_iter(m)) {
            vec3 g = facet_bary(m, f);

            vec3 gu = triangle_gradient(m, f, uv[m.corner(f, 0)].x,uv[m.corner(f, 1)].x,uv[m.corner(f, 2)].x);
            vec3 gv = triangle_gradient(m, f, uv[m.corner(f, 0)].y,uv[m.corner(f, 1)].y,uv[m.corner(f, 2)].y);

            int off = pl.nverts();
            pl.points.push_back(g);
            pl.points.push_back(g+gu.normalize()*avlen*.2);
            pl.points.push_back(g+gv.normalize()*avlen*.2);


            pl.vert(f*2+0, 0) = off;
            pl.vert(f*2+0, 1) = off+1;
            pl.vert(f*2+1, 0) = off;
            pl.vert(f*2+1, 1) = off+2;
        }
        write_geogram("field2.geogram", pl, {{}, {}});
    }
#endif
}


/*
int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }
    nlInitialize(argc, argv);

    Triangles m;
//    read_wavefront_obj(argv[1], m);
    read_geogram(argv[1], m);
    // TODO assert manifoldity

    CornerAttribute<vec2> uv(m);
    CornerAttribute<bool> feature_edge(m);
    CornerAttribute<int> layer_shift(m);
    compute_cross_field(m, feature_edge, uv, layer_shift, 5);


    write_geogram("pgp.geogram", m, {{}, {}, {{"uv", uv.ptr}, {"feature_edge", feature_edge.ptr}, {"layer_shift", layer_shift.ptr}}});

    return 0;

    return 0;
}
*/

