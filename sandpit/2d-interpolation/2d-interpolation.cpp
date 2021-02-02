#include <iostream>
#include <queue>
#include <cstdlib>
#include <ultimaille/mesh_io.h>
#include <ultimaille/surface.h>
#include <ultimaille/attributes.h>
#include <ultimaille/range.h>
#include <ultimaille/hboxes.h>

#include <OpenNL_psm/OpenNL_psm.h>

using namespace UM;

vec3 barycentric(const vec2 tri[3], const vec2 P) {
    mat<3,3> ABC = {{ {tri[0].x, tri[0].y, 1}, {tri[1].x, tri[1].y, 1}, {tri[2].x, tri[2].y, 1}}};
    if (ABC.det()<1e-5) return vec3(-1,1,1); // for a degenerate triangle generate negative coordinates, it will be thrown away by the rasterizator
    return ABC.invert_transpose() * vec3{P.x, P.y, 1};
}

double facet_area(Surface &m, int f) {
    double area = 0;
    for (int v=0; v<m.facet_size(f); v++) {
        vec3 a = m.points[m.vert(f, v)];
        vec3 b = m.points[m.vert(f, (v+1)%m.facet_size(f))];
        area += (b.y-a.y)*(b.x+a.x)/2;
    }
    return area;
}

int main(int argc, char** argv) {
    if (3>argc) {
        std::cerr << "Usage: " << argv[0] << " polyline.geogram mesh.geogram" << std::endl;
        return 1;
    }
    nlInitialize(argc, argv);

    PointSet pts;
    {
        Polygons pm;
        read_geogram(argv[1], pm);
        pts = {pm.points.data};
        for (vec3 &p : pts) p.z = 0; // make sure it is 2D
    }

    Triangles tri;
    {
        Polygons pm;
        read_geogram(argv[2], pm);
        pm.extract_triangles(tri);
        assert(tri.nfacets() == pm.nfacets());
        for (vec3 &p : tri.points) p.z = 0; // make sure it is 2D
    }

    PointAttribute<double> fct(tri.points);

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, tri.nverts());
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlBegin(NL_SYSTEM);
    nlBegin(NL_MATRIX);

    std::vector<BBox3> inboxes(tri.nfacets());
    for (int f : range(tri.nfacets()))
        for (int lv : range(3))
            inboxes[f].add(tri.points[tri.vert(f, lv)]);
    HBoxes hb(inboxes);

    std::vector<int> triid;
    for (vec3 &p : pts) {
        BBox3 b;
        b.add(p);
        hb.intersect(b, triid);
        for (int t : triid) {
            vec2 T[3];
            for (int lv : range(3))
                for (int d : range(2))
                    T[lv][d] = tri.points[tri.vert(t, lv)][d];
            vec3 b = barycentric(T, vec2{p.x, p.y});
            if (b.x<0 || b.y<0 || b.z<0) continue;
            nlRowScaling(100);
            nlBegin(NL_ROW);
            for (int lv : range(3)) {
                nlCoefficient(tri.vert(t, lv), b[lv]);
                std::cerr << tri.vert(t, lv) << " " << b[lv] << std::endl;
            }
            nlEnd(NL_ROW);
            break;
        }
    }

    const vec2 goal = {1, 0};
    for (int t : range(tri.nfacets())) {
        double A = facet_area(tri, t);
        for (int d : range(2)) {
            nlRowScaling(.1);
            nlBegin(NL_ROW);
            for (int lv : range(3)) {
                vec3 e = tri.points[tri.vert(t, (lv+1)%3)] - tri.points[tri.vert(t, lv)];
                vec2 n = {e.y, -e.x};
                nlCoefficient(tri.vert(t, (lv+2)%3), n[d]/(-2*A));
            }
            nlRightHandSide(goal[d]);
            nlEnd(NL_ROW);
        }
    }

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
    nlSolve();
    for (int v : range(tri.nverts())) {
        fct[v] = nlGetVariable(v);
    }

    nlDeleteContext(nlGetCurrent());

    write_geogram("fct.geogram", tri, { {{"fct", fct.ptr}}, {}, {} });
    return 0;
}

