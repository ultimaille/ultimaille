#include <iostream>
#include <cassert>
#include <limits>
#include "range.h"
#include "voronoi.h"

ConvexCell::ConvexCell(const vec2 min, const vec2 max) {
    status = success;

    clip[0] = vec3( 1.0,  0.0, -min.x);
    clip[1] = vec3( 0.0,  1.0, -min.y);
    clip[2] = vec3(-1.0,  0.0,  max.x);
    clip[3] = vec3( 0.0, -1.0,  max.y);
    nb_v = 4;

    tr[0] = 1;
    tr[1] = 2;
    tr[2] = 3;
    tr[3] = 0;
    nb_t = 4;

    for (int t : range(4))
        vpos[t] = vertex(t);

    head = 0;
}

void ConvexCell::clip_by_line(vec3 eqn) {
    if (status!=success) return;

    int beg = -1; // find the portion of the polygon to cut
    int end = -1;
    int cur = head;
    do {
        bool a = vpos[   cur ]*eqn<0;
        bool b = vpos[tr[cur]]*eqn<0;
        if (a) nb_t--;

        if (!a &&  b && beg!=-1) {
            status = inconsistent_boundary;
            return;
        }

        if (!a &&  b) beg = cur;
        if ( a && !b) end = cur;
        cur = tr[cur];
    } while (cur!=head);

    if (nb_t<1) {
        status = empty_cell;
        return;
    }

    assert((beg<0 && end<0) || (beg>=0 && end>=0));
    if (beg<0) return;

    if (nb_v >= _MAX_P_) {
        status = vertex_overflow;
        return;
    }

    tr[nb_v] = tr[end]; // update the linked list
    head = tr[tr[beg]] = nb_v;
    clip[nb_v++] = eqn;

    nb_t += 2; // compute 2 new vertices
    vpos[head]    = vertex(head);
    vpos[tr[beg]] = vertex(tr[beg]);
}

vec3 ConvexCell::vertex(int t) {
    return cross(clip[t], clip[tr[t]]);
}

void ConvexCell::export_verts(std::vector<vec2> &verts) {
    verts = std::vector<vec2>(nb_t);
    int cur = head;
    int cnt = 0;
    do {
        vec3 &p = vpos[cur];
        verts[cnt++] = vec2(p.x/p.z, p.y/p.z);
        cur = tr[cur];
    } while (cur!=head);
}

/*
bool voronoi_cell(const KNN &knn, ConvexCell &cc, std::vector<vec2> &verts, const std::vector<vec2> &pts, const int seed, const int k) {
    size_t iter = 0;
    size_t k_ = k;
    while (++iter<pts.size()) {
        ConvexCell cc_ = cc;

        { // actual clipping
            std::vector<size_t> nn = knn.query(pts[seed], k_);
            for (size_t i=1; i<k_; i++) { // clip by half-planes
                vec2 B = pts[nn[i]];
                vec2 dir = pts[seed] - B;
                vec3 eq = vec3(dir.x, dir.y, -((pts[seed]+B)*dir)/2.);
                cc_.clip_by_line(eq);
            }
            if (cc_.status==success) { // export verts and check the security radius
                cc_.export_verts(verts);
                double farthest = (pts[nn.back()] - pts[seed])*(pts[nn.back()] - pts[seed]);
                double maxlen = -std::numeric_limits<double>::max();
                for (auto const &v : verts)
                    maxlen = std::max(maxlen, (v-pts[seed])*(v-pts[seed]));
                if (maxlen*4>=farthest)
                    cc_.status = security_radius_not_reached;
            }
        }

        if (cc_.status==success || (cc_.status==security_radius_not_reached && k_==pts.size())) {
            cc = cc_;
            cc.status = success;
            return true;
        }
        if (cc_.status!=security_radius_not_reached) {
            cc = cc_;
            std::cerr << "Clipping error: " << cc_.status << std::endl;
            return false;
        }
        k_ = std::min(k_*2, pts.size());
//      std::cerr << "K: " << k_ << std::endl;
    }
    std::cerr << "Clipping error, normally this part of code should not be reachable" << std::endl;
    return false;
}

/*
void voronoi_diagram(const std::vector<vec2> &seeds, Mesh &voronoi, double colocate_tolerance) {
    vec2 min, max;
    get_bbox(seeds, min, max);
    voronoi = Mesh();

    std::cerr << "Computing the Voronoi diagram...";
    std::vector<std::vector<vec2> > vorverts(seeds.size());

    KNN knn(seeds);
#pragma omp parallel for
    for (size_t i=0; i<seeds.size(); i++) {
        ConvexCell cc(min, max);
        bool res = voronoi_cell(knn, cc, vorverts[i], seeds, i, 16);
        assert(res);
    }
    for (auto const &verts : vorverts) {
        size_t nverts = verts.size();
        size_t off_v = voronoi.create_vertices(nverts);
        size_t off_f = voronoi.create_facets(1, nverts);
        for (size_t j=0; j<nverts; j++) {
            voronoi.point(off_v+j) = verts[j];
            voronoi.vert(off_f, j) = off_v+j;
        }
    }
    glue_verts(voronoi, colocate_tolerance);
    std::cerr << "ok" << std::endl;
}

*/
