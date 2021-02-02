#include <iostream>
#include <cstdlib>

#include <ultimaille/knn.h>
#include <ultimaille/mesh_io.h>
#include <ultimaille/attributes.h>
#include <ultimaille/surface.h>
#include <ultimaille/hboxes.h>
#include <ultimaille/range.h>

using namespace UM;

double average_edge_size(const Surface &m) {
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

// Original idea from Martin Roberts: http://extremelearning.com.au/isotropic-blue-noise-point-sets/
// Reference Java implementation from Tommy Ettinger: https://github.com/tommyettinger/sarong/blob/master/src/test/java/sarong/PermutationEtc.java

constexpr int grid_size = 32;
constexpr int permx[32] = {11,3,13,28,14,8,17,7,5,11,0,1,6,19,18,15,2,16,27,31,26,29,20,4,12,24,9,25,21,23,30,22};
constexpr int permy[32] = {20,31,25,9,15,23,18,27,30,29,21,17,28,14,24,11,26,16,4,20,5,7,12,10,3,0,1,8,22,13,2,6};
vec2 permuted_grid(const int i, const int j, const int n, const double size) {
    return vec2(i/(double)n - (n-1-permy[j])/(double)(n*n), j/(double)n - permx[i]/(double)(n*n))*size;
}

void sample_exterior(const double size, Triangles &rocker) {
    PointSet seeds;
    double ave_len = average_edge_size(rocker);
    { // uniformly populate the square [0,size]^2 with blueish pointset
        int repeat = floor(size*(1./grid_size)/ave_len)+1;
        for (int rj : range(repeat))
            for (int ri : range(repeat))
                for (int j : range(grid_size))
                    for (int i : range(grid_size)) {
                        vec2 p = permuted_grid(i, j, grid_size, 1) + vec2(rj, ri);
                        seeds.push_back(vec3(p.x, p.y, 0));
                    }

        { // normalize point cloud between [0,size]^2
            vec3 min, max;
            seeds.bbox(min, max);
            float maxside = std::max(max.x-min.x, max.y-min.y);
            for (vec3 &p : seeds)
                p = (p-min)*size/maxside;
        }
    }

    { // delete points inside the domain and replace them with the domain points
        std::vector<BBox3> inboxes(rocker.nfacets());
        for (int f : range(rocker.nfacets())) {
            for (int lv : range(3))
                inboxes[f].add(rocker.points[rocker.vert(f, lv)]);
            inboxes[f].dilate(ave_len/4);
        }
        HBoxes hb(inboxes);

        std::vector<int> primitives;
        for (vec3 const &p : seeds) {
            BBox3 b;
            b.add(p);
            hb.intersect(b, primitives);
            if (!primitives.size()) rocker.points.push_back(p);
        }
    }
}

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }

    Triangles tri;
    {
        Polygons pm;
        read_wavefront_obj(argv[1], pm);
        pm.extract_triangles(tri);
        assert(tri.nfacets() == pm.nfacets());
    }

    for (vec3 &p : *tri.points.data) p.z = 0; // make sure it is 2D

    const double size = 10.;
    { // normalize the rocker, place it well inside the [0,side]^2 square
        vec3 min, max;
        tri.points.bbox(min, max);
        float maxside = std::max(max.x-min.x, max.y-min.y);
        for (vec3 &p : tri.points) {
            p = ((p - (max+min)/2.)/maxside)*size/1.3 + vec3(1,1,0)*size/2;
        }
    }

    sample_exterior(size, tri);

//  write_geogram_ascii("drop.geogram_ascii", pm, { {"rand", pid.ptr} }, { {"id", fid.ptr} }, {{"id", cid.ptr}});

    /*
    std::vector<bool> to_kill(pm.nfacets(), false);
    for (int i=0; i<pm.nfacets(); i++)
        to_kill[i] = !(rand()%5);
    pm.delete_facets(to_kill);
    */

    write_geogram("drop2.geogram", tri, { {}, {}, {} });


    return 0;
}

