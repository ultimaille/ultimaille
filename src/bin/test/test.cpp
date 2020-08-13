#include <iostream>
#include <cstdlib>
#include <cstdlib>
#include "ultimaille/knn.h"
#include "ultimaille/mesh_io.h"
#include "ultimaille/surface.h"

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }
    PolyMesh pm;
    read_wavefront_obj(argv[1], pm);

/*
    int npts = 50000000L;
    std::vector<vec3> pts(npts);
    for (vec3 &p : pts)
        for (int d=0; d<3; d++)
            p[d] = rand()/static_cast<double>(RAND_MAX);
    std::cerr << "making the query...";
    KNN knn(pts);
    std::vector<int> nn = knn.query(pts[13], 27);
    std::cerr << "ok!\n";

    for (int &i : nn)
        std::cout << pts[i] << std::endl;
*/
    write_wavefront_obj("drop.obj", pm);
    return 0;
}

