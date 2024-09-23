#include <iostream>
#include <ultimaille/all.h>
#include <ctime>
#include <cstdlib>

#include <chrono>

using Clock = std::chrono::high_resolution_clock;
  using namespace std::literals::chrono_literals;

double rand_range(double LO, double HI) {
    return LO + static_cast<double>(rand())/(static_cast<double>(RAND_MAX/(HI-LO)));
}

using namespace UM;

int main(int argc, char** argv) {
    srand (static_cast <unsigned> (time(0)));

    if (argc < 2) {
        std::cerr << "Filename arg was expected." << std::endl;
        exit(1);
    }

    Triangles m;
    read_by_extension(argv[1], m);
//  m.connect();




    NearestPointOnMesh npos(m);

    int nb = 3;
    std::vector<vec3> test_pts;
    {
        auto box = m.points.util.bbox();
        box.dilate(.3 * (box.max - box.min).norm());
        for (int i : range(nb)) {
            vec3 R;
            for (int d : {0,1,2})
                R[d] = rand_range(box.min[d], box.max[d]);
            test_pts.push_back(R);
        }
    }

    {
        auto starting_time = Clock::now();
        PolyLine pl;
        pl.create_edges(nb);
        pl.points.create_points(2 * nb);
        for (int sample : range(nb)) {
            vec3 R = test_pts[sample];
            auto P = npos.query(R);
            std::cerr << R << std::endl << P.f << " " << P.p << std::endl << std::endl;


            for (int i : {0, 1})
                pl.vert(sample, i) = 2 * sample + i;
            pl.points[pl.vert(sample, 0)] = R;
            pl.points[pl.vert(sample, 1)] = P;
        }

        std::cerr << "Running time: " << (Clock::now()-starting_time)/1.s << " seconds" << std::endl;

        PointAttribute<bool> start(pl.points);
        for (auto v : pl.iter_vertices())
            start[v] = v%2;
        write_by_extension("test.geogram", pl, {{{"start", start.ptr}}, {}});
    }



    return 0;
}

