#include <catch2/catch.hpp>
#include <ultimaille/all.h>
#include <cmath>

using namespace UM;

TEST_CASE("query one", "[k-NN]") {
    std::vector<vec3> pts = {{1, 0, 3}, {3, 1, 2}, {0, 2, 0}, {5, 4, 3}};
    KNN<3> knn(pts);
    std::vector<int> neigh = knn.query(vec3(0,1,0), 1);
    REQUIRE( neigh.size()==1 );
    CHECK( neigh[0] == 2 );
}

TEST_CASE("query multiple", "[k-NN]") {
    std::vector<vec3> pts = {{1, 2, 3}, {2, 1, 1}, {3, 0, 3}, {1, 3, 1}, {2, 2, 3}, {3, 1, 4}, {1, 0, 4}, {2, 3, 2}, {1, 1, 1}};
    KNN<3> knn(pts);
    std::vector<int> neigh = knn.query(vec3(0,1,0), 3);
    REQUIRE( neigh.size()==3 );
    CHECK( neigh[0] == 8 );
    CHECK( neigh[1] == 1 );
    CHECK( neigh[2] == 3 );
}

TEST_CASE("query many", "[k-NN]") {
    std::vector<vec3> pts = {
        {-22.58,-14.61,-41.04}, {-80.33,0.73,-97.58}, {42.48,-31.26,-72.12}, {-10.11,-69.45,58.85}, {-87.03,-38.55,90.76}, {-40.01,51.85,53.32},
        {-9.88,39.4,88.52}, {98.56,48.55,-1.05}, {-43.11,-31.89,91.65}, {-49.37,-16.84,42.65}, {1.31,-97.97,6.65}, {-55.78,-58.06,46.77},
        {-89.13,-8.26,-12.09}, {54.78,-55.75,72.38}, {-84.68,26.04,7.57}, {67.34,-24.62,-41.85}, {17.49,-6.15,-93.3}, {-62.44,97.48,-63.13}, {24.18,66.16,-43.69},
        {-15.52,-3.69,65.03}, {-9.47,54.79,8.91}, {-66.12,33.01,-26.22}, {-2.22,18.17,77.2}, {68.41,43.73,-58.54}, {90.74,-86.75,-51.48}, {-49.57,29.36,-76.84},
        {99.09,-17.6,7.54}, {55.18,80.72,-30.87}, {72.26,52.04,-74.93}, {2.27,-51.31,68.33}, {-64.53,0.02,-19.26}, {-63.97,24.51,-24.38}, {-87.99,43.85,39.45},
        {-6.89,79.39,64.71}, {-82.75,9.67,98.92}, {70.12,23.89,-21.59}, {10.94,14.83,-38.91}, {-37.09,-69.49,48.53}, {-7.88,-9.15,93.14}, {92.73,-85.2,91.72},
        {-71.91,75.52,6.88}, {3.79,-63.4,-31.08}, {35.29,-21.58,29.21}, {7.67,-40.22,-45.23}, {-36.83,-3.33,-70.86}, {64.22,15.32,-7.52}, {30.96,8.19,-28.58},
        {25.27,-8.97,75.98}, {-64.83,84.4,-94.6}, {53.04,9.49,8.68}
    };

    std::vector<vec3> queries = {
        {63.1,42.52,3.92},
        {30.16,81.1,-80.29},
        {51.24,24.44,-45.93},
        {78.61,-5.45,53.99},
        {-17.01,-26.54,41.71}
    };

    std::vector<std::vector<int>> expected = {
        {45,35,49,7,27,46,23,18,26,36},
        {18,28,27,23,36,16,46,35,25,17},
        {23,35,46,45,36,28,18,15,49,27},
        {26,42,49,47,13,45,7,35,22,39},
        {19,9,29,3,37,11,42,38,8,47}
    };

    KNN<3> knn(pts);

    for (int i=0; i<5; i++) {
        std::vector<int> neigh = knn.query(queries[i], 10);
        REQUIRE( neigh.size()==10 );
        for (int j=0; j<10; j++) {
            CHECK( neigh[j] == expected[i][j] );
        }
    }
}

TEST_CASE("query .geogram", "[k-NN]") {
    PointSet cloud, request;
    read_geogram(std::string(TEST_INPUT_DIR) +   "knn-cloud.geogram", cloud);
    read_geogram(std::string(TEST_INPUT_DIR) + "knn-request.geogram", request);

    {
        std::vector<bool> to_kill(request.size(), false);
        for (int i : range(request.size())) {
            for (int d : {0,1,2}) {
                if (std::isfinite(request[i][d])) continue;
                to_kill[i] = true;
            }
        }
        std::vector<int> old2new;
        request.delete_points(to_kill, old2new);
        write_geogram("knn-request.geogram", request);
    }
    {
        std::vector<bool> to_kill(cloud.size(), false);
        for (int i : range(cloud.size())) {
            for (int d : {0,1,2}) {
                if (std::isfinite(cloud[i][d])) continue;
                to_kill[i] = true;
            }
        }
        std::vector<int> old2new;
        cloud.delete_points(to_kill, old2new);
        write_geogram("knn-cloud.geogram", cloud);
    }

    KNN<3> knn(*cloud.data);

    for (vec3 &r: request) {
        vec3 n = cloud[knn.query(r)[0]];
        double d1 = (r - n).norm2();
        double d2 = std::numeric_limits<double>::max();

        for (vec3 &c: cloud)
            d2 = std::min(d2, (r - c).norm2());

        CHECK( d1 <= d2 );
    }
}
