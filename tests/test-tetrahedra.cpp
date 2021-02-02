//#define CATCH_CONFIG_RUNNER
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <iostream>
#include <ultimaille/all.h>

using namespace UM;

TEST_CASE( "Tetrahedra::Util::cell_volume" ) {
    Tetrahedra m;
    *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};
    m.cells = {0,1,2,3};
    REQUIRE( std::abs(m.util.cell_volume(0)-1./6.)<1e-13 );
}

/*
int main( int argc, char* argv[] ) {
    int result = Catch::Session().run( argc, argv );

    return result;
}


/*

int main() {
    {
        Tetrahedra m;
        *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};
        m.cells = {0,1,2,3};

        double vol = m.util.cell_volume(0);

        std::cerr << "Tet volume: " << vol << std::endl;
        assert(std::abs(vol-1./6.)<1e-14);

        for (int lf : range(4)) {
            vec3 n = m.util.facet_normal(0, lf);
            std::cerr << "Normal: " << n << std::endl;
            //          assert(n*(m.util.bary_facet(0, lf) - m.util.bary_verts(0))>0);
        }
        write_geogram("tet.geogram", m, {{}, {}, {}, {}});
    }

    {
        Hexahedra m;
        *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1} };
        m.cells = {0,1,2,3,4,5,6,7};

        double vol = m.util.cell_volume(0);
        std::cerr << "Hex volume: " << vol << std::endl;
        assert(std::abs(vol-1.)<1e-14);
        write_geogram("hex.geogram", m, {{}, {}, {}, {}});
    }

    return 0;
}

*/
