#include <catch2/catch.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <ultimaille/all.h>

using namespace UM;


TEST_CASE("VTK + attributes IO test", "[VTK]") {
    std::string vtk_str;
    static const std::string filename[2] = { std::string(TEST_INPUT_DIR) + "test-with-attribute.vtk", "ultimaille-test-hexme-out.vtk" };

    Tetrahedra m[2] = {};
    VolumeAttributes attr[2] = {};
    for (int i : range(2)) {
        attr[i] = read_by_extension(filename[i], m[i]);
        REQUIRE( m[i].nverts()==263 );
        REQUIRE( m[i].ncells()==886 );
        if (i) {
            OppositeFacet conn(m[0]);
            PointAttribute<bool> b1("boundary", attr[1], m[1]);
            PointAttribute<bool> b2(m[1], false);
            for (int f : range(m[0].nfacets())) {
                if (conn[f]>=0) continue;
                for (int lv : range(3))
                    b2[m[0].facet_vert(f/4, f%4, lv)] = true;
            }
            for (int v : vert_iter(m[1]))
                CHECK ( b1[v]==b2[v] );
        } else {
           write_by_extension(filename[1], m[0], attr[0]);
        }
    }
    write_geogram("gna.geogram", m[0], attr[0]);
}

