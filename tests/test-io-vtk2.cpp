#include <catch2/catch.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <ultimaille/all.h>

using namespace UM;


TEST_CASE("VTK + attributes IO test", "[VTK]") {
	std::string vtk_str;
    static const std::string filename[2] = { "ultimaille-test-hexme-in.vtk", "ultimaille-test-hexme-out.vtk" };
   // std::ifstream ofs(filename[0], std::ios::binary);
	std::cerr << TEST_INPUT_DIR << std::endl;
	
	//std::ofstream ofs(filename[0], std::ios::binary);
    //ofs << vtk_str;
    //ofs.close();

    /*
    Tetrahedra m[2] = {};
    for (int i : range(2)) {
        read_by_extension(filename[i], m[i]);

        REQUIRE( m[i].nverts()==4 );
        REQUIRE( m[i].ncells()==1 );
        CHECK( std::abs(m[i].util.cell_volume(0)-1./6.)<ftol );
        if (!i)
            write_by_extension(filename[1], m[0]);
    }
    */
}
