#include <iostream>
#include <ultimaille/all.h>
#include <chrono>
#include <fstream>

using namespace UM;


static const std::string poly_str =
R"(MeshVersionFormatted 2

Dimension
3

Vertices
5
 1         0        0 0
 0.309017  0.951057 0 1
-0.809017  0.587785 0 4
-0.809017 -0.587785 0 9
 0.309017 -0.951057 0 16

Triangles
1
1 2 3 0

Quadrilaterals
1
4 5 1 3 1

End
)";

static const std::string tet_str =
R"(MeshVersionFormatted 2

Dimension
3

Vertices
5
0. 0. 0. 0
1. 0. 0. 1
0. 1. 0. 4
0. 1. 1. 9
0. 0. 1. 16

Tetrahedra
2
1 2 3 4 0
1 3 4 5 1

End
)";


static const std::string pyramid_str =
R"(MeshVersionFormatted 2

Dimension
3

Vertices
5
0.  0.  0.  0
1.  0.  0.  1
1.  1.  0.  4
0.  1.  0.  9
0.5 0.5 0.5 16

Pyramids
1
1 2 3 4 5 42

End
)";

int main(int argc, char** argv) {

	// Triangles m;


	// m.create_facets(2);

	// m.vert(0, 0) = 0;
	// m.vert(0, 1) = 1;
	// m.vert(0, 2) = 2;
	// m.vert(1, 0) = 1;
	// m.vert(1, 1) = 3;
	// m.vert(1, 2) = 2;




	// write_by_extension(argv[1], m);


	// --------------------------------------------------------
	// ----------- BENCHMARK ON GET VERTEX ON VOLUME FACET ----
	
	// if (argc < 2) {
	// 	std::cerr << "Filename arg was expected." << std::endl;
	// 	exit(1);
	// }

	// std::string filename = argv[1];

	// Tetrahedra m;
	// read_by_extension(filename, m);
	// m.connect();
	// bool ok = true;
	// for (auto f : m.iter_facets()) {
	// 	ok &= f.__vertex(0) == f.vertex(0);
	// 	ok &= f.__vertex(1) == f.vertex(1);
	// 	ok &= f.__vertex(2) == f.vertex(2);

	// 	if (!ok)
	// 		break;
	// }

	// if (!ok)
	// 	std::cout << "fail." << std::endl;
	// else  
	// 	std::cout << "success." << std::endl;

	// // Bench
	// auto start = std::chrono::steady_clock::now();
	// for (int i = 0; i < 500; i++) {
	// 	for (auto f: m.iter_facets()) {
	// 		f.__vertex(0);
	// 	}
	// }
	// auto end = std::chrono::steady_clock::now();
	// std::cout << "time:" << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << std::endl;

	// start = std::chrono::steady_clock::now();
	// for (int i = 0; i < 500; i++) {
	// 	for (auto f: m.iter_facets()) {
	// 		f.vertex(0);
	// 	}
	// }
	// end = std::chrono::steady_clock::now();
	// std::cout << "time:" << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << std::endl;

	
    // Hexahedra m;
    // read_by_extension("/home/tex/Documents/Models/catorus_hex.geogram", m);
	// for (auto c : m.iter_cells()) {

	// 	// std::cout << "n corners: " << c.ncorners() << std::endl;
	// 	// auto c0 = c.corner(7);
	// 	// std::cout << c0 << std::endl;
	// 	std::cout << "cell: " << c << std::endl;
	// 	for (auto corner : c.iter_corners()) {
	// 		std::cout << corner << std::endl;
	// 	}
	// }

	// PolyLine p;
	// PolyLineAttributes attributes = read_by_extension("/home/tex/sandbox_5.486954/output/out.mesh", p);
	// EdgeAttribute<int> region_attr("region", attributes, p);
	// PointAttribute<int> vreg_attr("region", attributes, p);
	// write_by_extension("test.mesh", p, {{{"bob", vreg_attr.ptr}}, {{"bobi", region_attr.ptr}}});

	// Triangles m;
	// SurfaceAttributes attributes = read_by_extension("/home/tex/sandbox_5.486954/output/out.mesh", m);
	// PointAttribute<int> p_attr("region", attributes, m);
	// FacetAttribute<int> f_attr("region", attributes, m);
	// write_by_extension("test.mesh", m, {{{"region", p_attr.ptr}}, {{"region", f_attr.ptr}}, {}});


    // static const std::string filename[2] = { "ultimaille-test-polygons-in.mesh", "ultimaille-test-polygons-out.mesh" };
    // std::ofstream ofs(filename[0], std::ios::binary);
    // ofs << poly_str;
    // ofs.close();

    // Polyhedron m[2] = {};
    // for (int i : range(2)) {
    //     VolumeAttributes attrs = read_by_extension(filename[i], m[i]);
    //     PointAttribute<int> vint("region", attrs, m[i]);
    //     FacetAttribute<int> fint("region", attrs, m[i]);

	// 	std::cout << "----" << std::endl;
	// 	for (auto f : m[i].iter_facets())
	// 	{
	// 		std::cout << "f:" << f <<", fint: " << fint[f] << std::endl;
	// 	}

    //     if (!i)
    //         write_by_extension(filename[1], m[0], {{{"region", vint.ptr}}, {{"region", fint.ptr}}, {}});
    // }


	// std::cout << "end" << std::endl;

    static const std::string filename[2] = { "ultimaille-test-pyramids-in.mesh", "ultimaille-test-pyramids-out.mesh" };
    std::ofstream ofs(filename[0], std::ios::binary);
    ofs << pyramid_str;
    ofs.close();

    Pyramids m[2] = {};
    for (int i : range(2)) {
        VolumeAttributes attrs = read_by_extension(filename[i], m[i]);
        PointAttribute<int> vint("region", attrs, m[i]);
        CellAttribute<int> cint("region", attrs, m[i]);

        for (auto c : m[i].iter_cells()) {
            std::cout << "hello: " << cint[c] << std::endl;
        }

        if (!i)
            write_by_extension(filename[1], m[0], {{{"region", vint.ptr}}, {{"region", cint.ptr}}, {}, {}});
    }

	return 0;
}

