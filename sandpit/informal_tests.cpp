#include <iostream>
#include <ultimaille/all.h>
#include <chrono>

using namespace UM;


int main(int argc, char** argv) {

	if (argc < 2)
		exit(1);

	// std::cout << "hello " << std::endl;
    // Triangles m;

	// m.points.create_points(4);
	// m.points[0] = vec3(0,0,0);
	// m.points[1] = vec3(1,0,0);
	// m.points[2] = vec3(0.5,0.5,0);
	// m.points[3] = vec3(2,0.5,0);

	// m.create_facets(2);

	// m.vert(0, 0) = 0;
	// m.vert(0, 1) = 1;
	// m.vert(0, 2) = 2;
	// m.vert(1, 0) = 1;
	// m.vert(1, 1) = 3;
	// m.vert(1, 2) = 2;


    // write_by_extension(argv[1], m);

	std::string filename = argv[1];

	Tetrahedra m;
	read_by_extension(filename, m);
	m.connect();
	bool ok = true;
	for (auto f : m.iter_facets()) {
		ok &= f.__vertex(0) == f.vertex(0);
		ok &= f.__vertex(1) == f.vertex(1);
		ok &= f.__vertex(2) == f.vertex(2);

		if (!ok)
			break;
	}

	if (!ok)
		std::cout << "fail." << std::endl;
	else  
		std::cout << "success." << std::endl;

	// Bench
	auto start = std::chrono::steady_clock::now();
	for (int i = 0; i < 500; i++) {
		for (auto f: m.iter_facets()) {
			f.__vertex(0);
		}
	}
	auto end = std::chrono::steady_clock::now();
	std::cout << "time:" << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << std::endl;

	start = std::chrono::steady_clock::now();
	for (int i = 0; i < 500; i++) {
		for (auto f: m.iter_facets()) {
			f.vertex(0);
		}
	}
	end = std::chrono::steady_clock::now();
	std::cout << "time:" << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << std::endl;

    return 0;
}

