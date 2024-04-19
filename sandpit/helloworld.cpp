#include <iostream>
#include <ultimaille/all.h>

using namespace UM;

int main(int argc, char** argv) {

    if (argc < 2) {
        std::cerr << "Filename arg was expected." << std::endl;
        exit(1);
    }

    Triangles m;
    read_by_extension(argv[1], m);
    return 0;
}

