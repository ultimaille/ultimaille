#include <iostream>
#include <ultimaille/all.h>

using namespace UM;

int main(int argc, char** argv) {
    Triangles m;
    read_by_extension(argv[1], m);
    return 0;
}

