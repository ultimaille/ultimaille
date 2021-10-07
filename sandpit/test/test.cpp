#include <iostream>

#include <ultimaille/all.h>

using namespace UM;

int main(int argc, char** argv) {
    Polygons test;
    write_by_extension("gloglo.geogram", test);
    return 0;
}

