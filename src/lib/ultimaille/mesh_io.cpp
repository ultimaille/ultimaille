#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "surface.h"
#include "mesh_io.h"

// supposes .obj file to have "f " entries without slashes
// TODO: improve the parser
void read_wavefront_obj(const std::string filename, PolyMesh &m) {
    m = PolyMesh();
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            vec3 v;
            for (int i=0;i<3;i++) iss >> v[i];
            m.points.data->push_back(v);
        } else if (!line.compare(0, 2, "f ")) {
            int idx, cnt=0;
            iss >> trash;
            while (iss >> idx) {
                m.facets.push_back(--idx);  // in wavefront obj all indices start at 1, not zero
                cnt++;
            }
            m.offset.push_back(m.offset.back()+cnt);
        }
    }
    std::cerr << "#v: " << m.nverts() << " #f: "  << m.nfacets() << std::endl;
}

void write_wavefront_obj(const std::string filename, const Surface &m) {
    std::fstream out;
    out.open(filename, std::ios_base::out);
    out << std::fixed << std::setprecision(4);
    for (int v=0; v<m.nverts(); v++)
        out << "v " << m.points[v] << std::endl;
    for (int f=0; f<m.nfacets(); f++) {
        out << "f ";
        for (int v=0; v<m.facet_size(f); v++)
            out << (m.vert(f,v)+1) << " ";
        out << std::endl;
    }
    out.close();
}

