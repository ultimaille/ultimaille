#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "ultimaille/io/obj.h"

namespace UM {
    // supposes .obj file to have "f " entries without slashes
    // TODO: improve the parser
    // TODO: export vn and vt attributes
    void read_wavefront_obj(const std::string filename, Polygons &m) {
        m = Polygons();
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
                std::vector<int> f;
                int idx;
                iss >> trash;
                while (iss >> idx) {
                    f.push_back(--idx);  // in wavefront obj all indices start at 1, not zero
                }
                int off_f = m.create_facets(1, f.size());
                for (int i=0; i<static_cast<int>(f.size()); i++) {
                    m.vert(off_f, i) = f[i];
                }
            }
        }
        in.close();
        std::cerr << "#v: " << m.nverts() << " #f: "  << m.nfacets() << std::endl;
    }

    void read_wavefront_obj(const std::string filename, Triangles &m) {
        Polygons mpoly;
        read_wavefront_obj(filename, mpoly);

        std::vector<bool> to_kill(mpoly.nfacets(), false);
        for (int f=0; f<mpoly.nfacets(); f++)
            to_kill[f] = (3!=mpoly.facet_size(f));
        mpoly.delete_facets(to_kill);

        m.points = mpoly.points;
        m.facets = mpoly.facets;
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
}

