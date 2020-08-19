#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "surface.h"
#include "attributes.h"
#include "mesh_io.h"

// supposes .obj file to have "f " entries without slashes
// TODO: improve the parser
// TODO: export vn and vt attributes
void read_wavefront_obj(const std::string filename, PolyMesh &m) {
//    m = PolyMesh();
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

// Attention: only double attributes
typedef unsigned int index_t;
void write_geogram_ascii(const std::string filename, const Surface &m,
        std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > pattr,
        std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > fattr,
        std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > cattr) {
    std::fstream file;
    file.open(filename, std::ios_base::out);
    file << "[HEAD]\n\"GEOGRAM\"\n\"1.0\"\n[ATTS]\n\"GEO::Mesh::vertices\"\n";
    file << m.nverts() << "\n";
    file << "[ATTR]\n\"GEO::Mesh::vertices\"\n\"point\"\n\"double\"\n8 # this is the size of an element (in bytes)\n3 # this is the number of elements per item\n";
    file << std::fixed << std::setprecision(4);
    for (const vec3 &p : m.points)
        file << p.x << "\n" << p.y << "\n" << p.z << "\n";
    file << "[ATTS]\n\"GEO::Mesh::facets\"\n";
    file << m.nfacets() << "\n";
    file << "[ATTR]\n\"GEO::Mesh::facets\"\n\"GEO::Mesh::facets::facet_ptr\"\n\"index_t\"\n" << sizeof(index_t) << "\n1\n";
    for (int f=0; f<m.nfacets(); f++)
        file << m.facet_corner(f,0) << "\n";

    file << "[ATTS]\n\"GEO::Mesh::facet_corners\"\n";
    file << m.ncorners() << "\n";
    file << "[ATTR]\n\"GEO::Mesh::facet_corners\"\n\"GEO::Mesh::facet_corners::corner_vertex\"\n\"index_t\"\n" << sizeof(index_t) << "\n1\n";
    for (int f=0; f<m.nfacets(); f++)
        for (int v=0; v<m.facet_size(f); v++)
           file << m.vert(f,v) << "\n";
    file << "[ATTR]\n\"GEO::Mesh::facet_corners\"\n\"GEO::Mesh::facet_corners::corner_adjacent_facet\"\n\"index_t\"\n" << sizeof(index_t) << "\n1\n";
    MeshConnectivity fec(m);
    for (int c=0; c<m.ncorners(); c++) {
        int opp = fec.opposite(c);
        file << (opp < 0 ? index_t(-1) : fec.c2f[opp]) << "\n";
    }

    for (auto &it : pattr) {
        std::string name = it.first;
        file << "[ATTR]\n\"GEO::Mesh::vertices\"\n\"" << name << "\"\n";

        std::shared_ptr<GenericAttributeContainer> &ptr = it.second;
        if (auto aint = std::dynamic_pointer_cast<AttributeContainer<int> >(ptr); aint.get()!=nullptr) {
            file << "\"int\"\n" << sizeof(int) << " # this is the size of an element (in bytes)\n1 # this is the number of elements per item\n";
            for (int i=0; i<m.nverts(); i++)
                file << aint->data[i] << "\n";
        } else if (auto adouble = std::dynamic_pointer_cast<AttributeContainer<double> >(ptr); adouble.get()!=nullptr) {
            file << "\"double\"\n" << sizeof(double) << " # this is the size of an element (in bytes)\n1 # this is the number of elements per item\n";
            for (int i=0; i<m.nverts(); i++)
                file << adouble->data[i] << "\n";
        } else {
            assert(false);
        }
    }

    for (auto &it : fattr) {
        std::string name = it.first;
        file << "[ATTR]\n\"GEO::Mesh::facets\"\n\"" << name << "\"\n";

        std::shared_ptr<GenericAttributeContainer> &ptr = it.second;
        if (auto aint = std::dynamic_pointer_cast<AttributeContainer<int> >(ptr); aint.get()!=nullptr) {
            file << "\"int\"\n" << sizeof(int) << " # this is the size of an element (in bytes)\n1 # this is the number of elements per item\n";
            for (int i=0; i<m.nfacets(); i++)
                file << aint->data[i] << "\n";
        } else if (auto adouble = std::dynamic_pointer_cast<AttributeContainer<double> >(ptr); adouble.get()!=nullptr) {
            file << "\"double\"\n" << sizeof(double) << " # this is the size of an element (in bytes)\n1 # this is the number of elements per item\n";
            for (int i=0; i<m.nfacets(); i++)
                file << adouble->data[i] << "\n";
        } else {
            assert(false);
        }
    }

    for (auto &it : cattr) {
        std::string name = it.first;
        file << "[ATTR]\n\"GEO::Mesh::facet_corners\"\n\"" << name << "\"\n";

        std::shared_ptr<GenericAttributeContainer> &ptr = it.second;
        if (auto aint = std::dynamic_pointer_cast<AttributeContainer<int> >(ptr); aint.get()!=nullptr) {
            file << "\"int\"\n" << sizeof(int) << " # this is the size of an element (in bytes)\n1 # this is the number of elements per item\n";
            for (int i=0; i<m.ncorners(); i++)
                file << aint->data[i] << "\n";
        } else if (auto adouble = std::dynamic_pointer_cast<AttributeContainer<double> >(ptr); adouble.get()!=nullptr) {
            file << "\"double\"\n" << sizeof(double) << " # this is the size of an element (in bytes)\n1 # this is the number of elements per item\n";
            for (int i=0; i<m.ncorners(); i++)
                file << adouble->data[i] << "\n";
        } else {
            assert(false);
        }
    }

    file.close();
}
