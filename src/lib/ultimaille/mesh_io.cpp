#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "surface.h"
#include "attributes.h"
#include "mesh_io.h"

#include <zlib/zlib.h>


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
    file << "[ATTR]\n\"GEO::Mesh::facets\"\n\"GEO::Mesh::facets::facet_ptr\"\n\"index_t\"\n4\n1\n";
    for (int f=0; f<m.nfacets(); f++)
        file << m.facet_corner(f,0) << "\n";

    file << "[ATTS]\n\"GEO::Mesh::facet_corners\"\n";
    file << m.ncorners() << "\n";
    file << "[ATTR]\n\"GEO::Mesh::facet_corners\"\n\"GEO::Mesh::facet_corners::corner_vertex\"\n\"index_t\"\n4\n1\n";
    for (int f=0; f<m.nfacets(); f++)
        for (int v=0; v<m.facet_size(f); v++)
           file << m.vert(f,v) << "\n";
    file << "[ATTR]\n\"GEO::Mesh::facet_corners\"\n\"GEO::Mesh::facet_corners::corner_adjacent_facet\"\n\"index_t\"\n4\n1\n";
    MeshConnectivity fec(m);
    for (int c=0; c<m.ncorners(); c++) {
        int opp = fec.opposite(c);
        file << (opp < 0 ? index_t(-1) : fec.c2f[opp]) << "\n";
    }
    auto attr = pattr;
    attr.insert(attr.end(), fattr.begin(), fattr.end());
    attr.insert(attr.end(), cattr.begin(), cattr.end());

    int npattr = pattr.size();
    int nfattr = fattr.size();

    for (int i=0; i<static_cast<int>(attr.size()); i++) {
        std::string name = attr[i].first;
        std::shared_ptr<GenericAttributeContainer> ptr = attr[i].second;

        if (i<npattr)
            file << "[ATTR]\n\"GEO::Mesh::vertices\"\n\"" << name << "\"\n";
        else if (i<npattr+nfattr)
            file << "[ATTR]\n\"GEO::Mesh::facets\"\n\"" << name << "\"\n";
        else
            file << "[ATTR]\n\"GEO::Mesh::facet_corners\"\n\"" << name << "\"\n";

        if (auto aint = std::dynamic_pointer_cast<AttributeContainer<int> >(ptr); aint.get()!=nullptr) {
            file << "\"int\"\n4\n1\n";
            for (const auto &v : aint->data)
                file << v << "\n";
        } else if (auto adouble = std::dynamic_pointer_cast<AttributeContainer<double> >(ptr); adouble.get()!=nullptr) {
            file << "\"double\"\n8\n1\n";
            for (const auto &v : adouble->data)
                file << v << "\n";
        } else {
            assert(false);
        }
    }

    file.close();
}


namespace {
    struct GeogramGZWriter {
        GeogramGZWriter(std::string const &fName) : m_gzFile() {
            m_gzFile = gzopen(fName.c_str(), "wb");
            if (!m_gzFile)
                throw "Can not open file";
        }

        ~GeogramGZWriter() {
            if (m_gzFile)
                gzclose(m_gzFile);
        }

        void addFileHeader() {
            addHeader("HEAD");
            addU64(4+7+4+3);
            addString("GEOGRAM");
            addString("1.0");
        }

        void addAttributeSize(std::string const &wh, size_t n) {
            addHeader("ATTS");
            addU64(4+wh.length()+4);
            addString(wh);
            addU32(uint32_t(n));
        }

        template <typename T> void addAttribute(std::string const &wh, std::string const &name, std::string const &type, std::vector<T> const &data) {
            addHeader("ATTR");
            addU64((4+wh.length())+(4+name.length())+(4+type.length())+4+4+sizeof(T)*data.size());
            addString(wh);
            addString(name);
            addString(type);
            addU32(sizeof(T));
            addU32(1);
            addData(static_cast<void const *>(data.data()), sizeof(T)*data.size());
        }

        void addAttribute(std::string const &wh, std::string const &name, std::vector<vec3> const &data) {
            addHeader("ATTR");
            addU64((4+wh.length())+(4+name.length())+(4+6)+4+4+sizeof(double)*3*data.size());
            addString(wh);
            addString(name);
            addString("double");
            addU32(8);
            addU32(3);
            std::vector<double> values;
            for (const auto &p : data) {
                values.push_back(p.x);
                values.push_back(p.y);
                values.push_back(p.z);
            }
            addData(static_cast<void const *>(values.data()), sizeof(double)*values.size());
        }

        protected:
        void addU64(uint64_t size) {
            addData(&size, 8);
        }

        void addU32(uint32_t index) {
            addData(&index, 4);
        }

        void addHeader(std::string const &header) {
            if (header.length()!=4)
                throw "GeogramGZWriter::bad header";
            addData(static_cast<void const *>(header.c_str()), 4);
        }

        void addString(std::string const &str) {
            size_t len = str.length();
            addU32(uint32_t(len));
            addData(str.c_str(),len);
        }

        void addData(void const *data, size_t len) {
            int check = gzwrite(m_gzFile, data, (unsigned int)(len));
            if (size_t(check) != len)
                throw "Could not write attribute data";
        }

        gzFile m_gzFile;
    };
}

void write_geogram(const std::string filename, const Surface &m,
        std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > pattr,
        std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > fattr,
        std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > cattr) {
    try {
        GeogramGZWriter writer(filename);
        writer.addFileHeader();
        writer.addAttributeSize("GEO::Mesh::vertices", m.nverts());
        writer.addAttribute("GEO::Mesh::vertices", "point", *m.points.data);

        writer.addAttributeSize("GEO::Mesh::facets", m.nfacets());
        std::vector<index_t> facet_ptr;
        for (int f=0; f<m.nfacets(); f++) facet_ptr.push_back(m.facet_corner(f,0));
        writer.addAttribute("GEO::Mesh::facets", "GEO::Mesh::facets::facet_ptr", "index_t", facet_ptr);

        writer.addAttributeSize("GEO::Mesh::facet_corners", m.ncorners());
        std::vector<index_t> corner_vertex;
        for (int f=0; f<m.nfacets(); f++)
            for (int v=0; v<m.facet_size(f); v++)
                corner_vertex.push_back(m.vert(f,v));
        writer.addAttribute("GEO::Mesh::facet_corners", "GEO::Mesh::facet_corners::corner_vertex", "index_t", corner_vertex);

        std::vector<index_t> corner_adjacent_facet;
        MeshConnectivity fec(m);
        for (int c=0; c<m.ncorners(); c++) {
            int opp = fec.opposite(c);
            corner_adjacent_facet.push_back(opp < 0 ? index_t(-1) : fec.c2f[opp]);
        }
        writer.addAttribute("GEO::Mesh::facet_corners", "GEO::Mesh::facet_corners::corner_adjacent_facet", "index_t", corner_adjacent_facet);

        auto attr = pattr;
        attr.insert(attr.end(), fattr.begin(), fattr.end());
        attr.insert(attr.end(), cattr.begin(), cattr.end());

        int npattr = pattr.size();
        int nfattr = fattr.size();

        for (int i=0; i<static_cast<int>(attr.size()); i++) {
            std::string name = attr[i].first;
            std::shared_ptr<GenericAttributeContainer> ptr = attr[i].second;
            std::string place = "";

            if (i<npattr)
                place = "GEO::Mesh::vertices";
            else if (i<npattr+nfattr)
                place = "GEO::Mesh::facets";
            else
                place = "GEO::Mesh::facet_corners";

            if (auto aint = std::dynamic_pointer_cast<AttributeContainer<int> >(ptr); aint.get()!=nullptr) {
                writer.addAttribute(place, name, "int", aint->data);
            } else if (auto adouble = std::dynamic_pointer_cast<AttributeContainer<double> >(ptr); adouble.get()!=nullptr) {
                writer.addAttribute(place, name, "double", adouble->data);
            } else {
                assert(false);
            }
        }
    } catch (char const *error) {
        std::cerr << "Ooops: catch error=" << error << " when creating " << filename << "\n";
    }
}


