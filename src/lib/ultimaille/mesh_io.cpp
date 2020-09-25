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

// Attention: only int and double attributes
typedef unsigned int index_t;
void write_geogram_ascii(const std::string filename, const Surface &m, SurfaceAttributes attr) {
    auto [pattr, fattr, cattr] = attr;
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
    SurfaceConnectivity fec(m);
    for (int c=0; c<m.ncorners(); c++) {
        int opp = fec.opposite(c);
        file << (opp < 0 ? index_t(-1) : fec.c2f[opp]) << "\n";
    }

    // TODO ugly, repair it
    std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > >  A[3] = {std::get<0>(attr), std::get<1>(attr), std::get<2>(attr)};
    for (int i=0; i<3; i++) {
        auto att = A[i];

        for (int i=0; i<static_cast<int>(att.size()); i++) {
            std::string name = att[i].first;
            std::shared_ptr<GenericAttributeContainer> ptr = att[i].second;

            if (i==0)
                file << "[ATTR]\n\"GEO::Mesh::vertices\"\n\"" << name << "\"\n";
            else if (i==1)
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
    }

    file.close();
}


namespace {
    struct GeogramGZWriter {
        GeogramGZWriter(std::string const &fName) : file_() {
            file_ = gzopen(fName.c_str(), "wb");
            if (!file_)
                throw std::runtime_error("Can not open file");
        }

        ~GeogramGZWriter() {
            if (file_)
                gzclose(file_);
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

        void addAttribute(std::string const &wh, std::string const &name, std::vector<vec2> const &data) {
            addHeader("ATTR");
            addU64((4+wh.length())+(4+name.length())+(4+6)+4+4+sizeof(double)*2*data.size());
            addString(wh);
            addString(name);
            addString("double");
            addU32(8);
            addU32(2);
            std::vector<double> values;
            for (const auto &p : data) {
                values.push_back(p.x);
                values.push_back(p.y);
            }
            addData(static_cast<void const *>(values.data()), sizeof(double)*values.size());
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
                throw std::runtime_error("GeogramGZWriter::bad header");
            addData(static_cast<void const *>(header.c_str()), 4);
        }

        void addString(std::string const &str) {
            size_t len = str.length();
            addU32(uint32_t(len));
            addData(str.c_str(),len);
        }

        void addData(void const *data, size_t len) {
            int check = gzwrite(file_, data, (unsigned int)(len));
            if (size_t(check) != len)
                throw std::runtime_error("Could not write attribute data");
        }

        gzFile file_;
    };
}

void write_geogram(const std::string filename, const Surface &m, SurfaceAttributes attr) {
    auto [pattr, fattr, cattr] = attr;
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
        SurfaceConnectivity fec(m);
        for (int c=0; c<m.ncorners(); c++) {
            int opp = fec.opposite(c);
            corner_adjacent_facet.push_back(opp < 0 ? index_t(-1) : fec.c2f[opp]);
        }
        writer.addAttribute("GEO::Mesh::facet_corners", "GEO::Mesh::facet_corners::corner_adjacent_facet", "index_t", corner_adjacent_facet);



        // TODO ugly, repair it
        std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > >  A[3] = {std::get<0>(attr), std::get<1>(attr), std::get<2>(attr)};
        for (int z=0; z<3; z++) {
            auto att = A[z];

            for (int i=0; i<static_cast<int>(att.size()); i++) {
                std::string name = att[i].first;
                std::shared_ptr<GenericAttributeContainer> ptr = att[i].second;
                std::string place = "";

                if (z==0)
                    place = "GEO::Mesh::vertices";
                else if (z==1)
                    place = "GEO::Mesh::facets";
                else
                    place = "GEO::Mesh::facet_corners";

                if (auto aint = std::dynamic_pointer_cast<AttributeContainer<int> >(ptr); aint.get()!=nullptr) {
                    writer.addAttribute(place, name, "int", aint->data);
                } else if (auto adouble = std::dynamic_pointer_cast<AttributeContainer<double> >(ptr); adouble.get()!=nullptr) {
                    writer.addAttribute(place, name, "double", adouble->data);
                } else if (auto avec2 = std::dynamic_pointer_cast<AttributeContainer<vec2> >(ptr); avec2.get()!=nullptr) {
                    writer.addAttribute(place, name, avec2->data);
                } else if (auto avec3 = std::dynamic_pointer_cast<AttributeContainer<vec3> >(ptr); avec3.get()!=nullptr) {
                    writer.addAttribute(place, name, avec3->data);
                } else {
                    //      std::cerr << place << std::endl;
                    //        TODO : ignore adjacency
                    assert(false);
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Ooops: catch error= " << e.what() << " when creating " << filename << "\n";
    }
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct GeogramGZReader {
    GeogramGZReader(std::string const &filename) :
        file_(nullptr),
        current_chunk_class_("0000"),
        current_chunk_size_(0),
        current_chunk_file_pos_(0) {

//        std::cerr << "GeogramGZReader()\n";
        file_ = gzopen(filename.c_str(), "rb");
        if (!file_)
            throw std::runtime_error("Can not open file");

//        std::cerr << "read_chunk_header()\n";
        read_chunk_header();
        if (current_chunk_class_ != "HEAD") {
            throw std::runtime_error(filename + " Does not start with HEAD chunk");
        }

//        std::cerr << "read_string()\n";
        std::string magic = read_string();
        if(magic != "GEOGRAM") {
            throw std::runtime_error(filename + " is not a GEOGRAM file");
        }
        std::string version = read_string();
        //          Logger::out("I/O") << "GeoFile version: " << version << std::endl;
        check_chunk_size();
    }

    ~GeogramGZReader() {
        if (file_)
            gzclose(file_);
    }

    index_t read_int() {
        std::uint32_t result = 0;
        int check = gzread(file_, &result, sizeof(std::uint32_t));
        if (!check && gzeof(file_))
            result = std::uint32_t(-1);
        else if(size_t(check) != sizeof(std::uint32_t))
            throw std::runtime_error("Could not read integer from file");
        return result;
    }

    std::string read_string() {
        std::string result;
        index_t len = read_int();
        result.resize(len);
        if (len) {
            int check = gzread(file_, &result[0], (unsigned int)(len));
            if (index_t(check) != len)
                throw std::runtime_error("Could not read string data from file");
        }
        return result;
    }

    void check_chunk_size() {
        long chunk_size = gztell(file_) - current_chunk_file_pos_;
        if (current_chunk_size_ != chunk_size)
            throw std::runtime_error(std::string("Chunk size mismatch: ") + " expected " + std::to_string(current_chunk_size_) + "/ got " + std::to_string(chunk_size));
    }

    std::string read_chunk_class() {
        std::string result;
        result.resize(4,'\0');
        int check = gzread(file_, &result[0], 4);
        if (check == 0 && gzeof(file_))
            result = "EOFL";
        else if (check != 4)
            throw std::runtime_error("Could not read chunk class from file");
        return result;
    }

    size_t read_size() {
        std::uint64_t result = 0;
        int check = gzread(file_, &result, sizeof(std::uint64_t));
        if (check == 0 && gzeof(file_))
            result = size_t(-1);
        else
            if (size_t(check) != sizeof(std::uint64_t))
                throw std::runtime_error("Could not read size from file");
        return size_t(result);
    }

    void read_chunk_header() {
        current_chunk_class_ = read_chunk_class();
        if (gzeof(file_)) {
            gzclose(file_);
            file_ = nullptr;
            current_chunk_size_ = 0;
            current_chunk_class_ = "EOFL";
            return;
        }
        current_chunk_size_ = long(read_size());
        current_chunk_file_pos_ = gztell(file_);
    }

    void skip_chunk() {
        gzseek(file_, current_chunk_size_ + current_chunk_file_pos_, SEEK_SET);
    }

    const std::string& next_chunk() {
        // If the file pointer did not advance as expected
        // between two consecutive calls of next_chunk, it
        // means that the client code does not want to
        // read the current chunk, then it needs to be
        // skipped.
        if (gztell(file_) != current_chunk_file_pos_ + current_chunk_size_)
            skip_chunk();

        read_chunk_header();
        return current_chunk_class_;
    }

    void read_attribute(void* addr, size_t size) {
        assert(current_chunk_class_ == "ATTR");
        int check = gzread(file_, addr, (unsigned int)(size));
        if (size_t(check) != size)
            throw std::runtime_error("Could not read attribute  (" + std::to_string(check) + "/" + std::to_string(size) + " bytes read)");
        check_chunk_size();
    }

    gzFile file_;
    std::string current_chunk_class_;
    long current_chunk_size_;
    long current_chunk_file_pos_;
};

std::tuple<std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > >,
           std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > >,
           std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > > read_geogram(const std::string filename, Polygons &m) {
    std::vector<std::pair<std::string, std::shared_ptr<GenericAttributeContainer> > > pattr, fattr, cattr;
    m = Polygons();
    try {
        GeogramGZReader in(filename);
        std::string chunk_class;
        for (chunk_class=in.next_chunk(); chunk_class!="EOFL"; chunk_class=in.next_chunk()) {
            if (chunk_class == "ATTS") {
                std::string attribute_set_name = in.read_string();
                index_t nb_items = in.read_int();
                in.check_chunk_size();
                std::cerr << "ATTS " << attribute_set_name << " " << nb_items << std::endl;

                if (attribute_set_name == "GEO::Mesh::vertices") {
                    m.points.resize(nb_items);
                } else if (attribute_set_name == "GEO::Mesh::facets") {
                    m.offset.resize(nb_items+1);
                    for (index_t o=0; o<nb_items+1; o++) m.offset[o] = 3*o;
                } else if (attribute_set_name == "GEO::Mesh::facet_corners") {
                    m.facets.resize(nb_items);
                }
            } else if (chunk_class == "ATTR") {
                std::string attribute_set_name = in.read_string();
                int nb_items = 0;
                if (attribute_set_name == "GEO::Mesh::vertices") {
                    nb_items = m.points.size();
                } else if (attribute_set_name == "GEO::Mesh::facets") {
                    nb_items = m.offset.size()-1;
                } else if (attribute_set_name == "GEO::Mesh::facet_corners") {
                    nb_items = m.facets.size();
                } else {
                    continue;
                }
                assert(nb_items>0);

                std::string attribute_name = in.read_string();
                std::string element_type = in.read_string();
                index_t element_size = in.read_int();
                index_t dimension = in.read_int();
                size_t size = size_t(element_size) * size_t(dimension) * size_t(nb_items);

                std::cerr << "ATTR " << attribute_set_name << " " << attribute_name << " " << element_type << " " << element_size << " " << dimension << "\n";

                if (attribute_set_name == "GEO::Mesh::vertices" && attribute_name == "point") {
                    assert(element_size==8);
                    std::vector<double> raw(nb_items*dimension);
                    in.read_attribute((void *)raw.data(), size);
                    for (int v=0; v<nb_items; v++)
                        m.points[v] = {raw[v*3+0], raw[v*3+1], raw[v*3+2]};
                } else if (attribute_set_name == "GEO::Mesh::facets" && attribute_name == "GEO::Mesh::facets::facet_ptr") {
                    in.read_attribute((void *)m.offset.data(), size);
                } else if (attribute_set_name == "GEO::Mesh::facet_corners" && attribute_name == "GEO::Mesh::facet_corners::corner_vertex") {
                    in.read_attribute((void *)m.facets.data(), size);
                    std::cerr << "size " << size <<std::endl; 
                } else if (attribute_name!="GEO::Mesh::facet_corners::corner_adjacent_facet") {
                    std::shared_ptr<GenericAttributeContainer> P;
                    if (element_type=="int" || element_type=="signed_index_t") {
                        GenericAttribute<int> A(nb_items);
                        void *ptr = std::dynamic_pointer_cast<AttributeContainer<int> >(A.ptr)->data.data();
                        in.read_attribute(ptr, size);
                        P = A.ptr;
                    } else if (element_type=="double") {
                        assert(element_size==8);
                        std::vector<double> raw(nb_items*dimension);
                        in.read_attribute((void *)raw.data(), size);

                        if (1==dimension) {
                            GenericAttribute<double> A(nb_items);
                            for (int i=0; i<nb_items; i++)
                                A[i] = raw[i];
                            P = A.ptr;
                        } if (2==dimension) {
                            GenericAttribute<vec2> A(nb_items);
                            for (int i=0; i<nb_items; i++)
                                    A[i] = {raw[i*2+0], raw[i*2+1]};
                             P = A.ptr;
                       } if (3==dimension) {
                            GenericAttribute<vec3> A(nb_items);
                            for (int i=0; i<nb_items; i++)
                                    A[i] = {raw[i*3+0], raw[i*3+1], raw[i*3+2]};
                             P = A.ptr;
                       }

//                      GenericAttribute<double> A(nb_items);
//                      void *ptr = std::dynamic_pointer_cast<AttributeContainer<double> >(A.ptr)->data.data();
//                      in.read_attribute(ptr, size);
//                      P = A.ptr;
                    }

                    if (attribute_set_name == "GEO::Mesh::vertices") {
                        pattr.emplace_back(attribute_name, P);
                        m.points.attr.emplace_back(P);
                    } else if (attribute_set_name == "GEO::Mesh::facets") {
                        fattr.emplace_back(attribute_name, P);
                        m.attr_facets.emplace_back(P);
                    } else if (attribute_set_name == "GEO::Mesh::facet_corners") {
                        cattr.emplace_back(attribute_name, P);
                        m.attr_corners.emplace_back(P);
                    }
                }
            } // chunk_class = ATTR
        } // chunks
        m.offset.back() = m.facets.size();
    } catch (const std::exception& e) {
        std::cerr << "Ooops: catch error= " << e.what() << " when reading " << filename << "\n";
    }
    return std::make_tuple(pattr, fattr, cattr);
}

