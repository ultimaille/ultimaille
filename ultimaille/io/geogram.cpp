#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "ultimaille/io/geogram.h"
#include <zlib/zlib.h>

namespace UM {
    typedef unsigned int index_t;

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

        template <typename T> void addAttribute(std::string const &wh, std::string const &name, std::string const &type, T* const data, const int nb_items, const int dim) {
            addHeader("ATTR");
            addU64((4+wh.length())+(4+name.length())+(4+type.length())+4+4+sizeof(T)*dim*nb_items);
            addString(wh);
            addString(name);
            addString(type);
            addU32(sizeof(T));
            addU32(dim);
            addData(static_cast<void const *>(data), sizeof(T)*nb_items*dim);
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

    void write_geogram(const std::string filename, const PolyLine &pl, PolyLineAttributes attr) {
        try {
            GeogramGZWriter writer(filename);
            writer.addFileHeader();
            writer.addAttributeSize("GEO::Mesh::vertices", pl.nverts());
            writer.addAttribute("GEO::Mesh::vertices", "point", "double", reinterpret_cast<const double *>(pl.points.data->data()), pl.nverts(), 3);

            writer.addAttributeSize("GEO::Mesh::edges", pl.nsegments());
            {
                std::vector<index_t> segments;
                for (int s : pl.segments) segments.push_back(s);
                writer.addAttribute("GEO::Mesh::edges", "GEO::Mesh::edges::edge_vertex", "index_t", reinterpret_cast<const index_t *>(segments.data()), pl.nsegments(), 2);
            }

            std::vector<NamedContainer> A[2] = {std::get<0>(attr), std::get<1>(attr)};
            for (int z=0; z<2; z++) {
                auto &att = A[z];

                for (int i=0; i<static_cast<int>(att.size()); i++) {
                    std::string name = att[i].first;
                    std::shared_ptr<GenericAttributeContainer> ptr = att[i].second;
                    std::string place = "";

                    if (z==0)
                        place = "GEO::Mesh::vertices";
                    else if (z==1)
                        place = "GEO::Mesh::edges";

                    // TODO externalize that
                    if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<int> >(ptr); cont_ptr.get()!=nullptr) {
                        writer.addAttribute(place, name, "int", cont_ptr->data.data(), cont_ptr->data.size(), 1);
                    } else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<double>>(ptr); cont_ptr.get()!=nullptr) {
                        writer.addAttribute(place, name, "double", cont_ptr->data.data(), cont_ptr->data.size(), 1);
                    } else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<vec2>>(ptr); cont_ptr.get()!=nullptr) {
                        writer.addAttribute(place, name, "double", reinterpret_cast<const double *>(cont_ptr->data.data()), cont_ptr->data.size(), 2);
                    } else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<vec3>>(ptr); cont_ptr.get()!=nullptr) {
                        writer.addAttribute(place, name, "double", reinterpret_cast<const double *>(cont_ptr->data.data()), cont_ptr->data.size(), 3);
                    } else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<bool>>(ptr); cont_ptr.get()!=nullptr) {
                        std::vector<char> tmp(cont_ptr->data.size());
                        for (int i=0; i<(int)cont_ptr->data.size(); i++) tmp[i] = cont_ptr->data[i];
                        writer.addAttribute(place, name, "bool", reinterpret_cast<const char *>(tmp.data()), tmp.size(), 1);
                    } else {
                        assert(false);
                    }
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Ooops: catch error= " << e.what() << " when creating " << filename << "\n";
        }
    }

    void write_geogram(const std::string filename, const Surface &m, SurfaceAttributes attr) {
        try {
            GeogramGZWriter writer(filename);

            if (!m.nverts()) return;

            writer.addFileHeader();
            writer.addAttributeSize("GEO::Mesh::vertices", m.nverts());
            writer.addAttribute("GEO::Mesh::vertices", "point", "double", reinterpret_cast<const double *>(m.points.data->data()), m.nverts(), 3);

            writer.addAttributeSize("GEO::Mesh::facets", m.nfacets());
            std::vector<index_t> facet_ptr;
            for (int f=0; f<m.nfacets(); f++) facet_ptr.push_back(m.corner(f,0));
            writer.addAttribute("GEO::Mesh::facets", "GEO::Mesh::facets::facet_ptr", "index_t", facet_ptr.data(), m.nfacets(), 1);

            writer.addAttributeSize("GEO::Mesh::facet_corners", m.ncorners());
            std::vector<index_t> corner_vertex;
            for (int f=0; f<m.nfacets(); f++)
                for (int v=0; v<m.facet_size(f); v++)
                    corner_vertex.push_back(m.vert(f,v));
            writer.addAttribute("GEO::Mesh::facet_corners", "GEO::Mesh::facet_corners::corner_vertex", "index_t", corner_vertex.data(), m.ncorners(), 1);

            std::vector<index_t> corner_adjacent_facet;
            SurfaceConnectivity fec(m);
            for (int c=0; c<m.ncorners(); c++) {
                int opp = fec.opposite(c);
                corner_adjacent_facet.push_back(opp < 0 ? index_t(-1) : fec.c2f[opp]);
            }
            writer.addAttribute("GEO::Mesh::facet_corners", "GEO::Mesh::facet_corners::corner_adjacent_facet", "index_t", corner_adjacent_facet.data(), m.ncorners(), 1);

            std::vector<NamedContainer> A[3] = {std::get<0>(attr), std::get<1>(attr), std::get<2>(attr)};
            for (int z=0; z<3; z++) {
                auto &att = A[z];

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

                    // TODO externalize that
                    if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<int> >(ptr); cont_ptr.get()!=nullptr) {
                        writer.addAttribute(place, name, "int", cont_ptr->data.data(), cont_ptr->data.size(), 1);
                    } else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<double>>(ptr); cont_ptr.get()!=nullptr) {
                        writer.addAttribute(place, name, "double", cont_ptr->data.data(), cont_ptr->data.size(), 1);
                    } else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<vec2>>(ptr); cont_ptr.get()!=nullptr) {
                        writer.addAttribute(place, name, "double", reinterpret_cast<const double *>(cont_ptr->data.data()), cont_ptr->data.size(), 2);
                    } else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<vec3>>(ptr); cont_ptr.get()!=nullptr) {
                        writer.addAttribute(place, name, "double", reinterpret_cast<const double *>(cont_ptr->data.data()), cont_ptr->data.size(), 3);
                    } else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<bool>>(ptr); cont_ptr.get()!=nullptr) {
                        std::vector<char> tmp(cont_ptr->data.size());
                        for (int i=0; i<(int)cont_ptr->data.size(); i++) tmp[i] = cont_ptr->data[i];
                        writer.addAttribute(place, name, "bool", reinterpret_cast<const char *>(tmp.data()), tmp.size(), 1);
                    } else {
                        assert(false);
                    }
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Ooops: catch error= " << e.what() << " when creating " << filename << "\n";
        }
    }

    const int   geogram_nb_verts_per_cell_type[5] = {4, 8, 6, 5, 4}; // geogram types: tet, hex, prism, pyramid, connector
    const int  geogram_nb_facets_per_cell_type[5] = {4, 6, 5, 5, 3}; //
    const int geogram_nb_padding_per_cell_type[5] = {0, 2, 1, 0, 1}; // geogram cell facet padding. AARGH Bruno!

    template <typename T> void pad_attribute(const Volume &m, std::vector<T> &tmp) {
        int fct_per_cell = geogram_nb_facets_per_cell_type[m.cell_type()];
        int bruno_fct_per_cell = fct_per_cell + geogram_nb_padding_per_cell_type[m.cell_type()];
        int bruno_nfacets = m.ncells() * bruno_fct_per_cell;
//      std::cerr << "PADDING; bruno_nfacets: " << bruno_nfacets << " nfacets: " << m.nfacets() << std::endl;
        tmp.resize(bruno_nfacets);
        for (int c=m.ncells(); c--;)
            for (int f=fct_per_cell; f--;)
                tmp[c*bruno_fct_per_cell+f] = tmp[c*fct_per_cell+f];
    }

    void write_geogram(const std::string filename, const Volume &m, VolumeAttributes attr) {
        try {
            GeogramGZWriter writer(filename);

            if (!m.nverts()) return;

            writer.addFileHeader();
            writer.addAttributeSize("GEO::Mesh::vertices", m.nverts());
            writer.addAttribute("GEO::Mesh::vertices", "point", "double", reinterpret_cast<const double *>(m.points.data->data()), m.nverts(), 3);

            if (!m.ncells()) return;

            writer.addAttributeSize("GEO::Mesh::cells", m.ncells());

            std::vector<char> cell_type(m.ncells());
            for (int c=0; c<m.ncells(); c++) cell_type[c] = m.cell_type();
            writer.addAttribute("GEO::Mesh::cells", "GEO::Mesh::cells::cell_type", "char", cell_type.data(),  m.ncells(), 1);

            std::vector<index_t> cell_ptr(m.ncells());
            for (int c=1; c<m.ncells(); c++) cell_ptr[c] = cell_ptr[c-1]+m.nverts_per_cell();
            writer.addAttribute("GEO::Mesh::cells", "GEO::Mesh::cells::cell_ptr", "index_t", cell_ptr.data(), m.ncells(), 1);

            writer.addAttributeSize("GEO::Mesh::cell_corners", m.ncorners());
            std::vector<index_t> corner_vertex;
            for (int c=0; c<m.ncells(); c++)
                for (int v=0; v<m.nverts_per_cell(); v++)
                    corner_vertex.push_back(m.vert(c,v));
            writer.addAttribute("GEO::Mesh::cell_corners", "GEO::Mesh::cell_corners::corner_vertex", "index_t", corner_vertex.data(), corner_vertex.size(), 1);

            // TODO GEO::Mesh::cell_facets::adjacent_cell

            int bruno_nfacets = m.nfacets()+m.ncells()*geogram_nb_padding_per_cell_type[m.cell_type()]; // AARGH Bruno!
            writer.addAttributeSize("GEO::Mesh::cell_facets", bruno_nfacets);
            std::vector<NamedContainer> A[4] = {std::get<0>(attr), std::get<1>(attr), std::get<2>(attr), std::get<3>(attr)};
            for (int z=0; z<4; z++) {
                auto &att = A[z];

                for (int i=0; i<static_cast<int>(att.size()); i++) {
                    std::string name = att[i].first;
                    std::shared_ptr<GenericAttributeContainer> ptr = att[i].second;
                    std::string place = "";

                    if (z==0)
                        place = "GEO::Mesh::vertices";
                    else if (z==1)
                        place = "GEO::Mesh::cells";
                    else if (z==2)
                        place = "GEO::Mesh::cell_facets";
                    else
                        place = "GEO::Mesh::cell_corners";

                    if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<int> >(ptr); cont_ptr.get()!=nullptr) {
                        std::vector<int> tmp = cont_ptr->data;
                        if (2==z) pad_attribute(m, tmp);
                        writer.addAttribute(place, name, "int", tmp.data(), tmp.size(), 1);
                    } else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<double>>(ptr); cont_ptr.get()!=nullptr) {
                        std::vector<double> tmp = cont_ptr->data;
                        if (2==z) pad_attribute(m, tmp);
                        writer.addAttribute(place, name, "double", tmp.data(), tmp.size(), 1);
                    } else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<vec2>>(ptr); cont_ptr.get()!=nullptr) {
                        std::vector<vec2> tmp = cont_ptr->data;
                        if (2==z) pad_attribute(m, tmp);
                        writer.addAttribute(place, name, "double", reinterpret_cast<const double *>(tmp.data()), tmp.size(), 2);
                    } else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<vec3>>(ptr); cont_ptr.get()!=nullptr) {
                        std::vector<vec3> tmp = cont_ptr->data;
                        if (2==z) pad_attribute(m, tmp);
                        writer.addAttribute(place, name, "double", reinterpret_cast<const double *>(tmp.data()), tmp.size(), 3);
                    } else if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<bool>>(ptr); cont_ptr.get()!=nullptr) {
                        std::vector<char> tmp(cont_ptr->data.size());
                        for (int i=0; i<(int)cont_ptr->data.size(); i++) tmp[i] = cont_ptr->data[i];
                        if (2==z) pad_attribute(m, tmp);
                        writer.addAttribute(place, name, "bool", reinterpret_cast<const char *>(tmp.data()), tmp.size(), 1);
                    } else {
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

    const std::string attrib_set_names[7] = {"GEO::Mesh::vertices", "GEO::Mesh::edges", "GEO::Mesh::facets", "GEO::Mesh::facet_corners", "GEO::Mesh::cells", "GEO::Mesh::cell_facets", "GEO::Mesh::cell_corners"};
    void read_geogram(const std::string filename, std::vector<NamedContainer> attr[7]) {
        int set_size[7] = {-1, -1, -1, -1, -1, -1, -1};
        try {
            GeogramGZReader in(filename);
            std::string chunk_class;
            for (chunk_class=in.next_chunk(); chunk_class!="EOFL"; chunk_class=in.next_chunk()) {
                if (chunk_class == "ATTS") {
                    std::string attribute_set_name = in.read_string();
                    index_t nb_items = in.read_int();
                    in.check_chunk_size();
                    std::cerr << "ATTS " << attribute_set_name << " " << nb_items << std::endl;
                    for (int i=0; i<7; i++)
                        if (attribute_set_name == attrib_set_names[i])
                            set_size[i] = nb_items;
                } else if (chunk_class == "ATTR") {
                    std::string attribute_set_name = in.read_string();
                    int nb_items = -1;
                    for (int i=0; i<7; i++)
                        if (attribute_set_name == attrib_set_names[i])
                            nb_items = set_size[i];

                    assert(nb_items>0);

                    std::string attribute_name = in.read_string();
                    std::string element_type   = in.read_string();
                    index_t element_size = in.read_int();
                    index_t dimension    = in.read_int();
                    size_t size = size_t(element_size) * size_t(dimension) * size_t(nb_items);

                    if (attribute_name=="GEO::Mesh::cell_facets::adjacent_cell" || attribute_name=="GEO::Mesh::facet_corners::corner_adjacent_facet") continue;

                    std::cerr << "ATTR " << attribute_set_name << " " << attribute_name << " " << element_type << " " << element_size << " " << dimension << "\n";

                    std::shared_ptr<GenericAttributeContainer> P;
                    if (element_type=="char") {
                        assert(dimension == 1);
                        std::vector<char> tmp(nb_items);
                        in.read_attribute(tmp.data(), size);
                        GenericAttribute<int> A(nb_items);
                        for (int i=0; i<nb_items; i++) {
                            A[i] = tmp[i];
                        }
                        P = A.ptr;
                    } else if (element_type=="int" || element_type=="index_t" || element_type=="signed_index_t") {
                        GenericAttribute<int> A(nb_items);
                        if (attribute_name=="GEO::Mesh::edges::edge_vertex")
                            A.ptr->data.resize(nb_items*2); // TODO AARGH Bruno!
                        assert(dimension == 1 || (attribute_name=="GEO::Mesh::edges::edge_vertex" && dimension == 2));
                        void *ptr = std::dynamic_pointer_cast<AttributeContainer<int> >(A.ptr)->data.data();
                        in.read_attribute(ptr, size);
                        P = A.ptr;
                    } else if (element_type=="double" && 1==dimension) {
                        GenericAttribute<double> A(nb_items);
                        in.read_attribute(std::dynamic_pointer_cast<AttributeContainer<double> >(A.ptr)->data.data(), size);
                        P = A.ptr;
                    } else if ((element_type=="vec2" && 1==dimension) || (element_type=="double" && 2==dimension)) {
                        GenericAttribute<vec2> A(nb_items);
                        in.read_attribute(std::dynamic_pointer_cast<AttributeContainer<vec2> >(A.ptr)->data.data(), size);
                        P = A.ptr;
                    } else if ((element_type=="vec3" && 1==dimension) || (element_type=="double" && 3==dimension)) {
                        GenericAttribute<vec3> A(nb_items);
                        in.read_attribute(std::dynamic_pointer_cast<AttributeContainer<vec3> >(A.ptr)->data.data(), size);
                        P = A.ptr;
                    } else if (element_type=="bool" && 1==dimension) {
                        std::vector<char> tmp(nb_items, 0);
                        in.read_attribute(tmp.data(), size);
                        GenericAttribute<bool> A(nb_items);
                        for (int i=0; i<nb_items; i++) A[i] = tmp[i];
                        P = A.ptr;
                    } else {
                        continue;
                    }

                    for (int i=0; i<7; i++)
                        if (attribute_set_name == attrib_set_names[i])
                            attr[i].emplace_back(attribute_name, P);
                } // chunk_class = ATTR
            } // chunks
        } catch (const std::exception& e) {
            std::cerr << "Ooops: catch error= " << e.what() << " when reading " << filename << "\n";
        }
    }

    void parse_pointset_attributes(PointSet &pts, std::vector<NamedContainer> &attr) {
        for (int i=0; i<(int)attr.size(); i++) {
            if (attr[i].first != "point") continue;
            std::shared_ptr<AttributeContainer<vec3> > ptr = std::dynamic_pointer_cast<AttributeContainer<vec3> >(attr[i].second);
            pts.resize(ptr->data.size());
            for (int v=0; v<pts.size(); v++)
                pts[v] = ptr->data[v];
            attr.erase(attr.begin()+i);
            i--;
        }
    }

    void parse_int_array(const std::string &name, std::vector<int> &array, std::vector<NamedContainer> &attr) {
        for (int i=0; i<(int)attr.size(); i++) {
            if (attr[i].first != name) continue;
            std::shared_ptr<AttributeContainer<int> > ptr = std::dynamic_pointer_cast<AttributeContainer<int> >(attr[i].second);
            array = ptr->data;
            attr.erase(attr.begin()+i);
            i--;
        }
    }

    void parse_volume_data(const std::string filename, PointSet &pts, VolumeAttributes &va, std::vector<int> &corner_vertex, int type2keep) {
        std::vector<NamedContainer> attrib[7];
        read_geogram(filename, attrib);
        parse_pointset_attributes(pts, attrib[0]);

        std::vector<int> old_corner_vertex;
        parse_int_array("GEO::Mesh::cell_corners::corner_vertex", old_corner_vertex, attrib[6]);
        std::vector<int> cell_type(old_corner_vertex.size()/4, 0); // if cell_type and cell_ptr are not present in .geogram, create it (tet mesh)
        std::vector<int> cell_ptr(old_corner_vertex.size()/4, 0);
        for (int i=0; i<(int)cell_type.size(); i++) cell_ptr[i] = i*4;

        parse_int_array("GEO::Mesh::cells::cell_type", cell_type, attrib[4]);
        parse_int_array("GEO::Mesh::cells::cell_ptr",  cell_ptr,  attrib[4]);

        int ncells = cell_type.size();
        assert(cell_ptr.size()==cell_type.size());
        cell_ptr.push_back(old_corner_vertex.size());

        corner_vertex = std::vector<int>();
        corner_vertex.reserve(old_corner_vertex.size());

        std::vector<int> cells_old2new(ncells,  -1);
        int nfacets = 0;
        int ncorners = 0;
        for (int c=0; c<ncells; c++) {
            assert(cell_type[c]>=0 && cell_type[c]<=6);
            nfacets  += geogram_nb_facets_per_cell_type[cell_type[c]] + geogram_nb_padding_per_cell_type[cell_type[c]];
            ncorners +=  geogram_nb_verts_per_cell_type[cell_type[c]];
        }
        assert(ncorners == (int)old_corner_vertex.size());
        std::vector<int> corners_old2new(ncorners, -1);
        std::vector<int> facets_old2new(nfacets, -1);

        int new_ncells = 0;
        int new_ncorners = 0;
        int cur_corner = 0;
        int cur_facet = 0;
        int new_nfacets = 0;

        for (int c=0; c<ncells; c++) {
            if (type2keep!=cell_type[c]) continue;

            for (int v=cell_ptr[c]; v<cell_ptr[c+1]; v++)
                corner_vertex.push_back(old_corner_vertex[v]);

            cells_old2new[c] = new_ncells++;

            for (int i=0; i<geogram_nb_facets_per_cell_type[cell_type[c]]; i++)
                facets_old2new[cur_facet+i] = new_nfacets++;
            cur_facet += geogram_nb_facets_per_cell_type[cell_type[c]] + geogram_nb_padding_per_cell_type[cell_type[c]];

            for (int i=0; i<geogram_nb_verts_per_cell_type[cell_type[c]]; i++)
                corners_old2new[cur_corner+i] = new_ncorners++;
            cur_corner += geogram_nb_verts_per_cell_type[cell_type[c]];
        }

        for (auto &nc : attrib[4])
            (*nc.second).compress(cells_old2new);
        for (auto &nc : attrib[5])
            (*nc.second).compress(facets_old2new);
        for (auto &nc : attrib[6])
            (*nc.second).compress(corners_old2new);

        std::get<0>(va) = attrib[0];
        std::get<1>(va) = attrib[4];
        std::get<2>(va) = attrib[5];
        std::get<3>(va) = attrib[6];
    }

    VolumeAttributes read_geogram(const std::string filename, Tetrahedra &m) {
        m = Tetrahedra();
        VolumeAttributes va;
        std::vector<int> corner_vertex;
        parse_volume_data(filename, m.points, va, corner_vertex, 0);
        assert(corner_vertex.size()%4==0);

        int ntetra = corner_vertex.size()*4;
        m.create_cells(ntetra);
        m.cells = corner_vertex;

        for (auto &a : std::get<0>(va)) m.points.attr.emplace_back(a.second);
        for (auto &a : std::get<1>(va)) m.attr_cells.emplace_back(a.second);
        for (auto &a : std::get<2>(va)) m.attr_facets.emplace_back(a.second);
        for (auto &a : std::get<3>(va)) m.attr_corners.emplace_back(a.second);

        return va;
    }

    VolumeAttributes read_geogram(const std::string filename, Hexahedra &m) {
        m = Hexahedra();
        VolumeAttributes va;
        std::vector<int> corner_vertex;
        parse_volume_data(filename, m.points, va, corner_vertex, 1);
        assert(corner_vertex.size()%8==0);

        int nhexa = corner_vertex.size()*8;
        m.create_cells(nhexa);
        m.cells = corner_vertex;

        for (auto &a : std::get<0>(va)) m.points.attr.emplace_back(a.second);
        for (auto &a : std::get<1>(va)) m.attr_cells.emplace_back(a.second);
        for (auto &a : std::get<2>(va)) m.attr_facets.emplace_back(a.second);
        for (auto &a : std::get<3>(va)) m.attr_corners.emplace_back(a.second);

        return va;
    }

    SurfaceAttributes read_geogram(const std::string filename, Polygons &polygons) {
        polygons = Polygons();

        std::vector<NamedContainer> attrib[7];
        read_geogram(filename, attrib);
        parse_pointset_attributes(polygons.points, attrib[0]);

        parse_int_array("GEO::Mesh::facet_corners::corner_vertex", polygons.facets, attrib[3]);

        polygons.offset.resize(polygons.facets.size()/3);// if facet_ptr is not present in .geogram, create it (tri mesh)
        for (int i=0; i<(int)polygons.offset.size(); i++) polygons.offset[i] = i*3;
        parse_int_array("GEO::Mesh::facets::facet_ptr", polygons.offset, attrib[2]);
        polygons.offset.push_back(polygons.facets.size());

        for (auto &a : attrib[0]) polygons.points.attr.emplace_back(a.second);
        for (auto &a : attrib[2]) polygons.attr_facets.emplace_back(a.second);
        for (auto &a : attrib[3]) polygons.attr_corners.emplace_back(a.second);

        return make_tuple(attrib[0], attrib[2], attrib[3]);
    }

    SurfaceAttributes read_geogram(const std::string filename, Triangles &m) {
        Polygons mpoly;
        SurfaceAttributes polyattr = read_geogram(filename, mpoly);
        std::vector<bool> to_kill(mpoly.nfacets(), false);
        for (int f=0; f<mpoly.nfacets(); f++)
            to_kill[f] = (3!=mpoly.facet_size(f));
        mpoly.delete_facets(to_kill);
        m.points = mpoly.points;
        m.facets = mpoly.facets;
        m.attr_facets = mpoly.attr_facets;
        m.attr_corners = mpoly.attr_corners;

        return polyattr;
    }

    SurfaceAttributes read_geogram(const std::string filename, Quads &m) {
        Polygons mpoly;
        SurfaceAttributes polyattr = read_geogram(filename, mpoly);
        std::vector<bool> to_kill(mpoly.nfacets(), false);
        for (int f=0; f<mpoly.nfacets(); f++)
            to_kill[f] = (4!=mpoly.facet_size(f));
        mpoly.delete_facets(to_kill);
        m.points = mpoly.points;
        m.facets = mpoly.facets;
        m.attr_facets = mpoly.attr_facets;
        m.attr_corners = mpoly.attr_corners;

        return polyattr;
    }

    PolyLineAttributes read_geogram(const std::string filename, PolyLine &pl) {
        pl = PolyLine();

        std::vector<NamedContainer> attrib[7];
        read_geogram(filename, attrib);
        parse_pointset_attributes(pl.points, attrib[0]);

        parse_int_array("GEO::Mesh::edges::edge_vertex", pl.segments, attrib[1]);

        for (auto &a : attrib[0]) pl.points.attr.emplace_back(a.second);
        for (auto &a : attrib[1]) pl.attr.emplace_back(a.second);

        return make_tuple(attrib[0], attrib[1]);
    }
}

