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
//        std::cerr << "read_chunk_class()\n";
        std::string result;
        result.resize(4,'\0');
        int check = gzread(file_, &result[0], 4);
//        std::cerr << "read_chunk_class()\n";
        if (check == 0 && gzeof(file_))
            result = "EOFL";
        else if (check != 4)
            throw std::runtime_error("Could not read chunk class from file");
        return result;
    }

    size_t read_size() {
//        std::cerr << "read_size()\n";
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
//        std::cerr << "read_chunk_header()\n";
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

    gzFile file_;
    std::string current_chunk_class_;
    long current_chunk_size_;
    long current_chunk_file_pos_;
};

void read_geogram(const std::string filename, PolyMesh &m) {
    m = PolyMesh();
    try {
        GeogramGZReader in(filename);
        std::string chunk_class;
        for (chunk_class=in.next_chunk(); chunk_class!="EOFL"; chunk_class=in.next_chunk()) {
            if (chunk_class == "ATTS") {
                std::string attribute_set_name = in.read_string();
                index_t nb_items = in.read_int();
                in.check_chunk_size();
                std::cerr << "ATTS " << attribute_set_name << " " << nb_items << std::endl;
            } else if (chunk_class == "ATTR") {
                std::string attribute_set_name = in.read_string();
                std::string attribute_name = in.read_string();
                std::string element_type = in.read_string();
                index_t element_size = in.read_int();
                index_t dimension = in.read_int();

                std::cerr << "ATTR " << attribute_set_name << " " << attribute_name << " " << element_type << " " << element_size << " " << dimension << "\n";
            }
            /*
               if (chunk_class == "ATTS") {
               read_attribute_set(in, M, ioflags);
               } else if (chunk_class == "ATTR") {
               if (String::string_starts_with(in.current_attribute().name, "GEO::Mesh::")) {
               read_internal_attribute(in, M, ioflags);
               } else {
               read_user_attribute(in, M, ioflags);
               }
               }
             */
        }
    } catch (const std::exception& e) {
        std::cerr << "Ooops: catch error= " << e.what() << " when reading " << filename << "\n";
    }
}

#if 0
    class GEOGRAM_API GeogramIOHandler : public MeshIOHandler {
    public:


        /**
         * \brief Loads a mesh from a GeoFile ('.geogram' file format).
         * \details
         * Loads the contents of the InputGeoFile \p geofile and stores the
         * resulting mesh to \p M. This function can be used to load several
         * meshes that are stored in the same GeoFile.
         * \param[in] in a reference to the InputGeoFile
         * \param[out] M the loaded mesh
         * \param[in] ioflags specifies which attributes and 
         *  elements should be loaded
         * \return true on success, false otherwise.
         */
        virtual bool load(InputGeoFile& in, Mesh& M, const MeshIOFlags& ioflags = MeshIOFlags()) {
            M.clear();
            try {

                std::string chunk_class;
                for(
                        chunk_class = in.next_chunk();
                        chunk_class != "EOFL" && chunk_class != "SPTR";
                        chunk_class = in.next_chunk()
                   ) {

                    if(chunk_class == "ATTS") {
                        read_attribute_set(in, M, ioflags);
                    } else if(chunk_class == "ATTR") {
                        if(
                                String::string_starts_with(
                                    in.current_attribute().name, "GEO::Mesh::"
                                    )
                          ) {
                            read_internal_attribute(in, M, ioflags);
                        } else {
                            read_user_attribute(in, M, ioflags);
                        }
                    } 
                }

                // Create facet "sentry"
                if(!M.facets.are_simplices()) {
                    M.facets.facet_ptr_[M.facets.nb()] = M.facet_corners.nb();
                }

                // Create cell "sentry"
                if(!M.cells.are_simplices()) {
                    M.cells.cell_ptr_[M.cells.nb()] = M.cell_corners.nb();
                }

                //  This warning when loading a single mesh from a file that may
                // contain several meshes -> deactivated for now.	       
                //                if(chunk_class == "SPTR") {
                //                    Logger::out("GeoFile")
                //                        << "File may contain several objects"
                //                        << std::endl;
                //                }

            } catch(const GeoFileException& exc) {
                Logger::err("I/O") << exc.what() << std::endl;
                M.clear();
                return false;
            } catch(...) {
                Logger::err("I/O") << "Caught exception" << std::endl;
                M.clear();                
                return false;
            }


            return true;
        }


        bool load( const std::string& filename, Mesh& M, const MeshIOFlags& ioflags = MeshIOFlags()) override {
            bool result = true;
            try {
                InputGeoFile in(filename);
                result = load(in, M, ioflags);
            }  catch(const GeoFileException& exc) {
                Logger::err("I/O") << exc.what() << std::endl;
                result = false;
            } catch(...) {
                Logger::err("I/O") << "Caught exception" << std::endl;
                result = false;
            }
            return result;
        }

    protected:

        /**
         * \brief Reads an attribute set from a geogram file and
         *  creates the relevant elements in a mesh.
         * \param[in] in a reference to the InputGeoFile
         * \param[in] M a reference to the Mesh
         * \param[in] ioflags the MeshIOFlags that specify which
         *  attributes and mesh elements should be read
         */
        void read_attribute_set(
            InputGeoFile& in,
            Mesh& M,
            const MeshIOFlags& ioflags
        ) {
            const std::string& name =
                in.current_attribute_set().name;
            index_t nb_items = in.current_attribute_set().nb_items;

            if(
                name == "GEO::Mesh::vertices" &&
                ioflags.has_element(MESH_VERTICES)
            ) {
                M.vertices.resize_store(nb_items);
            } else if(
                name == "GEO::Mesh::edges" &&
                ioflags.has_element(MESH_EDGES)
            ) {
                M.edges.resize_store(nb_items);
            } else if(
                name == "GEO::Mesh::facets" &&
                ioflags.has_element(MESH_FACETS)
            ) {
                M.facets.resize_store(nb_items);
            } else if(
                name == "GEO::Mesh::facet_corners" &&
                ioflags.has_element(MESH_FACETS)
            ) {
                M.facet_corners.resize_store(nb_items);
            } else if(
                name == "GEO::Mesh::cells" &&
                ioflags.has_element(MESH_CELLS)
            ) {
                M.cells.resize_store(nb_items);
            } else if(
                name == "GEO::Mesh::cell_corners" &&
                ioflags.has_element(MESH_CELLS)
            ) {
                M.cell_corners.resize_store(nb_items);
            } else if(
                name == "GEO::Mesh::cell_facets" &&
                ioflags.has_element(MESH_CELLS)
            ) {
                M.cell_facets.resize_store(nb_items);
            }
        }

        /**
         * \brief Reads a user attribute from a geogram file and
         *  stores it in a mesh.
         * \param[in] in a reference to the InputGeoFile
         * \param[in] M a reference to the Mesh
         * \param[in] ioflags the MeshIOFlags that specify which
         *  attributes and mesh elements should be read
         */
        void read_user_attribute(
            InputGeoFile& in,
            Mesh& M,
            const MeshIOFlags& ioflags
        ) {
            const std::string& name =
                in.current_attribute().name;
            const std::string& set_name =
                in.current_attribute_set().name;
            if(set_name == "GEO::Mesh::vertices") {
                if(ioflags.has_element(MESH_VERTICES)) {
                    //   Vertex geometry is a special attribute, already
                    // created by the Mesh class, therefore we cannot use
                    // the generic read_attribute() function.
                    if(name == "point") {
                        M.vertices.set_double_precision();
                        M.vertices.set_dimension(
                            in.current_attribute().dimension
                        );
                        in.read_attribute(M.vertices.point_ptr(0));
                    } else if(name == "point_fp32") {
                        M.vertices.set_single_precision();
                        M.vertices.set_dimension(
                            in.current_attribute().dimension
                        );
                        in.read_attribute(
                            M.vertices.single_precision_point_ptr(0)
                        );                                    
                    } else {
                        read_attribute(in, M.vertices.attributes());
                    }
                } 
            } else if(set_name == "GEO::Mesh::edges") {
                if(ioflags.has_element(MESH_EDGES)) {
                    read_attribute(in, M.edges.attributes());
                } 
            } else if(set_name == "GEO::Mesh::facets") {
                if(ioflags.has_element(MESH_FACETS)) {
                    read_attribute(in, M.facets.attributes());
                } 
            } else if(set_name == "GEO::Mesh::facet_corners") {
                if(ioflags.has_element(MESH_FACETS)) {
                    read_attribute(
                        in, M.facet_corners.attributes()
                    );
                } 
            } else if(set_name == "GEO::Mesh::cells") {
                if(ioflags.has_element(MESH_CELLS)) {
                    read_attribute(in, M.cells.attributes());
                } 
            } else if(set_name == "GEO::Mesh::cell_corners") {
                if(ioflags.has_element(MESH_CELLS)) {
                    read_attribute(in, M.cell_corners.attributes());
                } 
            } else if(set_name == "GEO::Mesh::cell_facets") {
                if(ioflags.has_element(MESH_CELLS)) {
                    read_attribute(in, M.cell_facets.attributes());
                } 
            } 
        }

        /**
         * \brief Reads an internal attribute from a geogram file and
         *  stores it in a mesh.
         * \param[in] in a reference to the InputGeoFile
         * \param[in] M a reference to the Mesh
         * \param[in] ioflags the MeshIOFlags that specify which
         *  attributes and mesh elements should be read
         */
        void read_internal_attribute(
            InputGeoFile& in,
            Mesh& M,
            const MeshIOFlags& ioflags
        ) {
            const std::string& name =
                in.current_attribute().name;

            const std::string& set_name =
                in.current_attribute_set().name;

            if(!String::string_starts_with(name, set_name + "::")) {
                Logger::warn("I/O")
                    << "Invalid internal attribute (GEO::Mesh:: scoped): "
                    << name << " does not start with "
                    << set_name
                    << std::endl;
                return;
            }
            if(name == "GEO::Mesh::edges::edge_vertex") {
                if(ioflags.has_element(MESH_EDGES)) {
                    M.edges.edge_vertex_.resize(M.edges.nb()*2);
                    in.read_attribute(M.edges.edge_vertex_.data());
                }
            } else if(name == "GEO::Mesh::facets::facet_ptr") {
                if(ioflags.has_element(MESH_FACETS)) {
                    M.facets.is_simplicial_ = false;
                    M.facets.facet_ptr_.resize(M.facets.nb()+1);
                    in.read_attribute(M.facets.facet_ptr_.data());
                } 
            } else if(name == "GEO::Mesh::facet_corners::corner_vertex") {
                if(ioflags.has_element(MESH_FACETS)) {
                    in.read_attribute(M.facet_corners.corner_vertex_.data());
                } 
            } else if(
                name == "GEO::Mesh::facet_corners::corner_adjacent_facet"
            ) {
                if(ioflags.has_element(MESH_FACETS)) {
                    in.read_attribute(
                        M.facet_corners.corner_adjacent_facet_.data()
                    );
                } 
            } else if(name == "GEO::Mesh::cells::cell_type") {
                if(ioflags.has_element(MESH_CELLS)) {
                    M.cells.is_simplicial_ = false;
                    M.cells.cell_type_.resize(M.cells.nb());
                    in.read_attribute(M.cells.cell_type_.data());
                } 
            } else if(name == "GEO::Mesh::cells::cell_ptr") {
                if(ioflags.has_element(MESH_CELLS)) {
                    M.cells.is_simplicial_ = false;
                    M.cells.cell_ptr_.resize(M.cells.nb()+1);
                    in.read_attribute(M.cells.cell_ptr_.data());
                } 
            } else if(name == "GEO::Mesh::cell_corners::corner_vertex") {
                if(ioflags.has_element(MESH_CELLS)) {
                    in.read_attribute(M.cell_corners.corner_vertex_.data());
                } 
            } else if(name == "GEO::Mesh::cell_facets::adjacent_cell") {
                if(ioflags.has_element(MESH_CELLS)) {
                    in.read_attribute(M.cell_facets.adjacent_cell_.data());
                } 
            } 
        }

        /**
         * \brief Reads a user attribute from a geogram file and 
         *  stores it in an AttributesManager
         * \param[in] in a reference to the InputGeoFile
         * \param[in] attributes a reference to the AttributesManager
         *  where the read attribute should be stored
         */
        void read_attribute(
            InputGeoFile& in,
            AttributesManager& attributes
        ) {
            if(
                !AttributeStore::element_type_name_is_known(
                    in.current_attribute().element_type
                )
            ) {
                Logger::warn("I/O") << "Skipping attribute "
                                    << in.current_attribute().name
                                    << ":"
                                    << in.current_attribute().element_type
                                    << " (unknown type)"
                                    << std::endl;
                return;
            }
            AttributeStore* store =
                AttributeStore::create_attribute_store_by_element_type_name(
                    in.current_attribute().element_type,
                    in.current_attribute().dimension
                );
            attributes.bind_attribute_store(in.current_attribute().name,store);
            in.read_attribute(store->data());
        }

        /**
         * \brief Writes all the user attributes of an AttributesManager
         *  into a geogram file.
         * \param[out] out a reference to the OutputGeoFile
         * \param[in] attribute_set_name the name to be used for the attribute
         *  set in the geogram file
         * \param[in] attributes a reference to the AttributesManager
         */
        void save_attributes(
            OutputGeoFile& out,
            const std::string& attribute_set_name,
            AttributesManager& attributes
        ) {
            vector<std::string> attribute_names;
            attributes.list_attribute_names(attribute_names);
            for(index_t i=0; i<attribute_names.size(); ++i) {
                AttributeStore* store = attributes.find_attribute_store(
                    attribute_names[i]
                );
                if(
                    AttributeStore::element_typeid_name_is_known(
                        store->element_typeid_name()
                    )
                ) {
                    std::string element_type = 
                      AttributeStore::element_type_name_by_element_typeid_name(
                          store->element_typeid_name()
                      );

                    out.write_attribute(
                        attribute_set_name,
                        attribute_names[i],
                        element_type,
                        store->element_size(),
                        store->dimension(),
                        store->data()
                    );
                } else {
                    Logger::warn("I/O")
                        << "Skipping attribute: "
                        << attribute_names[i]
                        << " on "
                        << attribute_set_name
                        << std::endl;
                    Logger::warn("I/O")
                        << "Typeid "
                        << store->element_typeid_name()
                        << " unknown"
                        << std::endl;
                }
            }
        }
    };


#endif
