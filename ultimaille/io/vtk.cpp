#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <array>
#include "ultimaille/io/vtk.h"

#define FOR(i, n) for(int i = 0; i < static_cast<int>(n); i++)

namespace UM {

    struct LineInput {
        LineInput() = delete;

        LineInput(const std::string& filename) {
            in.open(filename, std::ifstream::in);
            if (in.fail())
                throw std::runtime_error("Failed to open " + filename);
            getline();
        }

        bool good() {
            return in.good();
        }

        void getline() {
            line = {};
            if (in.good()) {
                std::getline(in, line);
                parse_words();
            }
        }

        void getline_nonempty() {
            do {
                getline();
            } while (!words.size() && in.good());
        }

        /*
        int nfields() {
            std::istringstream iss(line);
            return std::distance(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>());
        }
        */

        void parse_words() {
            words.resize(0);
            std::istringstream iss(line);
            std::string str;
            while (iss>>str)
                words.push_back(str);
        }

        ~LineInput() {
            in.close();
        }

        std::ifstream in = {};
        std::string line = {};
        std::vector<std::string> words = {};
    };

    static bool starts_with(const std::string str, const std::string prefix) { // TODO part of C++20, remove after migration
        return ((prefix.size() <= str.size()) && std::equal(prefix.begin(), prefix.end(), str.begin()));
    }

    void read_vtk_format(const std::string& filename, const int celltype2keep, std::vector<vec3>& verts_, std::vector<int> &cells_, std::vector<NamedContainer> attr[2]) {
        LineInput li(filename);

        if (!starts_with(li.line, "# vtk DataFile Version"))
            throw std::runtime_error("This is not a valid VTK file");

        li.getline(); // any text, 256 characters maximum

        li.getline(); // ASCII or BINARY
        if (!starts_with(li.line, "ASCII"))
            throw std::runtime_error("Error: only ASCII VTK file format is supported");

        li.getline_nonempty(); // DATASET, can be one of the following: STRUCTURED_POINTS STRUCTURED_GRID UNSTRUCTURED_GRID POLYDATA STRUCTURED_POINTS RECTILINEAR_GRID FIELD
        if (!starts_with(li.line, "DATASET UNSTRUCTURED_GRID"))
            throw std::runtime_error("Error: only UNSTRUCTURED_GRID VTK files are supported");

        { // load the point cloud
            li.getline_nonempty();
            if (li.words.size()!=3 || li.words.front()!="POINTS")
                throw std::runtime_error("Error: unstructured grid data must begin with POINTS section");
            if (li.words.back()!="float" && li.words.back()!="double")
                throw std::runtime_error("Error: unsupported point data type");
            int nb_verts = std::stoi(li.words[1]);
            verts_.resize(nb_verts);

            int cnt = 0;
            do { // there can be multiple data points per line
                li.getline_nonempty();
                for (std::string w : li.words) {
                    verts_[cnt/3][cnt%3] = std::stod(w);
                    cnt++;
                }
            } while (li.good() && cnt<nb_verts*3);
        }

        std::vector<int> cells = {};
        std::vector<int> offset = { 0 };
        { // load raw connectivity
            li.getline_nonempty();
            if (li.in.eof()) return;
            if (li.words.size()!=3 || li.words.front()!="CELLS")
                throw std::runtime_error("Error: POINTS section must be followed by CELLS section");
            int nb_cells = std::stoi(li.words[1]);
            int data_len = std::stoi(li.words[2]);

            offset.reserve(nb_cells+1);
            cells.reserve(data_len);

            int cell_size = -1; // current cell size
            int cntc = 0; // number of indices loaded for the current cell
            int cntf = 0; // overall number of fields loaded
            while (li.good() && cntf<data_len) {
                li.getline_nonempty();
                for (std::string w : li.words) {
                    if (cell_size<0) {
                        cell_size = std::stoi(w);
                        cntc = 0;
                        offset.push_back(offset.back()+cell_size);
                    } else {
                        cells.push_back(std::stoi(w));
                        if (++cntc==cell_size) cell_size = -1;
                    }
                    cntf++;
                }
            }
        }

        std::vector<int> cell_types;
        { // dispatch the connectivity data to different arrays w.r.t the cell type
            li.getline_nonempty();
            if (li.words.size()!=2 || li.words.front()!="CELL_TYPES")
                throw std::runtime_error("Error: CELLS section must be followed by CELL_TYPES section");
            int nb_cells = std::stoi(li.words[1]);
            if (nb_cells+1u!=offset.size())
                throw std::runtime_error("Error: incoherent number of cells in the CELL_TYPES section");

            int cnt = 0;
            do {
                li.getline_nonempty();
                for (std::string w : li.words) {
                    cell_types.push_back(std::stoi(w));
                    cnt++;
                }
            } while (li.good() && cnt<nb_cells);
            um_assert(nb_cells == (int)cell_types.size());

            for (int i=0; i<nb_cells; i++) {
                if (cell_types[i]!=celltype2keep) continue;
                int cell_size = offset[i+1] - offset[i];
                constexpr int pixel2quad[4] = { 0,1,3,2 };
                constexpr int vtk2geo[8] = { 0,1,3,2,4,5,7,6 };
                for (int j=0; j<cell_size; j++) {
                    switch (cell_types[i]) {
                        case  3:
                        case  5:
                        case  9:
                        case 10:
                        case 11:
                        case 13:
                        case 14: cells_.push_back(cells[offset[i] + j]);             break;
                        case 12: cells_.push_back(cells[offset[i] + vtk2geo[j]]);    break;
                        case  8: cells_.push_back(cells[offset[i] + pixel2quad[j]]); break;
                        default: break;
                    }
                }
            }
        }

        {
            int place = -1;
            while (li.good()) {
                li.getline_nonempty();
                if (li.words.size()==2) {
                    if (li.words.front()=="CELL_DATA" ) { place = 1; li.getline_nonempty(); }
                    if (li.words.front()=="POINT_DATA") { place = 0; li.getline_nonempty(); }
                }
                if (place<0) continue;
//              std::cerr << "==" << place << std::endl;
                if (li.words.size()<3 || li.words.front()!="SCALARS") continue; // TODO NORMALS, VECTORS, TEXTURE_COORDINATES
                std::string name = li.words[1];
                std::string datatype = li.words[2];
                int dim = li.words.size()==4 ? std::stoi(li.words.back()) : 1;
                if (dim!=1)
                    std::cerr << "Warning: multi-component attributes are not supported" << std::endl;

                li.getline_nonempty();
                if (li.words.size()!=2 || li.words.front()!="LOOKUP_TABLE") continue;
                if (li.words.back()!="default")
                    std::cerr << "Warning: only default lookup table is supported" << std::endl;

//                std::cerr << name << std::endl;
                int nb = place ? offset.size()-1 : verts_.size() ;
                int cnt = 0;
                std::vector<std::string> data;
                do {
                    li.getline_nonempty();
                    for (std::string &w : li.words) {
                        if (!place || cell_types[cnt]==celltype2keep)
                            data.push_back(w);
                        cnt++;
                    }
//                  cnt += li.words.size();
//                  data.insert(std::end(data), std::begin(li.words), std::end(li.words));
                } while (li.good() && cnt<nb*dim);
//              um_assert(data.size()==nb*dim);

//              std::cerr << "zzz " << place << " " << datatype << std::endl;
                std::shared_ptr<GenericAttributeContainer> P;
                if (datatype=="bit") {
//                  std::cerr << "Gngaslsdfkgj"<<std::endl;
                    GenericAttribute<bool> A(data.size());
                    std::transform(data.begin(), data.end(), A.ptr->data.begin(), [](const std::string& str) -> bool { return std::stoi(str); });
//                  std::cerr << "a.data.size " << data.size() << " " << A.ptr->data.size() << std::endl;
                    P = A.ptr;
                } else if (datatype=="unsigned_short" || datatype=="short" || datatype=="unsigned_int" || datatype=="int" || datatype=="unsigned_long" || datatype=="long") {
                    GenericAttribute<int> A(data.size());
                    std::transform(data.begin(), data.end(), A.ptr->data.begin(), [](const std::string& str) { return std::stoi(str); });
                    P = A.ptr;
                } else if (datatype=="float" || datatype=="double") {
                    GenericAttribute<double> A(data.size());
                    std::transform(data.begin(), data.end(), A.ptr->data.begin(), [](const std::string& str) { return std::stod(str); });
                    P = A.ptr;
                } else {
                    std::cerr << "Warning: unsupported attribute data type" << std::endl;
                    continue;
                }
                attr[place].emplace_back(name, P);
            }
        }
    }

    void drop_attributes(const std::vector<NamedContainer> &nc, std::ofstream &out) {
        for (const auto &[name, genptr] : nc) {
//          std::cerr << "name " << name << std::endl;
            if (auto p = std::dynamic_pointer_cast<AttributeContainer<int> >(genptr); p.get()!=nullptr) {
                out << "SCALARS " << name << " int 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
                for (auto v : p->data)
                    out << v << " ";
                out << std::endl;
            } else if (auto p = std::dynamic_pointer_cast<AttributeContainer<double> >(genptr); p.get()!=nullptr) {
                out << "SCALARS " << name << " double 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
                for (auto v : p->data)
                    out << v << " ";
                out << std::endl;
            } else if (auto p = std::dynamic_pointer_cast<AttributeContainer<bool> >(genptr); p.get()!=nullptr) {
                out << "SCALARS " << name << " bit 1" << std::endl << "LOOKUP_TABLE default" << std::endl;
                for (auto v : p->data)
                    out << v << " ";
                out << std::endl;
            }
        }
    }

    template <class M, class A> void write_vtk_format(const std::string& filename, const M &m, const A a) {
        std::ofstream out;
        out.open(filename, std::ifstream::out);
        if (out.fail())
            throw std::runtime_error("Failed to open " + filename);
        out << std::setprecision(std::numeric_limits<double>::max_digits10);
        out << "# vtk DataFile Version 3.0" << std::endl;
        out << "Mesh saved with ultimaille: https://github.com/ssloy/ultimaille" << std::endl;
        out << "ASCII" << std::endl;
        out << "DATASET UNSTRUCTURED_GRID" << std::endl;

        if constexpr (std::is_same_v<M, PointSet>) {
            out << "POINTS " << m.size()  << " double" << std::endl;
            for (const vec3 &p : m) {
                for (int d=0; d<3; d++)
                    out << p[d] << " ";
                out << std::endl;
            }
        } else {
            out << "POINTS " << m.nverts()  << " double" << std::endl;
            for (const vec3 &p : m.points) {
                for (int d=0; d<3; d++)
                    out << p[d] << " ";
                out << std::endl;
            }
        }
        out << std::endl;

        if constexpr (std::is_same_v<M, PolyLine>) {
            out << std::endl << "CELLS " << m.nsegments() << " " << (m.nsegments()*(1+2)) << std::endl;
            for (int c=0; c<m.nsegments(); c++) {
                out << "2 ";
                for (int lv=0; lv<2; lv++)
                    out << m.vert(c, lv) << " ";
                out << std::endl;
            }
            out << std::endl << "CELL_TYPES " << m.nsegments() << std::endl;
        } else if constexpr (std::is_base_of_v<Surface, M>) {
            /*
            if (auto ptr = dynamic_cast<Triangles *>(&m)) {
            } else if (auto ptr = dynamic_cast<Quads *>(&m)) {
            } else if (auto ptr = dynamic_cast<Polygons *>(&m)) {
            }
            */
            int nqt = 0;
            int cnt = 0;
            for (int f=0; f<m.nfacets(); f++) {
                if (m.facet_size(f)!=3 && m.facet_size(f)!=4) continue;
                nqt++;
                cnt += 3 + int(m.facet_size(f)==4);
            }
            out << std::endl << "CELLS " << nqt << " " << (nqt+cnt) << std::endl;
            for (int f=0; f<m.nfacets(); f++) {
                if (m.facet_size(f)!=3 && m.facet_size(f)!=4) continue;
                out << m.facet_size(f) << " ";
                for (int lv=0; lv<m.facet_size(f); lv++)
                    out << m.vert(f, lv) << " ";
                out << std::endl;
            }
            out << std::endl << "CELL_TYPES " << nqt << std::endl;
        } /* else if constexpr (std::is_same_v<M, Triangles> || std::is_same_v<M, Quads>) {
            out << std::endl << "CELLS " << m.nfacets() << " " << (m.nfacets()*(1+m.facet_size(0))) << std::endl;
            for (int c=0; c<m.nfacets(); c++) {
                out << m.facet_size(0) << " ";
                for (int lv=0; lv<m.facet_size(0); lv++)
                    out << m.vert(c, lv) << " ";
                out << std::endl;
            }
            out << std::endl << "CELL_TYPES " << m.nfacets() << std::endl;
        } else if constexpr (std::is_same_v<M, Polygons>) {
            int nqt = 0;
            int cnt = 0;
            for (int f=0; f<m.nfacets(); f++) {
                if (m.facet_size(f)!=3 && m.facet_size(f)!=4) continue;
                nqt++;
                cnt += 3+ int(m.facet_size(f)==4);
            }
            out << std::endl << "CELLS " << nqt << " " << (nqt+cnt) << std::endl;
            for (int f=0; f<m.nfacets(); f++) {
                if (m.facet_size(f)!=3 && m.facet_size(f)!=4) continue;
                out << m.facet_size(f) << " ";
                for (int lv=0; lv<m.facet_size(f); lv++)
                    out << m.vert(f, lv) << " ";
                out << std::endl;
            }
            out << std::endl << "CELL_TYPES " << nqt << std::endl;
        }*/
        else if constexpr (std::is_base_of_v<Volume, M>) {
//        else if constexpr (std::is_same_v<M, Tetrahedra> || std::is_same_v<M, Hexahedra> || std::is_same_v<M, Wedges> || std::is_same_v<M, Pyramids>) {
            out << std::endl << "CELLS " << m.ncells() << " " << (m.ncells()*(1+m.nverts_per_cell())) << std::endl;
            for (int c=0; c<m.ncells(); c++) {
                out << m.nverts_per_cell() << " ";
                for (int lv=0; lv<m.nverts_per_cell(); lv++)
                    out << m.vert(c, lv) << " ";
                out << std::endl;
            }
            out << std::endl << "CELL_TYPES " << m.ncells() << std::endl;
        }

        if constexpr (std::is_same_v<M, PolyLine>) {
            for (int c=0; c<m.nsegments(); c++)
                out << "3 ";
            out << std::endl;
        } else if constexpr (std::is_base_of_v<Surface, M>) {
            for (int f=0; f<m.nfacets(); f++) {
                if (m.facet_size(f)==3)
                    out << "5 ";
                if (m.facet_size(f)==4)
                    out << "9 ";
            }
            out << std::endl;
        } else if constexpr (std::is_base_of_v<Volume, M>) {
            for (int c=0; c<m.ncells(); c++)
                switch (m.cell_type) {
                    case Volume::CELL_TYPE::TETRAHEDRON:
                        out << "10 "; break;
// for the CELL_TYPE 12 here is the permutation:
//          constexpr std::array<int, 8> vtk = { 0,1,3,2,4,5,7,6 };
//          FOR(i, 8) out << " " << hexes_[8 * h + vtk[i]];
                    case Volume::CELL_TYPE::HEXAHEDRON:
                        out << "11 "; break;
                    case Volume::CELL_TYPE::WEDGE:
                        out << "13 "; break;
                    case Volume::CELL_TYPE::PYRAMID:
                        out << "14 "; break;
                }
            out << std::endl;
        }

        /* if constexpr (std::is_same_v<M, Triangles>) {
            for (int c=0; c<m.nfacets(); c++)
                out << "5 ";
            out << std::endl;
        } else if constexpr (std::is_same_v<M, Quads>) {
            for (int c=0; c<m.nfacets(); c++)
                out << "9 ";
            out << std::endl;
        } else if constexpr (std::is_same_v<M, Polygons>) {
            for (int f=0; f<m.nfacets(); f++) {
                if (m.facet_size(f)==3)
                    out << "5 ";
                if (m.facet_size(f)==4)
                    out << "9 ";
            }
            out << std::endl;
        } else if constexpr (std::is_same_v<M, Tetrahedra>) {
            for (int c=0; c<m.ncells(); c++)
                out << "10 ";
            out << std::endl;
        } else if constexpr (std::is_same_v<M, Hexahedra>) {
             for (int c=0; c<m.ncells(); c++)
                out << "11 ";
            out << std::endl;
        } else if constexpr (std::is_same_v<M, Wedges>) {
             for (int c=0; c<m.ncells(); c++)
                out << "13 ";
            out << std::endl;
        } else if constexpr (std::is_same_v<M, Pyramids>) {
             for (int c=0; c<m.ncells(); c++)
                out << "14 ";
            out << std::endl;
        }*/

        if constexpr (std::is_same_v<A, PointSetAttributes>) {
        } else if constexpr (std::is_same_v<A, PolyLineAttributes>) {
            if (a.segments.size()) {
                out << std::endl << "CELL_DATA " << m.nsegments() << std::endl;
                drop_attributes(a.segments, out);
            }
        } else if constexpr (std::is_same_v<A, SurfaceAttributes>) {
            if (a.facets.size()) {
                out << std::endl << "CELL_DATA " << m.nfacets() << std::endl;
                drop_attributes(a.facets, out);
            }
        } else if constexpr (std::is_same_v<A, VolumeAttributes>) {
            if (a.cells.size()) {
                out << std::endl << "CELL_DATA " << m.ncells() << std::endl;
                drop_attributes(a.cells, out);
            }
        }

        if (a.points.size()) {
            if constexpr (std::is_same_v<M, PointSet>)
                out << std::endl << "POINT_DATA " << m.size() << std::endl;
            else
                out << std::endl << "POINT_DATA " << m.nverts() << std::endl;
            drop_attributes(a.points, out);
        }

        out.close();
    }

    void write_vtk(const std::string filename, const PointSet& ps, const PointSetAttributes attr) {
        write_vtk_format(filename, ps, attr);
    }

    void write_vtk(const std::string filename, const PolyLine& pl, const PolyLineAttributes attr) {
        write_vtk_format(filename, pl, attr);
    }

    void write_vtk(const std::string filename, const Surface &m, const SurfaceAttributes attr) {
        write_vtk_format(filename, m, attr);
    }

    void write_vtk(const std::string filename, const Volume &m, const VolumeAttributes attr) {
        write_vtk_format(filename, m, attr);
    }


/*
    void write_vtk(const std::string filename, const Triangles& m, const SurfaceAttributes attr) {
        write_vtk_format(filename, m, attr);
    }

    void write_vtk(const std::string filename, const Quads& m, const SurfaceAttributes attr) {
        write_vtk_format(filename, m, attr);
    }

    void write_vtk(const std::string filename, const Polygons& m, const SurfaceAttributes attr) {
        write_vtk_format(filename, m, attr);
    }

    void write_vtk(const std::string filename, const Tetrahedra& m, const VolumeAttributes attr) {
        write_vtk_format(filename, m, attr);
    }

    void write_vtk(const std::string filename, const Hexahedra& m, const VolumeAttributes attr) {
        write_vtk_format(filename, m, attr);
    }

    void write_vtk(const std::string filename, const Wedges& m, const VolumeAttributes attr) {
        write_vtk_format(filename, m, attr);
    }

    void write_vtk(const std::string filename, const Pyramids& m, const VolumeAttributes attr) {
        write_vtk_format(filename, m, attr);
    }
*/

    PointSetAttributes read_vtk(const std::string filename, PointSet   &m) {
        std::vector<vec3> verts;
        std::vector<int> cells;
        std::vector<NamedContainer> attrib[2];
        read_vtk_format(filename, -1, verts, cells, attrib);
        m = PointSet();
        m.create_points(verts.size());
        FOR(v, verts.size()) m[v] = verts[v];
        return { attrib[0] };
     }

    PolyLineAttributes read_vtk(const std::string filename, PolyLine& m) {
        std::vector<vec3> verts;
        std::vector<int> edges;
        std::vector<NamedContainer> attrib[2];
        read_vtk_format(filename, 3, verts, edges, attrib);
        m = PolyLine();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_segments(edges.size()/2);
        FOR(e, m.nsegments()) FOR(ev, 2) m.vert(e, ev) = edges[2 * e + ev];
        return { attrib[0], attrib[1] };
    }

    SurfaceAttributes read_vtk(const std::string filename, Triangles& m) {
        std::vector<vec3> verts;
        std::vector<int> tris;
        std::vector<NamedContainer> attrib[2];
        read_vtk_format(filename, 5, verts, tris, attrib);
        m = Triangles();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_facets(tris.size() / 3);
        FOR(t, m.nfacets()) FOR(tv, 3) m.vert(t, tv) = tris[3 * t + tv];
        return { attrib[0], attrib[1], {} };
    }

    void append_attribute(std::shared_ptr<GenericAttributeContainer> a, std::shared_ptr<GenericAttributeContainer> b) {
        if (auto A = std::dynamic_pointer_cast<AttributeContainer<int> >(a); A.get()!=nullptr) {
            auto B = std::dynamic_pointer_cast<AttributeContainer<int> >(b);
            um_assert(B.get()!=nullptr);
            A->data.insert(std::end(A->data), std::begin(B->data), std::end(B->data));
        } else if (auto A = std::dynamic_pointer_cast<AttributeContainer<double> >(a); A.get()!=nullptr) {
            auto B = std::dynamic_pointer_cast<AttributeContainer<double> >(b);
            um_assert(B.get()!=nullptr);
            A->data.insert(std::end(A->data), std::begin(B->data), std::end(B->data));
        } else if (auto A = std::dynamic_pointer_cast<AttributeContainer<bool> >(a); A.get()!=nullptr) {
            auto B = std::dynamic_pointer_cast<AttributeContainer<bool> >(b);
            um_assert(B.get()!=nullptr);
            A->data.insert(std::end(A->data), std::begin(B->data), std::end(B->data));
        } else {
            std::cerr << "Warning: unsupported attribute type" << std::endl;
        }
    }

    SurfaceAttributes read_vtk(const std::string filename, Quads& m) {
        std::vector<vec3> verts;
        std::vector<int> quads, pixel;
        std::vector<NamedContainer> attrib1[2], attrib2[2];
        read_vtk_format(filename, 9, verts, quads, attrib1);
        read_vtk_format(filename, 8, verts, quads, attrib2);
        um_assert(attrib2[1].size() == attrib1[1].size());

        FOR(i, attrib1[1].size())
            append_attribute(std::get<1>(attrib1[1][i]), std::get<1>(attrib2[1][i]));

        m = Quads();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_facets(quads.size() / 4);
        FOR(q, m.nfacets()) FOR(qv, 4) m.vert(q, qv) = quads[4 * q + qv];

        int off = m.create_facets(pixel.size() / 4);
        FOR(q, pixel.size()/4) FOR(qv, 4) m.vert(off+q, qv) = pixel[4 * q + qv];
        return { attrib1[0], attrib1[1], {} };
    }

    SurfaceAttributes read_vtk(const std::string filename, Polygons& m) {
        std::vector<vec3> verts;
        std::vector<int> tris, quads;
        std::vector<NamedContainer> attrib1[2], attrib2[2];
        read_vtk_format(filename, 5, verts, tris,  attrib1);
        read_vtk_format(filename, 9, verts, quads, attrib2);

        um_assert(attrib2[1].size() == attrib1[1].size());

        FOR(i, attrib1[1].size())
            append_attribute(std::get<1>(attrib1[1][i]), std::get<1>(attrib2[1][i]));

        m = Polygons();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_facets(tris.size() / 3, 3);
        FOR(t, m.nfacets()) FOR(tv, 3) m.vert(t, tv) = tris[3 * t + tv];

        int off = m.create_facets(quads.size() / 4, 4);
        FOR(q, quads.size() / 4) FOR(qv, 4) m.vert(off + q, qv) = quads[4 * q + qv];

        return { attrib1[0], attrib1[1], {} };
    }

    VolumeAttributes read_vtk(const std::string filename, Tetrahedra& m) {
        std::vector<vec3> verts;
        std::vector<int> tetra;
        std::vector<NamedContainer> attrib[2];
        read_vtk_format(filename, 10, verts, tetra, attrib);
        m = Tetrahedra();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_cells(tetra.size() / 4);
        FOR(t, m.ncells()) FOR(tv, 4) m.vert(t, tv) = tetra[4 * t + tv];
        return { attrib[0], attrib[1], {}, {} };
    }

    VolumeAttributes read_vtk(const std::string filename, Hexahedra& m) {
        std::vector<vec3> verts;
        std::vector<int> hexa, voxel;
        std::vector<NamedContainer> attrib1[2], attrib2[2];
        read_vtk_format(filename, 11, verts, hexa,  attrib1);
        read_vtk_format(filename, 12, verts, voxel, attrib2);

        um_assert(attrib2[1].size() == attrib1[1].size());

        FOR(i, attrib1[1].size())
            append_attribute(std::get<1>(attrib1[1][i]), std::get<1>(attrib2[1][i]));

        m = Hexahedra();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_cells(hexa.size() / 8);
        FOR(h, m.ncells()) FOR(hv, 8) m.vert(h, hv) = hexa[8 * h + hv];
        int off = m.create_cells(voxel.size() / 8);
        FOR(h, voxel.size()/8) FOR(hv, 8) m.vert(off+h, hv) = voxel[8 * h + hv];
        return { attrib1[0], attrib1[1], {}, {} };
    }

    VolumeAttributes read_vtk(const std::string filename, Wedges& m) {
        std::vector<vec3> verts;
        std::vector<int> wedges;
        std::vector<NamedContainer> attrib[2];
        read_vtk_format(filename, 13, verts, wedges, attrib);
        m = Wedges();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_cells(wedges.size() / 6);
        FOR(h, m.ncells()) FOR(hv, 6) m.vert(h, hv) = wedges[6 * h + hv];
        return { attrib[0], attrib[1], {}, {} };
    }

    VolumeAttributes read_vtk(const std::string filename, Pyramids& m) {
        std::vector<vec3> verts;
        std::vector<int> pyramids;
        std::vector<NamedContainer> attrib[2];
        read_vtk_format(filename, 14, verts, pyramids, attrib);
        m = Pyramids();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_cells(pyramids.size() / 5);
        FOR(h, m.ncells()) FOR(hv, 5) m.vert(h, hv) = pyramids[5 * h + hv];
        return { attrib[0], attrib[1], {}, {} };
    }
}

