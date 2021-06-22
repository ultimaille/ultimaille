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
                std::cerr << "Failed to open ASCII file " << filename << std::endl;
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

    void read_vtk_format(const std::string& filename, std::vector<vec3>& verts_, std::vector<int>& edges_, std::vector<int>& tris_, std::vector<int>& quads_, std::vector<int>& tets_, std::vector<int>& hexes_, std::vector<int>& wedges_, std::vector<int>& pyramids_) {
        LineInput li(filename);

        if (!starts_with(li.line, "# vtk DataFile Version")) {
            std::cerr << "This is not a valid VTK file" << std::endl;
            return;
        }

        li.getline(); // any text, 256 characters maximum

        li.getline(); // ASCII or BINARY
        if (!starts_with(li.line, "ASCII")) {
            std::cerr << "Error: only ASCII VTK file format is supported" << std::endl;
            return;
        }

        li.getline_nonempty(); // DATASET, can be one of the following: STRUCTURED_POINTS STRUCTURED_GRID UNSTRUCTURED_GRID POLYDATA STRUCTURED_POINTS RECTILINEAR_GRID FIELD
        if (!starts_with(li.line, "DATASET UNSTRUCTURED_GRID")) {
            std::cerr << "Error: only UNSTRUCTURED_GRID VTK files are supported" << std::endl;
            return;
        }

        { // load the point cloud
            li.getline_nonempty();
            if (li.words.size()!=3 || li.words.front()!="POINTS") {
                std::cerr << "Error: unstructured grid data must begin with POINTS section" << std::endl;
                return;
            }
            if (li.words.back()!="float" && li.words.back()!="double") {
                std::cerr << "Error: unsupported point data type" << std::endl;
                return;
            }
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
            if (li.words.size()!=3 || li.words.front()!="CELLS") {
                std::cerr << "Error: POINTS section must be followed by CELLS section" << std::endl;
                std::cerr << filename << std::endl;
                return;
            }
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

        { // dispatch the connectivity data to different arrays w.r.t the cell type
            li.getline_nonempty();
            if (li.words.size()!=2 || li.words.front()!="CELL_TYPES") {
                std::cerr << "Error: CELLS section must be followed by CELL_TYPES section" << std::endl;
                return;
            }
            int nb_cells = std::stoi(li.words[1]);
            if (nb_cells+1u!=offset.size()) {
                std::cerr << "Error: incoherent number of cells in the CELL_TYPES section" << std::endl;
                return;
            }

            int cnt = 0;
            do {
                li.getline_nonempty();
                for (std::string w : li.words) {
                    int cell_type = std::stoi(w);
                    int cell_size = offset[cnt+1] - offset[cnt];
                    constexpr int pixel2quad[4] = { 0,1,3,2 };
                    constexpr int vtk2geo[8] = { 0,1,3,2,4,5,7,6 };
                    for (int j=0; j<cell_size; j++) {
                        switch (cell_type) {
                            case  3:    edges_.push_back(cells[offset[cnt] + j]);             break;
                            case  5:     tris_.push_back(cells[offset[cnt] + j]);             break;
                            case  8:    quads_.push_back(cells[offset[cnt] + pixel2quad[j]]); break;
                            case  9:    quads_.push_back(cells[offset[cnt] + j]);             break;
                            case 10:     tets_.push_back(cells[offset[cnt] + j]);             break;
                            case 11:    hexes_.push_back(cells[offset[cnt] + j]);             break;
                            case 12:    hexes_.push_back(cells[offset[cnt] + vtk2geo[j]]);    break;
                            case 13:   wedges_.push_back(cells[offset[cnt] + j]);             break;
                            case 14: pyramids_.push_back(cells[offset[cnt] + j]);             break;
                            default: break;
                        }
                    }
                    cnt++;
                }
            } while (li.good() && cnt<nb_cells);
        }

        {
            int attr = -1;
            while (li.good()) {
                li.getline_nonempty();
                if (li.words.size()==2) {
                    if (li.words.front()=="CELL_DATA" ) { attr = 0; li.getline_nonempty(); }
                    if (li.words.front()=="POINT_DATA") { attr = 1; li.getline_nonempty(); }
                }
                if (attr<0) continue;
                if (li.words.size()<3 || li.words.front()!="SCALARS") continue; // TODO NORMALS, VECTORS, TEXTURE_COORDINATES
                std::string name = li.words[1];
                std::string datatype = li.words[2];
                int dim = li.words.size()==4 ? std::stoi(li.words.back()) : 1;
                li.getline_nonempty();
                if (li.words.size()!=2 || li.words.front()!="LOOKUP_TABLE") continue;
                if (li.words.back()!="default")
                    std::cerr << "Error: only default lookup table is supported" << std::endl;

                int nb = attr ? verts_.size() : offset.size()-1;
                int cnt = 0;
                do {
                    li.getline_nonempty();
                    for (std::string w : li.words) {
                        cnt++;
                    }
                } while (li.good() && cnt<nb*dim);
            }
        }
    }

    void write_vtk_format(const std::string& filename, const std::vector<vec3>& verts_, const std::vector<int>& edges_, const std::vector<int>& tris_, const std::vector<int>& quads_, const std::vector<int>& tets_, const  std::vector<int>& hexes_, const  std::vector<int>& wedges_, const  std::vector<int>& pyramids_) {
        std::ofstream out_f;
        out_f.open(filename, std::ifstream::out);
        if (out_f.fail()) {
            std::cerr << "Failed to open " << filename << std::endl;
            return;
        }
        std::stringstream out;
        out << std::setprecision(std::numeric_limits<double>::max_digits10);
        out << "# vtk DataFile Version 3.0" << std::endl;
        out << "Mesh saved with ultimaille: https://github.com/ssloy/ultimaille" << std::endl;
        out << "ASCII" << std::endl;
        out << "DATASET UNSTRUCTURED_GRID" << std::endl;


        out << "POINTS " << verts_.size()  << " double" << std::endl;
        FOR(v, verts_.size()) {
            FOR(d, 3) out << verts_[v][d] << " ";
            out << std::endl;
        }

        out << std::endl;
        int nb_cells = edges_.size() / 2 + tris_.size() / 3 + quads_.size() / 4 + tets_.size() / 4 + hexes_.size() / 8 + wedges_.size()/6 + pyramids_.size()/5;
        int data_size = (edges_.size() / 2) * 3 + (tris_.size() / 3) * 4 + (quads_.size() / 4) * 5 + (tets_.size() / 4) * 5 + (hexes_.size() / 8) * 9 + (wedges_.size() / 6) * 7 + (pyramids_.size() / 5) * 6;
        out << "CELLS " << nb_cells << " " << data_size << std::endl;


        FOR(e, edges_.size() / 2) {
            out << "2";
            FOR(i, 2) out << " " << edges_[2 * e + i];
            out << std::endl;
        }
        FOR(t, tris_.size() / 3) {
            out << "3";
            FOR(i, 3) out << " " << tris_[3 * t + i];
            out << std::endl;
        }
        FOR(q, quads_.size() / 4) {
            out << "4";
            FOR(i, 4) out << " " << quads_[4 * q + i];
            out << std::endl;
        }
        FOR(t, tets_.size() / 4) {
            out << "4";
            FOR(i, 4) out << " " << tets_[4 * t + i];
            out << std::endl;
        }
        FOR(h, hexes_.size() / 8) {
            out << "8";
            constexpr std::array<int, 8> vtk = { 0,1,3,2,4,5,7,6 };
            FOR(i, 8) out << " " << hexes_[8 * h + vtk[i]];
            out << std::endl;
        }
        FOR(t, wedges_.size() / 6) {
            out << "6";
            FOR(i, 6) out << " " << wedges_[6 * t + i];
            out << std::endl;
        }
        FOR(t, pyramids_.size() / 5) {
            out << "5";
            FOR(i, 5) out << " " << pyramids_[5 * t + i];
            out << std::endl;
        }
        out << std::endl;

        out << "CELL_TYPES  " << nb_cells << std::endl;
        FOR(e, edges_.size() / 2) {
            out << "3" << std::endl;
        }
        FOR(t, tris_.size() / 3) {
            out << "5" << std::endl;
        }
        FOR(q, quads_.size() / 4) {
            out << "9" << std::endl;
        }
        FOR(t, tets_.size() / 4) {
            out << "10" << std::endl;
        }
        FOR(h, hexes_.size() / 8) {
            out << "12" << std::endl;
        }
        FOR(h, wedges_.size() / 6) {
            out << "13" << std::endl;
        }
        FOR(h, pyramids_.size() / 5) {
            out << "14" << std::endl;
        }


        out_f << out.rdbuf();
        out_f.close();
    }



    void write_vtk(const std::string filename, const PointSet& ps) {
        std::vector<vec3> verts(ps.size());
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        FOR(v, ps.size()) verts[v] = ps[v];
        write_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
    }

    void write_vtk(const std::string filename, const PolyLine& pl) {
        std::vector<vec3> verts(pl.nverts());
        std::vector<int> edges(2 * pl.nsegments());
        std::vector<int> tris, quads, tets, hexes, wedges, pyramids;
        FOR(v, pl.nverts()) verts[v] = pl.points[v];
        FOR(e, pl.nsegments()) FOR(ev, 2) edges[2 * e + ev] = pl.vert(e, ev);
        write_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
    }

    void write_vtk(const std::string filename, const Surface& m) {
        std::vector<vec3> verts(m.nverts());
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        FOR(v, m.nverts()) verts[v] = m.points[v];
        FOR(f, m.nfacets()) {
            if (m.facet_size(f) == 3) {
                FOR(i, 3) tris.push_back(m.vert(f, i));
            }
            else if (m.facet_size(f) == 4) {
                FOR(i, 4) quads.push_back(m.vert(f, i));
            }
            else {
                std::cerr << "Polygon are not supported in our vtk writer" << std::endl;
            }
        }
        write_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
    }

    void write_vtk(const std::string filename, const Volume& m) {
        std::vector<vec3> verts(m.nverts());;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        FOR(v, m.nverts()) verts[v] = m.points[v];
        if (m.cell_type() == Volume::TETRAHEDRON) {
            tets.resize(4 * m.ncells());
            FOR(t, m.ncells()) FOR(tv, 4) tets[4 * t + tv] = m.vert(t, tv);
        } else if (m.cell_type() == Volume::HEXAHEDRON) {
            hexes.resize(8 * m.ncells());
            FOR(h, m.ncells()) FOR(hv, 8) hexes[8 * h + hv] = m.vert(h, hv);
        } else if (m.cell_type() == Volume::WEDGE) {
            wedges.resize(6 * m.ncells());
            FOR(h, m.ncells()) FOR(hv, 6) wedges[6 * h + hv] = m.vert(h, hv);
        } else if (m.cell_type() == Volume::PYRAMID) {
            pyramids.resize(5 * m.ncells());
            FOR(h, m.ncells()) FOR(hv, 5) pyramids[5 * h + hv] = m.vert(h, hv);
        } else {
            std::cerr << "Volume type : " << m.cell_type() << "; not supported in our vtk writer" << std::endl;
        }

        write_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
    }

    PointSetAttributes read_vtk(const std::string filename, PointSet   &m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        read_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        m = PointSet();
        m.create_points(verts.size());
        FOR(v, verts.size()) m[v] = verts[v];
        return {};
     }

    PolyLineAttributes read_vtk(const std::string filename, PolyLine& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        read_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        m = PolyLine();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_segments(edges.size()/2);
        FOR(e, m.nsegments()) FOR(ev, 2) m.vert(e, ev) = edges[2 * e + ev];
        return {};
    }

    SurfaceAttributes read_vtk(const std::string filename, Triangles& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        read_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        m = Triangles();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_facets(tris.size() / 3);
        FOR(t, m.nfacets()) FOR(tv, 3) m.vert(t, tv) = tris[3 * t + tv];
        return {};
    }

    SurfaceAttributes read_vtk(const std::string filename, Quads& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        read_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        m = Quads();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_facets(quads.size() / 4);
        FOR(q, m.nfacets()) FOR(qv, 4) m.vert(q, qv) = quads[4 * q + qv];
        return {};
    }

    SurfaceAttributes read_vtk(const std::string filename, Polygons& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        read_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        m = Polygons();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_facets(tris.size() / 3, 3);
        FOR(t, m.nfacets()) FOR(tv, 3) m.vert(t, tv) = tris[3 * t + tv];

        int off = m.create_facets(quads.size() / 4, 4);
        FOR(q, quads.size() / 4) FOR(qv, 4) m.vert(off + q, qv) = quads[4 * q + qv];
        return {};
    }

    VolumeAttributes read_vtk(const std::string filename, Tetrahedra& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        read_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        m = Tetrahedra();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_cells(tets.size() / 4);
        FOR(t, m.ncells()) FOR(tv, 4) m.vert(t, tv) = tets[4 * t + tv];
        return {};
    }

    VolumeAttributes read_vtk(const std::string filename, Hexahedra& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        read_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        m = Hexahedra();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_cells(hexes.size() / 8);
        FOR(h, m.ncells()) FOR(hv, 8) m.vert(h, hv) = hexes[8 * h + hv];
        return{};
    }

    VolumeAttributes read_vtk(const std::string filename, Wedges& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        read_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        m = Wedges();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_cells(wedges.size() / 6);
        FOR(h, m.ncells()) FOR(hv, 6) m.vert(h, hv) = wedges[6 * h + hv];
        return{};
    }

    VolumeAttributes read_vtk(const std::string filename, Pyramids& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        read_vtk_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids);
        m = Pyramids();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_cells(pyramids.size() / 5);
        FOR(h, m.ncells()) FOR(hv, 5) m.vert(h, hv) = pyramids[5 * h + hv];
        return{};
    }

}

