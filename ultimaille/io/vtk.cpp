#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <array>
#include "ultimaille/io/vtk.h"

#define FOR(i, n) for(int i = 0; i < static_cast<int>(n); i++)

namespace UM {
    inline void file_must_no_be_at_end(std::ifstream& f, const std::string& reason = " should'nt") {
        if (f.eof()) {
            f.close();
            std::cout << "File ended to soon while : " << reason << std::endl;
            exit(1);
        }
    }

    inline static bool string_start(const std::string& string, const std::string& start_of_string) {
        size_t start = 0;
        FOR(i, string.size()) if (string[i] != ' ' && string[i] != '\t') {
            start = (size_t)i;
            break;
        }
        std::string copy_without_space(string.begin() + start, string.end());
        if (copy_without_space.size() < start_of_string.size()) return false;
        return (std::string(copy_without_space.begin(), copy_without_space.begin() + (long int)start_of_string.size()) == start_of_string);
    }

    void read_vtk_format(const std::string& filename, std::vector<vec3>& verts_, std::vector<int>& edges_, std::vector<int>& tris_, std::vector<int>& quads_, std::vector<int>& tets_, std::vector<int>& hexes_, std::vector<int>& wedges_, std::vector<int>& pyramids_) {
        std::ifstream in;
        in.open(filename, std::ifstream::in);
        if (in.fail()) {
            std::cerr << "Failed to open " << filename << std::endl;
            exit(1);
        }

        std::string line;
        std::getline(in, line);
        if (!string_start(line, "# vtk DataFile Version")) {
            std::cerr << "This vtk files seems incomplete, still trying... " << std::endl;
        }
        std::getline(in, line);
//      std::cerr << "File loaded info :  " << line << std::endl;
        std::getline(in, line);
        if (line.find("ASCII") == std::string::npos) {
            std::cerr << "This reader only reads standard ASCII VTK, failed reading " << filename << std::endl;
            exit(1);
        }
        std::getline(in, line);
        if (line.find("UNSTRUCTURED_GRID") == std::string::npos) {
            std::cerr << "This reader only reads UNSTRUCTURED_GRID standard ASCII VTK, failed reading " << filename << std::endl;
        }
        in >> line;
        if (!string_start(line, "POINTS")) {
            std::cerr << "This reader only reads UNSTRUCTURED_GRID standard ASCII VTK, failed reading " << filename << std::endl;
            exit(1);
        }
        int nb_vertices = 0; in >> nb_vertices;
        std::string uselessstring; in >> uselessstring;
        verts_.resize(nb_vertices);
        FOR(v, nb_vertices) FOR(d, 3) {
            file_must_no_be_at_end(in, "parsing VTK");
            in >> verts_[v][d];
        }
        std::getline(in, line);
        file_must_no_be_at_end(in, "parsing VTK");
        in >> line;
        if (!string_start(line, "CELLS")) {
            std::cerr << "Error in vtk while reading : " << filename << std::endl;
            exit(1);
        }
        file_must_no_be_at_end(in, "parsing VTK");
        int nb_cells = 0; in >> nb_cells;
        file_must_no_be_at_end(in, "parsing VTK");
        int nb_data = 0; in >> nb_data;
        std::vector<int> start_of_cell(nb_cells);
        std::vector<int> cell_content(nb_data - nb_cells);
        int actual_start = 0;
        FOR(i, nb_cells) {
            file_must_no_be_at_end(in, "parsing VTK");
            int cell_size = 0; in >> cell_size;
            start_of_cell[i] = actual_start;
            FOR(j, cell_size) {
                file_must_no_be_at_end(in, "parsing VTK");
                int id = 0; in >> id;
                cell_content[actual_start + j] = id;
            }
            actual_start += cell_size;
        }
        std::getline(in, line);
        file_must_no_be_at_end(in, "parsing VTK");
        in >> line;
        if (!string_start(line, "CELL_TYPES")) {
            std::cerr << "Error in vtk while reading : " << filename << std::endl;
            exit(1);
        }
        file_must_no_be_at_end(in, "parsing VTK");
        int new_nb_cells = 0; in >> new_nb_cells;
        if (new_nb_cells != nb_cells) {
            std::cerr << "Incoherent nb of cell : " << filename << std::endl;
            exit(1);
        }
        constexpr int pixel2quad[4] = { 0,1,3,2 };
        constexpr int vtk2geo[8] = { 0,1,3,2,4,5,7,6 };

        FOR(i, nb_cells) {
            file_must_no_be_at_end(in, "parsing VTK");
            int cell_type = 0; in >> cell_type;
            switch (cell_type)
            {
            case 3:
                FOR(j, 2) edges_.push_back(cell_content[start_of_cell[i] + j]);
                break;
            case 5:
                FOR(j, 3) tris_.push_back(cell_content[start_of_cell[i] + j]);
                break;
            case 8:
                FOR(j, 3) quads_.push_back(cell_content[start_of_cell[i] + pixel2quad[j]]);
                break;
            case 9:
                FOR(j, 4) quads_.push_back(cell_content[start_of_cell[i] + j]);
                break;
            case 10:
                FOR(j, 4) tets_.push_back(cell_content[start_of_cell[i] + j]);
                break;
            case 11:
                FOR(j, 8) hexes_.push_back(cell_content[start_of_cell[i] + j]);
                break;
            case 12:
                FOR(j, 8) hexes_.push_back(cell_content[start_of_cell[i] + vtk2geo[j]]);
                break;
            case 13:
                FOR(j, 6) wedges_.push_back(cell_content[start_of_cell[i] + j]);
                break;
            case 14:
                FOR(j, 5) pyramids_.push_back(cell_content[start_of_cell[i] + j]);
                break;
            default:
                break;
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
        out << std::fixed << std::setprecision(4);
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

