#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <array>
#include "ultimaille/io/medit.h"
#define FOR(i, n) for(int i = 0; i < n; i++)

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
    void read_medit_format(const std::string& filename, std::vector<vec3>& verts_, std::vector<int>& edges_, std::vector<int>& tris_, std::vector<int>& quads_, std::vector<int>& tets_, std::vector<int>& hexes_) {

        std::ifstream in;
        in.open(filename, std::ifstream::in);
        if (in.fail()) {
            std::cerr << "Failed to open " << filename << std::endl;
            return;
        }

        std::string firstline;

        while (!in.eof()) {
            std::getline(in, firstline);
            if (string_start(firstline, "Vertices")) {
                std::string line;
                int nb_of_vertices = 0;
                {
                    file_must_no_be_at_end(in, "parsing vertices");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    iss >> nb_of_vertices;
                }
                verts_.resize(nb_of_vertices);
                FOR(v, nb_of_vertices) {
                    file_must_no_be_at_end(in, "parsing vertices");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    FOR(i, 3)  iss >> verts_[v][i];
                }
            }
            if (string_start(firstline, "Edges")) {
                std::string line;
                int nb_of_edges = 0;
                {
                    file_must_no_be_at_end(in, "parsing Edges");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    iss >> nb_of_edges;
                }
                edges_.resize(2 * nb_of_edges);
                FOR(e, nb_of_edges) {
                    file_must_no_be_at_end(in, "parsing Edges");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    FOR(i, 2) {
                        int a = 0;
                        iss >> a;
                        edges_[2 * e + i] = a - 1;
                    }
                }
            }
            if (string_start(firstline, "Triangles")) {
                std::string line;
                int nb_of_tri = 0;
                {
                    file_must_no_be_at_end(in, "parsing Triangles");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    iss >> nb_of_tri;
                }
                tris_.resize(3 * nb_of_tri);
                FOR(t, nb_of_tri) {
                    file_must_no_be_at_end(in, "parsing Triangles");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    FOR(i, 3) {
                        int a = 0;
                        iss >> a;
                        tris_[3 * t + i] = a - 1;
                    }
                }
            }
            if (string_start(firstline, "Quadrilaterals")) {
                std::string line;
                int nb_of_quads = 0;
                {
                    file_must_no_be_at_end(in, "parsing Quadrilaterals");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    iss >> nb_of_quads;
                }
                quads_.resize(4 * nb_of_quads);
                FOR(q, nb_of_quads) {
                    file_must_no_be_at_end(in, "parsing Quadrilaterals");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    FOR(i, 4) {
                        int a = 0;
                        iss >> a;
                        quads_[4 * q + i] = a - 1;
                    }
                }
            }
            if (string_start(firstline, "Tetrahedra")) {
                std::string line;
                int nb_of_tets = 0;
                {
                    file_must_no_be_at_end(in, "parsing Tetrahedra");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    iss >> nb_of_tets;
                }
                tets_.resize(4 * nb_of_tets);
                FOR(t, nb_of_tets) {
                    file_must_no_be_at_end(in, "parsing Tetrahedra");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    FOR(i, 4) {
                        int a = 0;
                        iss >> a;
                        tets_[4 * t + i] = a - 1;
                    }
                }
            }
            if (string_start(firstline, "Hexahedra")) {
                std::string line;
                int nb_of_hexs = 0;
                {
                    file_must_no_be_at_end(in, "parsing Hexahedra");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    iss >> nb_of_hexs;
                }
                hexes_.resize(8 * nb_of_hexs);
                FOR(h, nb_of_hexs) {
                    file_must_no_be_at_end(in, "parsing Hexahedra");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    FOR(i, 8) {
                        int a = 0;
                        iss >> a;
                        constexpr std::array<int, 8> medit = { 1,0,2,3,5,4,6,7 };

                        hexes_[8 * h + medit[i]] = a - 1;
                    }
                }
            }

        }

    }


    void write_medit_format(const std::string& filename, const std::vector<vec3>& verts_, const std::vector<int>& edges_, const std::vector<int>& tris_, const std::vector<int>& quads_, const std::vector<int>& tets_, const  std::vector<int>& hexes_, const bool GMSH_numerotation) {
        std::ofstream out_f;
        out_f.open(filename, std::ifstream::out);
        if (out_f.fail()) {
            std::cerr << "Failed to open " << filename << std::endl;
            return;
        }
        std::stringstream out;
        out << std::fixed << std::setprecision(4);
        out << "MeshVersionFormatted 2" << std::endl << std::endl;
        out << "Dimension" << std::endl << "3" << std::endl << std::endl;

        out << "Vertices" << std::endl;
        out << verts_.size() << std::endl;
        FOR(v, verts_.size()) {
            FOR(d, 3) out << verts_[v][d] << " ";
            out << "1" << std::endl;
        }
        out << std::endl;
        if (edges_.size() > 0) {
            out << "Edges" << std::endl;
            out << edges_.size() / 2 << std::endl;
            FOR(e, edges_.size() / 2) {
                FOR(i, 2) out << edges_[2 * e + i] + 1 << " ";
                out << "1" << std::endl;
            }
            out << "end" << std::endl;
        }
        if (tris_.size() > 0) {
            out << "Triangles" << std::endl;
            out << tris_.size() / 3 << std::endl;
            FOR(t, tris_.size() / 3) {
                FOR(i, 3) out << tris_[3 * t + i] + 1 << " ";
                out << "1" << std::endl;
            }
            out << "end" << std::endl;
        }
        if (quads_.size() > 0) {
            out << "Quadrilaterals" << std::endl;
            out << quads_.size() / 4 << std::endl;
            FOR(q, quads_.size() / 4) {
                FOR(i, 4) out << quads_[4 * q + i] + 1 << " ";
                out << "1" << std::endl;
            }
            out << "end" << std::endl;
        }
        if (tets_.size() > 0) {
            out << "Tetrahedra" << std::endl;
            out << tets_.size() / 4 << std::endl;
            FOR(t, tets_.size() / 4) {
                FOR(i, 4) out << tets_[4*t +i] + 1 << " ";
                out << "1" << std::endl;
            }
            out << "end" << std::endl;
        }

        if (hexes_.size() > 0) {
            out << "Hexahedra" << std::endl;
            out << hexes_.size() / 8 << std::endl;
            FOR(h, hexes_.size() / 8) {
                // geogram convention -> GMSH + Medit
                constexpr std::array<int, 8> sign_of_det = { 0,2,1,3,4,6,5,7 };
                constexpr std::array<int, 8> medit = { 1,0,2,3,5,4,6,7 };
                if (GMSH_numerotation)
                    FOR(i, 8) out << hexes_[8 * h + sign_of_det[medit[i]]] + 1 << " ";
                else
                    FOR(i, 8) out << hexes_[8 * h + medit[i]] + 1 << " ";

                out << "1" << std::endl;
            }
            out << "end" << std::endl;
        }


        out_f << out.rdbuf();
        out_f.close();
    }



    void write_medit(const std::string filename, const PolyLine& pl) {
        std::vector<vec3> verts(pl.nverts());
        std::vector<int> edges(2 * pl.nsegments());
        std::vector<int> tris;
        std::vector<int> quads;
        std::vector<int> tets;
        std::vector<int> hexes;
        FOR(v, pl.nverts()) verts[v] = pl.points[v];
        FOR(e, pl.nsegments()) FOR(ev, 2) edges[2 * e + ev] = pl.vert(e, ev);
        write_medit_format(filename, verts, edges, tris, quads, tets, hexes, false);
    }
    void writ_medit(const std::string filename, const Surface& m) {
        std::vector<vec3> verts(m.nverts());
        std::vector<int> edges;
        std::vector<int> tris;
        std::vector<int> quads;
        std::vector<int> tets;
        std::vector<int> hexes;
        FOR(v, m.nverts()) verts[v] = m.points[v];
        FOR(f, m.nfacets()) {
            if (m.facet_size(f) == 3) {
                FOR(i, 3) tris.push_back(m.vert(f, i));
            }
            else if (m.facet_size(f) == 4) {
                FOR(i, 4) quads.push_back(m.vert(f, i));
            }
            else {
                std::cerr << "Polygon are not supported in our MEDIT writer";
            }
        }
        write_medit_format(filename, verts, edges, tris, quads, tets, hexes, false);
    }
    void write_medit(const std::string filename, const Volume& m, bool hexes_GMSH_numerotation) {
        std::vector<vec3> verts(m.nverts());;
        std::vector<int> edges;
        std::vector<int> tris;
        std::vector<int> quads;
        std::vector<int> tets;
        std::vector<int> hexes;
        FOR(v, m.nverts()) verts[v] = m.points[v];
        if (m.cell_type() == 0) {
            tets.resize(4 * m.ncells());
            FOR(t, m.ncells()) FOR(tv, 4) tets[4 * t + tv] = m.vert(t, tv);
        }
        else if (m.cell_type() == 1) {
            hexes.resize(8 * m.ncells());
            FOR(h, m.ncells()) FOR(hv, 8) hexes[8 * h + hv] = m.vert(h, hv);
        }
        else {
            std::cerr << "Volume type : " << m.cell_type() << "; not supported in our MEDIT writer";
        }

        write_medit_format(filename, verts, edges, tris, quads, tets, hexes, hexes_GMSH_numerotation);
    }


    void read_medit(const std::string filename, PolyLine& m) {
        std::vector<vec3> verts;
        std::vector<int> edges;
        std::vector<int> tris;
        std::vector<int> quads;
        std::vector<int> tets;
        std::vector<int> hexes;
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes);
        m = PolyLine();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_segments(edges.size()/2);
        FOR(e, m.nsegments()) FOR(ev, 2) m.vert(e, ev) = edges[2 * e + ev];

    }
    void read_medit(const std::string filename, Triangles& m) {
        std::vector<vec3> verts;
        std::vector<int> edges;
        std::vector<int> tris;
        std::vector<int> quads;
        std::vector<int> tets;
        std::vector<int> hexes;
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes);
        m = Triangles();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_facets(tris.size() / 3);
        FOR(t, m.nfacets()) FOR(tv, 3) m.vert(t, tv) = tris[3 * t + tv];
    }
    void read_medit(const std::string filename, Quads& m) {
        std::vector<vec3> verts;
        std::vector<int> edges;
        std::vector<int> tris;
        std::vector<int> quads;
        std::vector<int> tets;
        std::vector<int> hexes;
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes);
        m = Quads();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_facets(quads.size() / 4);
        FOR(q, m.nfacets()) FOR(qv, 4) m.vert(q, qv) = quads[4 * q + qv];
    }
    void read_medit(const std::string filename, Polygons& m) {
        std::vector<vec3> verts;
        std::vector<int> edges;
        std::vector<int> tris;
        std::vector<int> quads;
        std::vector<int> tets;
        std::vector<int> hexes;
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes);
        m = Polygons();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_facets(tris.size() / 3, 3);
        FOR(t, m.nfacets()) FOR(tv, 3) m.vert(t, tv) = tris[3 * t + tv];

        int off = m.create_facets(quads.size() / 4, 4);
        FOR(q, quads.size() / 4) FOR(qv, 4) m.vert(off + q, qv) = quads[4 * q + qv];

    }
    void read_medit(const std::string filename, Tetrahedra& m) {
        std::vector<vec3> verts;
        std::vector<int> edges;
        std::vector<int> tris;
        std::vector<int> quads;
        std::vector<int> tets;
        std::vector<int> hexes;
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes);
        m = Tetrahedra();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_cells(tets.size() / 4);
        FOR(t, m.ncells()) FOR(tv, 4) m.vert(t, tv) = tets[4 * t + tv];
    }
    void read_medit(const std::string filename, Hexahedra& m) {
        std::vector<vec3> verts;
        std::vector<int> edges;
        std::vector<int> tris;
        std::vector<int> quads;
        std::vector<int> tets;
        std::vector<int> hexes;
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes);
        m = Hexahedra();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_cells(tets.size() / 8);
        FOR(h, m.ncells()) FOR(hv, 8) m.vert(h, hv) = hexes[8 * h + hv];
    }


}

