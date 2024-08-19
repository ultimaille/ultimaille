#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <array>
#include "ultimaille/io/medit.h"
#include "ultimaille/syntactic-sugar/assert.h"
#define FOR(i, n) for(int i = 0; i < static_cast<int>(n); i++)

constexpr std::array<int, 8> hex_medit2geogram = { 0,1,3,2,4,5,7,6 };

namespace UM {
    inline void file_must_no_be_at_end(std::ifstream& f, const std::string& reason = " should not") {
        if (f.eof()) {
            f.close();
            throw std::runtime_error("File ended to soon while " + reason);
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

    const std::string attrib_set_names[4] = {"GEO::Mesh::vertices", "GEO::Mesh::edges", "GEO::Mesh::facets", "GEO::Mesh::cells"};

    void read_medit_format(const std::string& filename, std::vector<vec3>& verts_, std::vector<int>& edges_, std::vector<int>& tris_, std::vector<int>& quads_, std::vector<int>& tets_, std::vector<int>& hexes_, std::vector<int>& wedges_, std::vector<int>& pyramids_, std::vector<NamedContainer> attr[4]) {
        std::ifstream in;
        in.open(filename, std::ifstream::in);
        if (in.fail())
            throw std::runtime_error("Failed to open " + filename);

        const char attr_name[] = "region";

        std::string firstline;

        int nb_of_tri = 0, nb_of_quads = 0;
        std::vector<int> facet_region;

        int nb_of_tets = 0, nb_of_hexs = 0, nb_of_prisms = 0, nb_of_pyramids = 0;
        std::vector<int> cell_region;

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

                // Create vertices region attribute
                GenericAttribute<int> vertices_region(nb_of_vertices);

                verts_.resize(nb_of_vertices);
                FOR(v, nb_of_vertices) {
                    file_must_no_be_at_end(in, "parsing vertices");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    FOR(i, 3)  iss >> verts_[v][i];
                    iss >> vertices_region[v];
                }
                attr[0].emplace_back(attr_name, vertices_region.ptr);
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

                // Create edges region attribute
                GenericAttribute<int> edges_region(nb_of_edges);

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
                    iss >> edges_region[e];
                }
                attr[1].emplace_back(attr_name, edges_region.ptr);
            }           

            if (string_start(firstline, "Triangles")) {
                std::string line;
                {
                    file_must_no_be_at_end(in, "parsing Triangles");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    iss >> nb_of_tri;
                }
                
                facet_region.resize(nb_of_tri + nb_of_quads);

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
                    iss >> facet_region[t];
                }
            }
            if (string_start(firstline, "Quadrilaterals")) {
                std::string line;
                {
                    file_must_no_be_at_end(in, "parsing Quadrilaterals");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    iss >> nb_of_quads;
                }

                facet_region.resize(nb_of_tri + nb_of_quads);

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
                    iss >> facet_region[nb_of_tri + q];
                }
            }

            if (string_start(firstline, "Tetrahedra")) {
                std::string line;
                {
                    file_must_no_be_at_end(in, "parsing Tetrahedra");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    iss >> nb_of_tets;
                }

                cell_region.resize(nb_of_tets + nb_of_hexs + nb_of_prisms + nb_of_pyramids);

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
                    iss >> cell_region[t];
                }
            }
            if (string_start(firstline, "Hexahedra")) {
                std::string line;
                {
                    file_must_no_be_at_end(in, "parsing Hexahedra");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    iss >> nb_of_hexs;
                }

                cell_region.resize(nb_of_tets + nb_of_hexs + nb_of_prisms + nb_of_pyramids);
                
                hexes_.resize(8 * nb_of_hexs);
                FOR(h, nb_of_hexs) {
                    file_must_no_be_at_end(in, "parsing Hexahedra");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    FOR(i, 8) {
                        int a = 0;
                        iss >> a;
                        hexes_[8 * h + hex_medit2geogram[i]] = a - 1;
                    }
                    iss >> cell_region[nb_of_tets + h];
                }
            }
            if (string_start(firstline, "Prisms")) {
                std::string line;
                {
                    file_must_no_be_at_end(in, "parsing Prisms");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    iss >> nb_of_prisms;
                }

                cell_region.resize(nb_of_tets + nb_of_hexs + nb_of_prisms + nb_of_pyramids);

                wedges_.resize(6 * nb_of_prisms);
                FOR(h, nb_of_prisms) {
                    file_must_no_be_at_end(in, "parsing Prisms");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    FOR(i, 6) {
                        int a = 0;
                        iss >> a;
                        wedges_[6 * h + i] = a - 1;
                    }
                    iss >> cell_region[nb_of_tets + nb_of_hexs + h];
                }
            }
            if (string_start(firstline, "Pyramids")) {
                std::string line;
                {
                    file_must_no_be_at_end(in, "parsing Pyramids");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    iss >> nb_of_pyramids;
                }

                cell_region.resize(nb_of_tets + nb_of_hexs + nb_of_prisms + nb_of_pyramids);

                pyramids_.resize(5 * nb_of_pyramids);
                FOR(h, nb_of_pyramids) {
                    file_must_no_be_at_end(in, "parsing Pyramids");
                    std::getline(in, line);
                    std::istringstream iss(line.c_str());
                    FOR(i, 5) {
                        int a = 0;
                        iss >> a;
                        pyramids_[5 * h + i] = a - 1;
                    }
                    iss >> cell_region[nb_of_tets + nb_of_hexs + nb_of_prisms + h];
                }
            }
        }

        // Create facet region attribute
        const int nb_facets = nb_of_tri + nb_of_quads;
        if (nb_facets > 0) {    
            GenericAttribute<int> facet_region_attr(nb_facets);
            for (size_t i = 0; i < facet_region.size(); i++) {
                facet_region_attr[i] = facet_region[i];
            }

            attr[2].emplace_back(attr_name, facet_region_attr.ptr);
        }
        // Create cell region attribute
        const int nb_cells = nb_of_tets + nb_of_hexs + nb_of_prisms + nb_of_pyramids;
        if (nb_cells > 0) {    
            GenericAttribute<int> cell_region_attr(nb_cells);
            for (size_t i = 0; i < cell_region.size(); i++) {
                cell_region_attr[i] = cell_region[i];
            }

            attr[3].emplace_back(attr_name, cell_region_attr.ptr);
        }

    }

    std::vector<int> load_attr(std::vector<NamedContainer> &container, int size) {
        
        if (container.size() <= 0) {
            return std::move(std::vector<int>(size, 1));
        }

        std::shared_ptr<GenericAttributeContainer> point_attr_ptr = container[0].second;
        if (auto point_cont_ptr = std::dynamic_pointer_cast<AttributeContainer<int>>(point_attr_ptr); point_cont_ptr.get()!=nullptr) {
            return std::move(point_cont_ptr->data);
        } else {
            assert(false);
        }
    }

    void write_medit_format(const std::string& filename, const std::vector<vec3>& verts_, const std::vector<int>& edges_, const std::vector<int>& tris_, const std::vector<int>& quads_, const std::vector<int>& tets_, const  std::vector<int>& hexes_, const std::vector<int>& wedges_, const  std::vector<int>& pyramids_, std::vector<NamedContainer> attr[4]) {
        std::ofstream out_f;
        out_f.open(filename, std::ifstream::out);
        if (out_f.fail())
            throw std::runtime_error("Failed to open " + filename);
        std::stringstream out;
        out << std::setprecision(std::numeric_limits<double>::max_digits10);
        out << "MeshVersionFormatted 2" << std::endl << std::endl;
        out << "Dimension" << std::endl << "3" << std::endl << std::endl;

        out << "Vertices" << std::endl;
        out << verts_.size() << std::endl;
        auto point_attr = load_attr(attr[0], verts_.size());
        FOR(v, verts_.size()) {
            FOR(d, 3) out << verts_[v][d] << " ";
            out << point_attr[v] << std::endl;
        }
        out << std::endl;
        if (edges_.size() > 0) {
            out << "Edges" << std::endl;
            out << edges_.size() / 2 << std::endl;
            auto edge_attr = load_attr(attr[1], edges_.size());
            FOR(e, edges_.size() / 2) {
                FOR(i, 2) out << edges_[2 * e + i] + 1 << " ";
                out << edge_attr[e] << std::endl;
            }
            out << std::endl << "End" << std::endl;
        }
        
        // Load facet attributes
        auto facet_attr = load_attr(attr[2], tris_.size() / 3 + quads_.size() / 4);

        if (tris_.size() > 0) {
            out << "Triangles" << std::endl;
            out << tris_.size() / 3 << std::endl;
            FOR(t, tris_.size() / 3) {
                FOR(i, 3) out << tris_[3 * t + i] + 1 << " ";
                out << facet_attr[t] << std::endl;
            }
            out << std::endl << "End" << std::endl;
        }
        if (quads_.size() > 0) {
            out << "Quadrilaterals" << std::endl;
            out << quads_.size() / 4 << std::endl;
            FOR(q, quads_.size() / 4) {
                FOR(i, 4) out << quads_[4 * q + i] + 1 << " ";
                out << facet_attr[tris_.size() / 3 + q] << std::endl;
            }
            out << std::endl << "End" << std::endl;
        }

        // Load cell attributes
        auto cell_attr = load_attr(attr[3], tets_.size() / 4 + hexes_.size() / 8 + wedges_.size() / 6 + pyramids_.size() / 5);

        if (tets_.size() > 0) {
            out << "Tetrahedra" << std::endl;
            out << tets_.size() / 4 << std::endl;
            FOR(t, tets_.size() / 4) {
                FOR(i, 4) out << tets_[4*t +i] + 1 << " ";
                out << cell_attr[t] << std::endl;
            }
            out << std::endl << "End" << std::endl;
        }
        if (hexes_.size() > 0) {
            out << "Hexahedra" << std::endl;
            out << hexes_.size() / 8 << std::endl;
            FOR(h, hexes_.size() / 8) {
                FOR(i, 8) out << hexes_[8 * h + hex_medit2geogram[i]] + 1 << " ";
                out << cell_attr[tets_.size() / 4 + h] << std::endl;
            }
            out << std::endl << "End" << std::endl;
        }
        if (wedges_.size() > 0) {
            out << "Prisms" << std::endl;
            out << wedges_.size() / 6 << std::endl;
            FOR(h, wedges_.size() / 6) {
                FOR(i, 6) out << wedges_[6 * h + i] + 1 << " ";
                out << cell_attr[tets_.size() / 4 + hexes_.size() / 8 + h] << std::endl;
            }
            out << std::endl << "End" << std::endl;
        }
        if (pyramids_.size() > 0) {
            out << "Pyramids" << std::endl;
            out << pyramids_.size() / 5 << std::endl;
            FOR(h, pyramids_.size() / 5) {
                FOR(i, 5) out << pyramids_[5 * h + i] + 1 << " ";
                out << cell_attr[tets_.size() / 4 + hexes_.size() / 8 + wedges_.size() / 6 + h] << std::endl;
            }
            out << std::endl << "End" << std::endl;
        }

        out_f << out.rdbuf();
        out_f.close();
    }

    void write_medit(const std::string filename, const PointSet& ps, PointSetAttributes attr) {
        std::vector<vec3> verts(ps.size());
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        FOR(v, ps.size()) verts[v] = ps[v];
        std::vector<NamedContainer> A[4] = {attr.points, {}, {}, {}};
        write_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, A);
    }

    void write_medit(const std::string filename, const PolyLine& pl, PolyLineAttributes attr) {
        std::vector<vec3> verts(pl.nverts());
        std::vector<int> edges(2 * pl.nedges());
        std::vector<int> tris, quads, tets, hexes, wedges, pyramids;
        FOR(v, pl.nverts()) verts[v] = pl.points[v];
        FOR(e, pl.nedges()) FOR(ev, 2) edges[2 * e + ev] = pl.vert(e, ev);
        std::vector<NamedContainer> A[4] = {attr.points, attr.edges, {}, {}};
        write_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, A);
    }

    void write_medit(const std::string filename, const Surface& m, SurfaceAttributes attr) {
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
            else
                std::cerr << "Warning: Polygons are not supported in our MEDIT writer";
        }
        std::vector<NamedContainer> A[4] = {attr.points, {}, attr.facets, {}};
        write_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, A);
    }

    void write_medit(const std::string filename, const Volume& m, VolumeAttributes attr) {
        std::vector<vec3> verts(m.nverts());;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        FOR(v, m.nverts()) verts[v] = m.points[v];
        if (m.cell_type == Volume::TETRAHEDRON) {
            tets.resize(4 * m.ncells());
            FOR(t, m.ncells()) FOR(tv, 4) tets[4 * t + tv] = m.vert(t, tv);
        } else if (m.cell_type == Volume::HEXAHEDRON) {
            hexes.resize(8 * m.ncells());
            FOR(h, m.ncells()) FOR(hv, 8) hexes[8 * h + hv] = m.vert(h, hv);
        } else if (m.cell_type == Volume::WEDGE) {
            wedges.resize(6 * m.ncells());
            FOR(h, m.ncells()) FOR(hv, 6) wedges[6 * h + hv] = m.vert(h, hv);
        } else if (m.cell_type == Volume::PYRAMID) {
            pyramids.resize(5 * m.ncells());
            FOR(h, m.ncells()) FOR(hv, 5) pyramids[5 * h + hv] = m.vert(h, hv);
        } else {
            std::cerr << "Warning: Volume type : " << m.cell_type << "; not supported in our MEDIT writer";
        }
        std::vector<NamedContainer> A[4] = {attr.points, {}, attr.cell_facets, attr.cells};
        write_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, A);
    }


    PointSetAttributes read_medit(const std::string filename, PointSet& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        std::vector<NamedContainer> attrib[4];
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, attrib);
        m = PointSet();
        m.create_points(verts.size());
        FOR(v, verts.size()) m[v] = verts[v];
        return {attrib[0]};
    }

    PolyLineAttributes read_medit(const std::string filename, PolyLine& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        std::vector<NamedContainer> attrib[4];
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, attrib);
        m = PolyLine();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_edges(edges.size()/2);
        FOR(e, m.nedges()) FOR(ev, 2) m.vert(e, ev) = edges[2 * e + ev];
        return {attrib[0], attrib[1]};
    }

    SurfaceAttributes read_medit(const std::string filename, Triangles& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        std::vector<NamedContainer> attrib[4];
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, attrib);
        m = Triangles();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_facets(tris.size() / 3);
        FOR(t, m.nfacets()) FOR(tv, 3) m.vert(t, tv) = tris[3 * t + tv];
        return {attrib[0], attrib[2], {}};
    }

    SurfaceAttributes read_medit(const std::string filename, Quads& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        std::vector<NamedContainer> attrib[4];
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, attrib);
        m = Quads();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_facets(quads.size() / 4);
        FOR(q, m.nfacets()) FOR(qv, 4) m.vert(q, qv) = quads[4 * q + qv];
        return {attrib[0], attrib[2], {}};
    }

    SurfaceAttributes read_medit(const std::string filename, Polygons& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        std::vector<NamedContainer> attrib[4];
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, attrib);
        m = Polygons();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_facets(tris.size() / 3, 3);
        FOR(t, m.nfacets()) FOR(tv, 3) m.vert(t, tv) = tris[3 * t + tv];

        int off = m.create_facets(quads.size() / 4, 4);
        FOR(q, quads.size() / 4) FOR(qv, 4) m.vert(off + q, qv) = quads[4 * q + qv];
        return {attrib[0], attrib[2], {}};
    }

    VolumeAttributes read_medit(const std::string filename, Tetrahedra& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        std::vector<NamedContainer> attrib[4];
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, attrib);
        m = Tetrahedra();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];
        m.create_cells(tets.size() / 4);
        FOR(t, m.ncells()) FOR(tv, 4) m.vert(t, tv) = tets[4 * t + tv];
        return {attrib[0], attrib[3], attrib[2], {}};
    }

    VolumeAttributes read_medit(const std::string filename, Hexahedra& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        std::vector<NamedContainer> attrib[4];
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, attrib);
        m = Hexahedra();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_cells(hexes.size() / 8);
        FOR(h, m.ncells()) FOR(hv, 8) m.vert(h, hv) = hexes[8 * h + hv];
        return {attrib[0], attrib[3], attrib[2], {}};
    }

    VolumeAttributes read_medit(const std::string filename, Wedges& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        std::vector<NamedContainer> attrib[4];
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, attrib);
        m = Wedges();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_cells(wedges.size() / 6);
        FOR(h, m.ncells()) FOR(hv, 6) m.vert(h, hv) = wedges[5 * h + hv];
        return {attrib[0], attrib[3], attrib[2], {}};
    }

    VolumeAttributes read_medit(const std::string filename, Pyramids& m) {
        std::vector<vec3> verts;
        std::vector<int> edges, tris, quads, tets, hexes, wedges, pyramids;
        std::vector<NamedContainer> attrib[4];
        read_medit_format(filename, verts, edges, tris, quads, tets, hexes, wedges, pyramids, attrib);
        m = Pyramids();
        m.points.create_points(verts.size());
        FOR(v, verts.size()) m.points[v] = verts[v];

        m.create_cells(pyramids.size() / 5);
        FOR(h, m.ncells()) FOR(hv, 5) m.vert(h, hv) = pyramids[5 * h + hv];

        return {attrib[0], attrib[3], attrib[2], {}};
    }

}

