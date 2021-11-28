#include <iostream>
#include <limits>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "ultimaille/io/obj.h"

namespace UM {
    // supposes .obj file to have "f " entries without slashes
    // TODO: improve the parser
    // TODO: export vn (corner and point) and vt (corner) attributes
    SurfaceAttributes read_wavefront_obj(const std::string filename, Polygons &m) {
        SurfaceAttributes sa;
        std::vector<vec3> VN;
        std::vector<vec2> VT;
        std::vector<std::vector<int>> VTID;
        std::vector<std::vector<int>> VNID;

        m = Polygons();
        std::ifstream in;
        in.open (filename, std::ifstream::in);
        if (in.fail())
            throw std::runtime_error("Failed to open " + filename);
        std::string line;
        while (!in.eof()) {
            std::getline(in, line);
            std::istringstream iss(line.c_str());
            std::string type_token;
            iss >> type_token;

            if (type_token=="v") {
                vec3 v;
                for (int i=0;i<3;i++) iss >> v[i];
                m.points.data->push_back(v);
            } else if (type_token=="vn") {
                vec3 v;
                for (int i=0;i<3;i++) iss >> v[i];
                VN.push_back(v);
            } else if (type_token=="vt") {
                vec2 v;
                for (int i=0;i<2;i++) iss >> v[i];
                VT.push_back(v);
            } else if (type_token=="f") {
                std::vector<int> vid;
                std::vector<int> vnid;
                std::vector<int> vtid;
                int tmp;
                while (1) { // in wavefront obj all indices start at 1, not zero
                    while (!iss.eof() && !std::isdigit(iss.peek())) iss.get(); // skip (esp. trailing) white space
                    if (iss.eof()) break;
                    iss >> tmp;
                    vid.push_back(tmp-1);
                    if (iss.peek() == '/') {
                        iss.get();
                        if (iss.peek() == '/') {
                            iss.get();
                            iss >> tmp;
                            vnid.push_back(tmp-1);
                        } else {
                            iss >> tmp;
                            vtid.push_back(tmp-1);
                            if (iss.peek() == '/') {
                                iss.get();
                                iss >> tmp;
                                vnid.push_back(tmp-1);
                            }
                        }
                    }
                }
                VTID.push_back(vtid);
                VNID.push_back(vnid);

                int off_f = m.create_facets(1, vid.size());
                for (int i=0; i<static_cast<int>(vid.size()); i++)
                    m.vert(off_f, i) = vid[i];
            }
        }

        bool vt_pt_attr = ((int)VTID.size()==m.nfacets() && (int)VT.size()==m.nverts()); // check whether tex_coord is a PointAttribute
        if (vt_pt_attr) for (int f=0; f<m.nfacets(); f++) {
            vt_pt_attr = vt_pt_attr && ((int)VTID[f].size()==m.facet_size(f));
            if (vt_pt_attr) for (int v=0; v<m.facet_size(f); v++) {
                vt_pt_attr = vt_pt_attr && (m.vert(f, v)==VTID[f][v]);
            }
        }

        if (vt_pt_attr) {
            PointAttribute<vec2> tex_coord(m.points);
            for (int v=0; v<m.nverts(); v++)
                tex_coord[v] = VT[v];
            sa.points.emplace_back("tex_coord", tex_coord.ptr);
        } else {
            bool vt_c_attr = ((int)VTID.size()==m.nfacets() && (int)VT.size()>0); // check whether tex_coord is a CornerAttribute
            if (vt_c_attr) for (int f=0; f<m.nfacets(); f++) {
                vt_c_attr = vt_c_attr && ((int)VTID[f].size()==m.facet_size(f));
            }
            if (vt_c_attr) {
                CornerAttribute<vec2> tex_coord(m);
                for (int f=0; f<m.nfacets(); f++) {
                    for (int v=0; v<m.facet_size(f); v++) {
                        um_assert(VTID[f][v]>=0 && VTID[f][v]<(int)VT.size());
                        tex_coord[m.corner(f, v)] = VT[VTID[f][v]];
                    }
                }
                sa.corners.emplace_back("tex_coord", tex_coord.ptr);
            }
        }

        in.close();
//      std::cerr << "#v: " << m.nverts() << " #f: "  << m.nfacets() << std::endl;
        return sa;
    }

    SurfaceAttributes read_wavefront_obj(const std::string filename, Triangles &m) {
        Polygons mpoly;
        SurfaceAttributes sa = read_wavefront_obj(filename, mpoly);

        std::vector<bool> to_kill(mpoly.nfacets(), false);
        for (int f=0; f<mpoly.nfacets(); f++)
            to_kill[f] = (3!=mpoly.facet_size(f));
        mpoly.delete_facets(to_kill);

        m.points = mpoly.points;
        m.facets = mpoly.facets;
        return sa;
    }

    SurfaceAttributes read_wavefront_obj(const std::string filename, Quads &m) {
        Polygons mpoly;
        SurfaceAttributes sa = read_wavefront_obj(filename, mpoly);

        std::vector<bool> to_kill(mpoly.nfacets(), false);
        for (int f=0; f<mpoly.nfacets(); f++)
            to_kill[f] = (4!=mpoly.facet_size(f));
        mpoly.delete_facets(to_kill);

        m.points = mpoly.points;
        m.facets = mpoly.facets;
        return sa;
    }

    void write_wavefront_obj(const std::string filename, const Surface &m, const SurfaceAttributes attr) {
        std::fstream out;
        out.open(filename, std::ios_base::out);
        if (out.fail())
            throw std::runtime_error("Failed to open " + filename);
//      out << std::fixed << std::setprecision(14);
        out << std::setprecision(std::numeric_limits<double>::max_digits10);
        for (int v=0; v<m.nverts(); v++)
            out << "v " << m.points[v] << std::endl;

        for (auto &pair : attr.points) { // export tex_coord per vertex
            if (pair.first!="tex_coord") continue;
            std::shared_ptr<GenericAttributeContainer> ptr = pair.second;
            if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<vec2>>(ptr); cont_ptr.get()!=nullptr) {
                std::vector<vec2> tmp = cont_ptr->data;
                um_assert((int)tmp.size()==m.nverts());
                for (int v=0; v<m.nverts(); v++)
                    out << "vt " << tmp[v] << std::endl;
                for (int f=0; f<m.nfacets(); f++) {
                    out << "f ";
                    for (int v=0; v<m.facet_size(f); v++)
                        out << (m.vert(f,v)+1) << "/" << (m.vert(f,v)+1) << " ";
                    out << std::endl;
                }
                goto clear;
            }
        }

        for (auto &pair : attr.corners) { // export tex_coord per corner
            if (pair.first!="tex_coord") continue;
            std::shared_ptr<GenericAttributeContainer> ptr = pair.second;
            if (auto cont_ptr = std::dynamic_pointer_cast<AttributeContainer<vec2>>(ptr); cont_ptr.get()!=nullptr) {
                std::vector<vec2> tmp = cont_ptr->data;
                um_assert((int)tmp.size()==m.ncorners());
                for (int v=0; v<m.ncorners(); v++)
                    out << "vt " << tmp[v] << std::endl;
                for (int f=0; f<m.nfacets(); f++) {
                    out << "f ";
                    for (int v=0; v<m.facet_size(f); v++)
                        out << (m.vert(f,v)+1) << "/" << (m.corner(f,v)+1) << " ";
                    out << std::endl;
                }
                goto clear;
            }
        }

        for (int f=0; f<m.nfacets(); f++) {
            out << "f ";
            for (int v=0; v<m.facet_size(f); v++)
                out << (m.vert(f,v)+1) << " ";
            out << std::endl;
        }
clear:
        out.close();
    }
}

