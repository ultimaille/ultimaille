#include <iostream>
#include <algorithm>
#include <cassert>
#include "attributes.h"
#include "polyline.h"
#include "surface.h"
#include "volume.h"
#include "attr_binding.h"

namespace UM {

    void Surface::connect() {
        if (!conn) conn = std::make_unique<Connectivity>(*this);
        else conn->init();
    }

    void Surface::disconnect() {
        conn.reset();
    }

    void Surface::compact(bool delete_isolated_vertices) {
        if (!conn) return;
        um_assert(conn->active.ptr!=nullptr);
        std::vector<bool> to_kill = conn->active.ptr->data;
        to_kill.flip();
        disconnect(); // happy assert(conn==nullptr)!
        delete_facets(to_kill);
        if (delete_isolated_vertices)
            Surface::delete_isolated_vertices();
        connect();
    }

    bool Surface::Vertex::on_boundary() {
        assert(m.connected());
        for (Halfedge &he : iter_halfedges())
            if (he.on_boundary()) return true;

        return false;
    }

    Surface::Connectivity::Connectivity(Surface &m) : m(m), v2c(m, -1), c2f(m, -1), c2c(m, -1), active(m, true) {
        init();
    }

    void Surface::Connectivity::init() {
        active.fill(true);
        v2c.fill(-1);

        for (int f = 0; f < m.nfacets(); f++)
            for (int fc = 0; fc < m.facet_size(f); fc++) {
                int c = m.corner(f, fc);
                int v = m.vert(f, fc);
                c2f[c] = f;
                c2c[c] = v2c[v];
                v2c[v] = c;
            }
    }

    Surface::Facet Surface::Connectivity::create_facet(int *verts, int size) {
        int off = -1;
        auto tmp_conn = std::move(m.conn); // happy assert(conn==nullptr)!
        if (auto *mesh = dynamic_cast<Triangles*>(&m)) {
            assert(size==3);
            off = mesh->create_facets(1);
        } else if (auto *mesh = dynamic_cast<Quads*>(&m)){
            assert(size==4);
            off = mesh->create_facets(1);
        } else if (auto *mesh = dynamic_cast<Polygons*>(&m)){
            off = mesh->create_facets(1, size);
        } else
            um_assert(false);

        m.conn = std::move(tmp_conn);
        Facet f(m, off);
        assert(f.active());
        for (int lv = 0; lv < size; lv++) {
            m.vert(f, lv) = verts[lv];
            int h = m.corner(f, lv);
            c2f[h] = f;
            c2c[h] = v2c[verts[lv]];
            v2c[verts[lv]] = h;
        }
        return f;
    }

    Surface::Facet Surface::Connectivity::create_facet(std::initializer_list<int> verts) {
        std::vector<int> tmp = verts; // verts.begin() isn't necessarily a pointer to a contiguous memory chunk
        return create_facet(tmp.data(), tmp.size());
    }

    void Surface::Connectivity::change_from(Surface::Halfedge he, int new_vertex_id) {
        // edit the surface
        auto old_v = he.from();
        auto f = he.facet();
        int lh = he.id_in_facet();
        m.vert(f, lh) = new_vertex_id;

        // update the connectivity:
        // first detach he from old_v
        if (v2c[old_v] == he)
            v2c[old_v] = c2c[he];
        else {
            int it = v2c[old_v];
            while (c2c[it] != he) it = c2c[it];
            c2c[it] = c2c[he];
        }
        // then attach he to new_v
        c2c[he] = v2c[new_vertex_id];
        v2c[new_vertex_id] = he;
    }

    void Surface::delete_isolated_vertices()  {
        std::vector<bool> to_kill(nverts(), true);
        for (int f=0; f<nfacets(); f++)
            for (int lv=0; lv<facet_size(f); lv++)
                to_kill[vert(f, lv)] = false;
        delete_vertices(to_kill);
    }

    void Surface::resize_attrs() {
        for (auto &wp : attr_facets)  if (auto spt = wp.lock())
            spt->resize(nfacets());
        for (auto &wp : attr_corners) if (auto spt = wp.lock())
            spt->resize(ncorners());
    }

    void Surface::compress_attrs(const std::vector<bool> &facets_to_kill) {
        assert(facets_to_kill.size()==(size_t)nfacets());
        std::vector<int>  facets_old2new(nfacets(),  -1);
        std::vector<int> corners_old2new(ncorners(), -1);

        int new_nb_facets  = 0;
        int new_nb_corners = 0;

        for (int f=0; f<nfacets(); f++) {
            if (facets_to_kill[f]) continue;
            for (int lv=0; lv<facet_size(f); lv++)
                corners_old2new[corner(f, lv)] = new_nb_corners++;
            facets_old2new[f] = new_nb_facets++;
        }

        std::erase_if(attr_facets,  [](std::weak_ptr<GenericAttributeContainer> ptr) { return ptr.lock()==nullptr; }); // remove dead attributes
        std::erase_if(attr_corners, [](std::weak_ptr<GenericAttributeContainer> ptr) { return ptr.lock()==nullptr; });
        for (auto &wp : attr_facets) { // compress attributes
            auto spt = wp.lock();
            assert(spt!=nullptr);
            spt->compress(facets_old2new);
        }
        for (auto &wp : attr_corners) {
            auto spt = wp.lock();
            assert(spt!=nullptr);
            spt->compress(corners_old2new);
        }
    }

    void Surface::delete_vertices(const std::vector<bool> &to_kill) {
        assert(to_kill.size()==(size_t)nverts());
        std::vector<int> old2new;
        points.delete_points(to_kill, old2new); // conn.v2c is a PointAttribute, it is automatically updated here
        for (int &v : facets) {
            assert(old2new[v]>=0);
            v = old2new[v];
        }
    }

    void Surface::delete_facets(const std::vector<bool> &to_kill) {
        assert(!connected());
        assert(to_kill.size()==(size_t)nfacets());
        compress_attrs(to_kill);  // TODO: if to_kill comes from an attribute, compressing the attribute compromises the code below
        int new_nb_corners = 0;
        for (int f=0; f<nfacets(); f++) {
            if (to_kill[f]) continue;
            for (int lv=0; lv<facet_size(f); lv++)
                facets[new_nb_corners++] = vert(f, lv);
        }
        facets.resize(new_nb_corners);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    int Triangles::create_facets(const int n) {
        assert(!connected());
        facets.resize(facets.size()+n*3);
        resize_attrs();
        return nfacets()-n;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    int Quads::create_facets(const int n) {
        assert(!connected());
        facets.resize(facets.size()+n*4);
        resize_attrs();
        return nfacets()-n;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    int Polygons::create_facets(const int n, const int size) {
        assert(!connected());
        for (int i=0; i<n*size; i++)
            facets.push_back(0);
        for (int i=0; i<n; i++)
            offset.push_back(offset.back()+size);
        resize_attrs();
        return nfacets()-n;
    }

    void Polygons::delete_facets(const std::vector<bool> &to_kill) {
        assert(!connected());
        Surface::delete_facets(to_kill); // TODO: if to_kill comes from an attribute, Surface::delete_facets compacts it, thus compromising the code below
        int new_nb_facets = 0;
        for (int f=0; f<nfacets(); f++) {
            if (to_kill[f]) continue;
            offset[new_nb_facets+1] = offset[new_nb_facets] + facet_size(f);
            ++new_nb_facets;
        }
        offset.resize(new_nb_facets+1);
    }
}

