#ifndef __ITERATOR_H__
#define __ITERATOR_H__
#include <iterator>

namespace UM::Iterator {

    template<typename T>
    struct VertexIterator {
        T& m;
        VertexIterator(T& _m) : m(_m) {}
        
        struct iterator {

            typename T::Vertex v;
            void operator++() { ++v.id; }
            bool operator!=(const iterator& rhs) const { return v != rhs.v; }
            typename T::Vertex& operator*() { return v; }
        };
    
        iterator begin() { return VertexIterator::iterator{ {m,0} }; }
        iterator end() { return iterator{ {m,m.nverts()} }; }
    };

    template<typename T>
    struct HalfedgeIterator {
        T& m;
        HalfedgeIterator(T& _m) : m(_m) {}
        
        struct iterator {

            typename T::Halfedge he;
            
            void operator++() {                 
                ++he.id;
                if (he.m.connected())
                    while (he.id < he.m.ncorners() && !he.active()) ++he.id; 
            }

            bool operator!=(const iterator& rhs) const { return he != rhs.he; }
            typename T::Halfedge& operator*() { return he; }
        };
    
        iterator begin() { 
            iterator res{ {m,0} };
            if (m.connected())
                if (m.ncorners() > 0 && !res.he.active()) ++res;
            return res;
        }

        iterator end() { return iterator{ {m,m.ncorners()} }; }
    };

    template<typename T>
    struct VertexHalfEdgeIterator {
        typename T::Vertex v;

        VertexHalfEdgeIterator(typename T::Vertex& _v) : v(_v) { /* TODO CHECK THAT assert(m.connected());*/ }
        
        struct iterator {

            typename T::Halfedge he;
            
            void operator++() {                 
                const auto &c2c = he.m.conn->c2c;
                do {
                    he = c2c[he];
                    if (he == -1) return;
                } while (!he.active());
            }

            bool operator!=(const iterator& rhs) const { return he != rhs.he; }
            typename T::Halfedge& operator*() { return he; }
        };
    
        iterator begin() { return { v.halfedge() }; }
        iterator end() { return { {v.m, -1} }; }
    };

    template<typename T>
    struct FacetHalfEdgeIterator {

        typename T::Facet f;

        FacetHalfEdgeIterator(typename T::Facet& _f) : f(_f) { /* TODO CHECK THAT assert(m.connected());*/ }
        
        struct iterator {

            typename T::Halfedge he;
            
            void operator++() {                 
                auto f = he.facet();
                if (++he == f.halfedge()+f.size()) he = -1;
            }

            bool operator!=(const iterator& rhs) const { return he != rhs.he; }
            typename T::Halfedge& operator*() { return he; }
        };
    
        iterator begin() { return { f.halfedge() }; }
        iterator end() { return { {f.m,-1} }; }
    };


}

#endif //__ITERATOR_H__