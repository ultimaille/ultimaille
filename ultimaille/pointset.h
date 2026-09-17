#ifndef __POINTSET_H__
#define __POINTSET_H__
#include <vector>
#include <memory>
#include "attribute_base.h"
#include "algebra/vec.h"

namespace UM {


struct PolyLine;
struct Surface;
struct Volume;

struct PointSetAttributes;
struct PolyLineAttributes;
struct SurfaceAttributes; 
struct VolumeAttributes;

    struct PointSet {
        PointSet() : data(new std::vector<vec3>()), attr(new std::vector<std::weak_ptr<ContainerBase>>()) {}
//      PointSet(std::shared_ptr<std::vector<vec3> > ext) : data(ext) {}

        PointSet(const PointSet &p)            = default; // We need to be able to share point sets, therefore we allow copying of the pointers
        PointSet(PointSet &&p)                 = default; // N.B. attrs pointer is also shared
        PointSet& operator=(const PointSet& p) = default;

        int size() const { return data->size(); }
        vec3& operator[](const int i) { return data->at(i); }
        const vec3& operator[](const int i) const { return data->at(i); }
        int use_count() { return data.use_count(); }

        void resize(const int n);
        int push_back(const vec3 &p);

        template <typename T> void delete_points(const T &to_kill);
        template <typename T> void delete_points(const T &to_kill, std::vector<int> &old2new);
        int create_points(const int n);

        using       iterator = std::vector<vec3>::iterator;
        using const_iterator = std::vector<vec3>::const_iterator;

        iterator begin() { return data->begin(); }
        iterator end()   { return data->end();   }
        const_iterator begin() const { return data->begin(); }
        const_iterator end()   const { return data->end();   }


        template <typename T> struct Attribute : GenericAttribute<T> {
            Attribute(T def = T());
            Attribute(PointSet &pts, T def = T());
            Attribute(PolyLine &m,   T def = T());
            Attribute(Surface &m,    T def = T());
            Attribute(Volume  &m,    T def = T());
            Attribute(const PointSet &pts, T def = T());
            Attribute(const PolyLine &m,   T def = T());
            Attribute(const Surface &m,    T def = T());
            Attribute(const Volume  &m,    T def = T());

            Attribute(std::string name, PointSetAttributes &attributes, PointSet &ps,  T def = T());
            Attribute(std::string name, PolyLineAttributes &attributes, PolyLine &seg, T def = T());
            Attribute(std::string name, SurfaceAttributes  &attributes, Surface  &m,   T def = T());
            Attribute(std::string name, VolumeAttributes   &attributes, Volume   &m,   T def = T());

            bool bind(std::string name, PointSetAttributes &attributes, PointSet &ps  );
            bool bind(std::string name, PolyLineAttributes &attributes, PolyLine &seg );
            bool bind(std::string name, SurfaceAttributes  &attributes, Surface  &m   );
            bool bind(std::string name, VolumeAttributes   &attributes, Volume   &m   );

            virtual AttributeBase::TYPE kind() const { return AttributeBase::POINTS; }
        };

        void resize_attrs();
        void compress_attrs(const std::vector<int> &old2new);

        std::shared_ptr<std::vector<vec3>> data = nullptr;
        std::shared_ptr<std::vector<std::weak_ptr<ContainerBase>>> attr = nullptr;
    };

    template <typename T> using PointAttribute = PointSet::Attribute<T>;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    template <typename T> void PointSet::delete_points(const T &to_kill, std::vector<int> &old2new) {
        constexpr bool invocable = std::is_invocable_r_v<bool, T, int>;
        if constexpr (!invocable)
            assert(to_kill.size()==(size_t)size());
        old2new = std::vector<int>(size(),  -1);

        int new_nb_pts = 0;
        for (int v=0; v<size(); v++) {
            if constexpr (invocable) {
                if (to_kill(v)) continue;
            } else {
                if (to_kill[v]) continue;
            }
            data->at(new_nb_pts) = data->at(v);
            old2new[v] = new_nb_pts++;
        }
        data->resize(new_nb_pts);
        compress_attrs(old2new);
    }

    template <typename T> void PointSet::delete_points(const T &to_kill) {
        std::vector<int> old2new;
        delete_points(to_kill, old2new);
    }
}

#endif //__POINTSET_H__

