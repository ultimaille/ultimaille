#ifndef __ATTRIBUTES_H__
#define __ATTRIBUTES_H__
#include <vector>
#include <memory>
#include <cassert>
#include "syntactic-sugar/assert.h"
#include "attribute_base.h"

namespace UM {
    struct PolyLine;
    struct Surface;
    struct Volume;

    struct SurfaceAttributes;
    struct PointSetAttributes;
    struct PolyLineAttributes;
    struct VolumeAttributes;

    struct NamedContainer {
        std::string name;
        std::shared_ptr<ContainerBase> ptr;
    };

    struct NamedAttribute {
        std::string name;
        AttributeBase& attribute;
    };

    struct PointSetAttributes {
        PointSetAttributes() = default;
        PointSetAttributes(PointSetAttributes& p)        = default;
        PointSetAttributes(PointSetAttributes&& p)       = default;
        PointSetAttributes(const PointSetAttributes& p)  = default;
        PointSetAttributes& operator=(const PointSetAttributes& p)  = default;

        PointSetAttributes(std::vector<NamedContainer> list) : points(list) {}

        PointSetAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto na : list) {
                um_assert(na.attribute.kind()==AttributeBase::POINTS);
                points.emplace_back(na.name, na.attribute.get_ptr());
            }
        }

        std::vector<NamedContainer> points;
    };

    struct PolyLineAttributes {
        PolyLineAttributes() = default;
        PolyLineAttributes(PolyLineAttributes& p)        = default;
        PolyLineAttributes(PolyLineAttributes&& p)       = default;
        PolyLineAttributes(const PolyLineAttributes& p)  = default;
        PolyLineAttributes& operator=(const PolyLineAttributes& p)  = default;

        PolyLineAttributes(std::vector<NamedContainer> points, std::vector<NamedContainer> edges) : points(points), edges(edges) {}

        PolyLineAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto na : list) {
                switch (na.attribute.kind()) {
                    case AttributeBase::POINTS:  points.emplace_back(na.name, na.attribute.get_ptr()); break;
                    case AttributeBase::EDGES:    edges.emplace_back(na.name, na.attribute.get_ptr()); break;
                    default: um_assert(false);
                }
            }
        }

        std::vector<NamedContainer> points = {}, edges = {};
    };

    struct SurfaceAttributes {
        SurfaceAttributes() = default;
        SurfaceAttributes(SurfaceAttributes& p)        = default;
        SurfaceAttributes(SurfaceAttributes&& p)       = default;
        SurfaceAttributes(const SurfaceAttributes& p)  = default;
        SurfaceAttributes& operator=(const SurfaceAttributes& p)  = default;

        SurfaceAttributes(std::vector<NamedContainer> points, std::vector<NamedContainer> facets, std::vector<NamedContainer> corners) : points(points), facets(facets), corners(corners) {}

        SurfaceAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto na : list) {
                switch (na.attribute.kind()) {
                    case AttributeBase::POINTS:   points.emplace_back(na.name, na.attribute.get_ptr()); break;
                    case AttributeBase::FACETS:   facets.emplace_back(na.name, na.attribute.get_ptr()); break;
                    case AttributeBase::CORNERS: corners.emplace_back(na.name, na.attribute.get_ptr()); break;
                    default: um_assert(false);
                }
            }
        }

        std::vector<NamedContainer> points = {}, facets = {}, corners = {};
    };

    struct VolumeAttributes {
        VolumeAttributes() = default;
        VolumeAttributes(VolumeAttributes& p)        = default;
        VolumeAttributes(VolumeAttributes&& p)       = default;
        VolumeAttributes(const VolumeAttributes& p)  = default;
        VolumeAttributes& operator=(const VolumeAttributes& p)  = default;

        VolumeAttributes(std::vector<NamedContainer> points,
                         std::vector<NamedContainer> cells,
                         std::vector<NamedContainer> cell_facets,
                         std::vector<NamedContainer> cell_corners) : points(points),
                                                                     cells(cells),
                                                                     cell_facets(cell_facets),
                                                                     cell_corners(cell_corners) {
        }

        VolumeAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto na : list) {
                switch (na.attribute.kind()) {
                    case AttributeBase::POINTS:            points.emplace_back(na.name, na.attribute.get_ptr()); break;
                    case AttributeBase::CELLS:              cells.emplace_back(na.name, na.attribute.get_ptr()); break;
                    case AttributeBase::CELLFACETS:   cell_facets.emplace_back(na.name, na.attribute.get_ptr()); break;
                    case AttributeBase::CELLCORNERS: cell_corners.emplace_back(na.name, na.attribute.get_ptr()); break;
                    default: um_assert(false);
                }
            }
        }

        std::vector<NamedContainer> points = {}, cells = {}, cell_facets = {}, cell_corners = {};
    };
}

#endif //__ATTRIBUTES_H__

