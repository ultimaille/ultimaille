#ifndef __ATTRIBUTES_H__
#define __ATTRIBUTES_H__
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

    struct PointSetAttributes {
        PointSetAttributes() = default;
        PointSetAttributes(PointSetAttributes& p)        = default;
        PointSetAttributes(PointSetAttributes&& p)       = default;
        PointSetAttributes(const PointSetAttributes& p)  = default;
        PointSetAttributes& operator=(const PointSetAttributes& p)  = default;

        PointSetAttributes(AttributeMap points) : points(std::move(points)) {}

        PointSetAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto& attribute : list)
                add_attribute(points, attribute, AttributeBase::POINTS);
        }

        AttributeMap points;
    };

    struct PolyLineAttributes {
        PolyLineAttributes() = default;
        PolyLineAttributes(PolyLineAttributes& p)        = default;
        PolyLineAttributes(PolyLineAttributes&& p)       = default;
        PolyLineAttributes(const PolyLineAttributes& p)  = default;
        PolyLineAttributes& operator=(const PolyLineAttributes& p)  = default;

        PolyLineAttributes(AttributeMap points, AttributeMap edges) : points(std::move(points)), edges(std::move(edges)) {}

        PolyLineAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto& attribute : list) {
                switch (attribute.attribute.kind()) {
                    case AttributeBase::POINTS: add_attribute(points, attribute, AttributeBase::POINTS); break;
                    case AttributeBase::EDGES:  add_attribute(edges, attribute, AttributeBase::EDGES);   break;
                    default:
                        um_assert(false);
                }
            }
        }
        AttributeMap points = {};
        AttributeMap edges  = {};
    };

    struct SurfaceAttributes {
        SurfaceAttributes() = default;
        SurfaceAttributes(SurfaceAttributes& p)        = default;
        SurfaceAttributes(SurfaceAttributes&& p)       = default;
        SurfaceAttributes(const SurfaceAttributes& p)  = default;
        SurfaceAttributes& operator=(const SurfaceAttributes& p)  = default;

        SurfaceAttributes(AttributeMap points, AttributeMap facets, AttributeMap corners) : points(std::move(points)), facets(std::move(facets)), corners(std::move(corners)) {}

        SurfaceAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto& attribute : list) {
                switch (attribute.attribute.kind()) {
                    case AttributeBase::POINTS:  add_attribute(points, attribute, AttributeBase::POINTS);   break;
                    case AttributeBase::FACETS:  add_attribute(facets, attribute, AttributeBase::FACETS);   break;
                    case AttributeBase::CORNERS: add_attribute(corners, attribute, AttributeBase::CORNERS); break;
                    default:
                        um_assert(false);
                }
            }
        }

        AttributeMap points  = {};
        AttributeMap facets  = {};
        AttributeMap corners = {};
    };

    struct VolumeAttributes {
        VolumeAttributes() = default;
        VolumeAttributes(VolumeAttributes& p)        = default;
        VolumeAttributes(VolumeAttributes&& p)       = default;
        VolumeAttributes(const VolumeAttributes& p)  = default;
        VolumeAttributes& operator=(const VolumeAttributes& p)  = default;

        VolumeAttributes(
                AttributeMap points,
                AttributeMap cells,
                AttributeMap cell_facets,
                AttributeMap cell_corners) :
            points(std::move(points)),
            cells(std::move(cells)),
            cell_facets(std::move(cell_facets)),
            cell_corners(std::move(cell_corners)) {}

        VolumeAttributes(std::initializer_list<NamedAttribute> list) {
            for (const auto& attribute : list) {
                switch (attribute.attribute.kind()) {
                    case AttributeBase::POINTS:      add_attribute(points, attribute, AttributeBase::POINTS);            break;
                    case AttributeBase::CELLS:       add_attribute(cells, attribute, AttributeBase::CELLS);              break;
                    case AttributeBase::CELLFACETS:  add_attribute(cell_facets, attribute, AttributeBase::CELLFACETS);   break;
                    case AttributeBase::CELLCORNERS: add_attribute(cell_corners, attribute, AttributeBase::CELLCORNERS); break;
                    default:
                        um_assert(false);
                }
            }
        }

        AttributeMap points       = {};
        AttributeMap cells        = {};
        AttributeMap cell_facets  = {};
        AttributeMap cell_corners = {};
    };
}

#endif //__ATTRIBUTES_H__

