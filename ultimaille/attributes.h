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
        PointSetAttributes(PointSetAttributes& p)       = default;
        PointSetAttributes(PointSetAttributes&& p)      = default;
        PointSetAttributes(const PointSetAttributes& p) = default;
        PointSetAttributes& operator=(const PointSetAttributes& p) = default;

        PointSetAttributes(AttributeMap points) : points(std::move(points)) {}

        PointSetAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto& na : list)
                points[na.name] = na.attribute.get_ptr();
        }

        AttributeMap points;
    };

    struct PolyLineAttributes {
        PolyLineAttributes() = default;
        PolyLineAttributes(PolyLineAttributes& p)       = default;
        PolyLineAttributes(PolyLineAttributes&& p)      = default;
        PolyLineAttributes(const PolyLineAttributes& p) = default;
        PolyLineAttributes& operator=(const PolyLineAttributes& p) = default;

        PolyLineAttributes(AttributeMap points, AttributeMap edges) : points(std::move(points)), edges(std::move(edges)) {}

        PolyLineAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto& na : list) {
                switch (na.attribute.kind()) {
                    case AttributeBase::POINTS: points[na.name] = na.attribute.get_ptr(); break;
                    case AttributeBase::EDGES:   edges[na.name] = na.attribute.get_ptr(); break;
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
        SurfaceAttributes(SurfaceAttributes& p)       = default;
        SurfaceAttributes(SurfaceAttributes&& p)      = default;
        SurfaceAttributes(const SurfaceAttributes& p) = default;
        SurfaceAttributes& operator=(const SurfaceAttributes& p) = default;

        SurfaceAttributes(AttributeMap points, AttributeMap facets, AttributeMap corners) : points(std::move(points)), facets(std::move(facets)), corners(std::move(corners)) {}

        SurfaceAttributes(std::initializer_list<NamedAttribute> list) {
            for (auto& na : list) {
                switch (na.attribute.kind()) {
                    case AttributeBase::POINTS:   points[na.name] = na.attribute.get_ptr(); break;
                    case AttributeBase::FACETS:   facets[na.name] = na.attribute.get_ptr(); break;
                    case AttributeBase::CORNERS: corners[na.name] = na.attribute.get_ptr(); break;
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
        VolumeAttributes(VolumeAttributes& p)       = default;
        VolumeAttributes(VolumeAttributes&& p)      = default;
        VolumeAttributes(const VolumeAttributes& p) = default;
        VolumeAttributes& operator=(const VolumeAttributes& p) = default;

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
            for (const auto& na : list) {
                switch (na.attribute.kind()) {
                    case AttributeBase::POINTS:            points[na.name] = na.attribute.get_ptr(); break;
                    case AttributeBase::CELLS:              cells[na.name] = na.attribute.get_ptr(); break;
                    case AttributeBase::CELLFACETS:   cell_facets[na.name] = na.attribute.get_ptr(); break;
                    case AttributeBase::CELLCORNERS: cell_corners[na.name] = na.attribute.get_ptr(); break;
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

