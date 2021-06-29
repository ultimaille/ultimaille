#include <iostream>
#include <algorithm>
#include <cassert>
#include "polyline.h"
#include "attributes.h"

namespace UM {
    int PolyLine::nverts() const {
        return points.size();
    }

    int PolyLine::nsegments() const {
        assert(segments.size()%2==0);
        return segments.size()/2;
    }

    int PolyLine::vert(const int s, const int lv) const {
        assert(s>=0 && s<nsegments() && lv>=0 && lv<2);
        return segments[s*2 + lv];
    }

    int &PolyLine::vert(const int s, const int lv) {
        assert(s>=0 && s<nsegments() && lv>=0 && lv<2);
        return segments[s*2 + lv];
    }

    int PolyLine::create_segments(const int n) {
        segments.resize(segments.size()+n*2);
        resize_attrs();
        return nsegments()-n;
    }

    void PolyLine::resize_attrs() {
        for (auto &wp : attr)  if (auto spt = wp.lock())
            spt->resize(nsegments());
    }

    void PolyLine::compress_attrs(const std::vector<bool> &segments_to_kill) {
        assert(segments_to_kill.size()==(size_t)nsegments());
        std::vector<int>  segments_old2new(nsegments(),  -1);

        int new_nb_segments = 0;

        for (int s=0; s<nsegments(); s++)
            if (!segments_to_kill[s])
                segments_old2new[s] = new_nb_segments++;

        for (auto &wp : attr) if (auto spt = wp.lock())
            spt->compress(segments_old2new);
    }

    void PolyLine::delete_vertices(const std::vector<bool> &to_kill) {
        assert(to_kill.size()==(size_t)nverts());
        std::vector<bool> segments_to_kill(nsegments(), false);

        for (int s=0; s<nsegments(); s++)
            for (int lv=0; lv<2; lv++)
                if (to_kill[vert(s, lv)])
                    segments_to_kill[s] = true;

        delete_segments(segments_to_kill);

        std::vector<int> old2new;
        points.delete_points(to_kill, old2new);
        for (int &v : segments)
            v = old2new[v];
    }

    void PolyLine::delete_segments(const std::vector<bool> &to_kill) {
        assert(to_kill.size()==(size_t)nsegments());
        compress_attrs(to_kill);

        int new_nb_endpoints = 0;
        for (int s=0; s<nsegments(); s++)
            if (!to_kill[s])
                for (int lv=0; lv<2; lv++)
                    segments[new_nb_endpoints++] = vert(s, lv);
        segments.resize(new_nb_endpoints);
    }
}

