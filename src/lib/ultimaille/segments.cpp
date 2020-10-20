#include <iostream>
#include <algorithm>
#include <cassert>
#include "segments.h"
#include "attributes.h"

int Segments::nverts() const {
    return points.size();
}

int Segments::nsegments() const {
    assert(segments.size()%2==0);
    return segments.size()/2;
}

int Segments::vert(const int s, const int lv) const {
    assert(s>=0 && s<nsegments() && lv>=0 && lv<2);
    return segments[s*2 + lv];
}

int &Segments::vert(const int s, const int lv) {
    assert(s>=0 && s<nsegments() && lv>=0 && lv<2);
    return segments[s*2 + lv];
}

int Segments::create_segments(const int n) {
    segments.resize(segments.size()+n*2);
    resize_attrs();
    return nsegments()-n;
}

void Segments::resize_attrs() {
    for (auto &wp : attr)  if (auto spt = wp.lock())
        spt->resize(nsegments());
}

