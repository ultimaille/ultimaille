#include <iostream>
#include <ultimaille/all.h>

using namespace UM;

int main(int argc, char** argv) {
    Tetrahedra m;
    std::cerr << "Reading... ";
    read_by_extension(argv[1], m);
    std::cerr << "Ok; connecting... ";
    m.connect();


    double ave = 0;
    std::cerr << "done" << std::endl;
    for(auto h:m.iter_halfedges()) ave += (h.to().pos()-h.from().pos()).norm(); 
    ave/=3*m.nfacets();

    std::cerr << "ave: " << ave << std::endl;

    CellCornerAttribute<double> dist(m,1e20);
    auto cmp = [&](int left, int right) { return dist[left] > dist[right]; };
    std::priority_queue<int, std::vector<int>, decltype(cmp)> queue(cmp);
    queue.push(0);
    dist[0]=0;
    while(!queue.empty()){
        Volume::Corner cur (m,queue.top());
        queue.pop();
        Volume::Halfedge h(m,-1);

        for (auto f: cur.cell().iter_facets()) {
            for (auto cir:f.iter_halfedges()) {
                if (cir.from_corner() == cur) h = cir;
            }
        }
        um_assert(h.active());

        for(int lh=0;lh<3;lh++){
            auto opp = h.opposite_c();

            if (dist[opp.from_corner()]>dist[cur]+(cur.vertex().pos()-opp.from_corner().vertex().pos()).norm()){
                dist[opp.from_corner()]=dist[cur]+(cur.vertex().pos()-opp.from_corner().vertex().pos()).norm();
                queue.push(opp.from_corner());
            }

            if (opp.active() && dist[opp.to_corner()]>dist[cur]){
                dist[opp.to_corner()]=dist[cur];
                queue.push(opp.to_corner());
            }
            h = h.prev().opposite_f();
        }
        if (queue.empty()){
            int far_corner=0;
            for(auto c:m.iter_corners()) if (dist[far_corner]<dist[c]) far_corner  = c;
            if (dist[far_corner]>10*ave){
                dist[far_corner] = 0 ;
                queue.push(far_corner);
            }
        }
    }
    write_by_extension("out.geogram",m,{{},{},{},{{"dist",dist.ptr}}});
}

