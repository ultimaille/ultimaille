#include <iostream>
#include <algorithm>
#include <queue>
#include <cstdlib>

#include <ultimaille/disjointset.h>
#include <ultimaille/mesh_io.h>
#include <ultimaille/surface.h>
#include <ultimaille/attributes.h>
#include <ultimaille/range.h>

#include <OpenNL_psm/OpenNL_psm.h>

using namespace UM;

double average_edge_length(const Surface &m) {
    double sum = 0;
    int nb = 0;
    for (int f : range(m.nfacets())) {
        for (int lv : range(m.facet_size(f))) {
            vec3 a = m.points[m.vert(f, lv)];
            vec3 b = m.points[m.vert(f, (lv+1)%m.facet_size(f))];
            sum += (a-b).norm();
            nb++;
        }
    }
    return sum/nb;
}


    inline int next_mod(int i, int imax) {
        return (i + 1) % imax;
    }

    inline int prev_mod(int i, int imax) {
        return (i + int(int(imax) - 1)) % imax;
    }

int main(int argc, char** argv) {
    if (2>argc) {
        std::cerr << "Usage: " << argv[0] << " model.obj" << std::endl;
        return 1;
    }
    Polygons m;
    SurfaceAttributes attributes = read_geogram(argv[1], m);
    CornerAttribute<vec2> uv("uti", attributes, m);

    for (int c : range(m.ncorners()))
        for (int d : range(2)) {
            uv[c][d] *= 2;
            if (std::abs(uv[c][d] - round(uv[c][d])) < .05) uv[c][d] = round(uv[c][d]);
        }


std::cerr << "ncorners: " << m.ncorners() << std::endl;

    int nb_init_facets = m.nfacets();
    PointAttribute<int> added_vertices(m.points);
    FacetAttribute<int> orig_tri_fid(m);
    PointAttribute<int> resp_facet(m.points);

    typedef std::pair<int /*v_id*/, vec2 /*uv*/> NewCorner;
    std::vector<std::vector<NewCorner> > new_corners(m.ncorners());



    SurfaceConnectivity fec(m);
    for (int c : range(m.ncorners())) {
        int opp = fec.opposite(c);

//      if (singular[fec.facet(h)]) continue; // the face is responsible for the edge if opp<h or if its opposite is singular
//        if (-1 != opp && !singular[fec.facet(opp)] && opp > h) continue;
          if (opp>=0 && /*!singular[fec.facet(opp)] &&*/ opp > c) continue;

        vec2 lU[2] = { uv[c], uv[fec.next(c)] };

        std::vector<double> coeff;
        for(int coord : range(2)) { // find barycentric coordinates of both u and v integer isos
            double v[2] = { lU[0][coord], lU[1][coord] };
            double from = floor(std::min(v[0], v[1])) + 1;
            double to   = std::max(v[0], v[1]);
            if (to-from > 1000) continue;
            for (double iso = from; iso < to; iso += 1.) {
                double c = (iso - v[0]) / (v[1] - v[0]); // v[0] is far from v[1] (U was pre-snapped to integers with .05 tolerance)
                if (/*!Numeric::is_nan(c) && */c > 0 && c < 1) {
                    // vec2 u = lU[0] + c*(lU[1] - lU[0]);
                    // u[coord] = iso;
                    coeff.push_back(c);
                }
            }
        }
        std::sort(coeff.begin(), coeff.end(), std::less<double>()); // we need to sort in order to merge close values


                std::vector<vec3> pts;
                std::vector<vec2> lu;
                std::vector<vec2> luopp;

                for(int i: range(coeff.size())) {
                    // a + c*(b-a) returns exactly a when a==b and no NaNs involved
                    // no need to worry about the cases when c==0 and c==1, because of the presnapping of the parameterization
                    vec3 pt = m.points[fec.from(c)] + coeff[i] * (m.points[fec.to(c)] - m.points[fec.from(c)]);

                    vec2 u  = lU[0] + coeff[i] * (lU[1] - lU[0]);
                    vec2 uopp;
                    if (opp>=0) uopp = uv[fec.next(opp)] + coeff[i] * (uv[opp] - uv[fec.next(opp)]);

                    for (int d : range(2)) { // we must guarantee that new vertices have (at least) one integer component
                        if (std::abs(u[d] - round(u[d])) < 1e-10) {
                            u[d] = round(u[d]);
                        }

                        if (opp>=0 && std::abs(uopp[d] - round(uopp[d])) < 1e-10) {
                            uopp[d] = round(uopp[d]);
                        }
                    }

                    pts.push_back(pt);
                    lu.push_back(u);
                    if (opp>=0) luopp.push_back(uopp);
                }

                // create vertices
                int off = m.nverts();
                m.points.resize(off+pts.size());

                //m->vertices.create_vertices(pts.size());
                for(int i : range(pts.size())) {
                    added_vertices[i+off] = true;
                    resp_facet[i+off] = orig_tri_fid[fec.c2f[c]];
                }
                for (int i : range(pts.size())) {
                    m.points[off + i] = pts[i];
                    new_corners[c].push_back(NewCorner(off + i, lu[i]));
                    if (opp>=0)
                        new_corners[opp].push_back(NewCorner(off + i, luopp[i]));
                }
                if (opp>=0)
                    std::reverse(new_corners[opp].begin(), new_corners[opp].end());


    }


    for( int f : range( nb_init_facets)) {
        std::vector<int> polyV;
        std::vector<vec2> poly_uv;
        for( int fc : range(m.facet_size(f))) {
            polyV.push_back(m.vert(f, fc));
            int c = m.facet_corner(f, fc);
            poly_uv.push_back(uv[c]);
            for( int i : range( new_corners[c].size())) {
                polyV.push_back(new_corners[c][i].first);
                poly_uv.push_back(new_corners[c][i].second);
            }
        }
        int nf = m.create_facets(1, polyV.size());
        for (int i:range(polyV.size()))
            m.vert(nf, i) = polyV[i];
        //            m->facets.attributes().copy_item(nf ,f);
        for( int fc : range(m.facet_size(nf)))
            uv[m.facet_corner(nf, fc)] = poly_uv[fc];
    }

{
std::cerr << "ncorners: " << m.ncorners() << std::endl;
        std::vector<bool> to_kill(nb_init_facets, true); // kill old (pre-split) facets
        to_kill.resize(m.nfacets(), false);
      m.delete_facets(to_kill);
}


        std::vector<bool> to_kill(m.nfacets(), false);
        for (int f : range(m.nfacets())) {
            //            if (singular[f]) continue;
            int nbc = m.facet_size(f);


            // STEP 1: find a couple (coordinate, iso value) to split the facet
            int coord = int(-1);
            double iso = 0;
            for (int test_coord : range(2)) {
                double min_v =  1e20;
                double max_v = -1e20;
                for(int fc : range( nbc)) {
                    min_v = std::min(min_v, uv[m.facet_corner(f, fc)][test_coord]);
                    max_v = std::max(max_v, uv[m.facet_corner(f, fc)][test_coord]);
                }
                if (floor(min_v)+1 < max_v) {  // floor(min_v)+1 is the first integer value strictly superior to min_v
                    coord = test_coord;
                    iso = floor(min_v)+1;
                    break;
                }
            }
            if (coord == int(-1)) continue;

            // STEP 2: if STEP 1 succedeed, compute the extremities of the new edge
            int cut[2] = { -1, -1 };

            for(int fc : range(nbc)) {
                if (uv[m.facet_corner(f, fc)][coord] != iso) continue;
                if (cut[0] == -1) {
                    cut[0] = fc;
                } else if (fc != next_mod(cut[0], nbc) && next_mod(fc, nbc) != cut[0]) { // no biangles
                    cut[1] = fc;
                }
            }
            if (cut[1] == -1) continue;

            { // let us check that the cut separates the facet
                bool inf1=true, sup1=true, inf2=true, sup2=true;
                for (int fc = cut[0]+1; fc<cut[1]; fc++) {
                    double tex = uv[m.facet_corner(f, fc)][coord];
                    inf1 = inf1 && (tex<iso);
                    sup1 = sup1 && (tex>iso);
                }
                for (int fc = cut[1]+1; fc<cut[0]+nbc; fc++) {
                    double tex = uv[m.facet_corner(f, fc%nbc)][coord];
                    inf2 = inf2 && (tex<iso);
                    sup2 = sup2 && (tex>iso);
                }
                if (!( (inf1 && sup2) || (sup1 && inf2) )) {
                    //                    plop("WARNING: facet cannot be properly cut by the iso");
                    continue;
                }
            }


            to_kill[f] = true;



            // STEP 3: compute new vertices (iso-integer value of uv) on the new edge
            std::vector<vec3> nv_pts;
            std::vector<vec2> nv_uv;
            {
                vec3 lX[2];
                for(int i : range( 2)) lX[i] = m.points[m.vert(f, cut[i])];
                vec2 lU[2];
                for(int i: range(2)) lU[i] = uv[m.facet_corner(f, cut[i])];
                std::vector<double> coeff;
                double v[2] = { lU[0][(coord+1)%2], lU[1][(coord+1)%2] }; // recall that coord is the cutting dimension

                for (double cur_iso = ceil(std::min(v[0], v[1])); cur_iso < std::max(v[0], v[1]); cur_iso += 1.0) {
                    double c = (cur_iso - v[0]) / (v[1] - v[0]);
                    if (/*!Numeric::is_nan(c) && */c > 0 && c < 1)
                        coeff.push_back(c);
                }

                std::sort(coeff.begin(), coeff.end(), std::less<double>());
                for (int i : range( coeff.size())) {
                    vec3 x = lX[0] + coeff[i] * (lX[1] - lX[0]); // it guarantees x==lX[0] when lX[0]==lX[1]
                    vec2 u = lU[0] + coeff[i] * (lU[1] - lU[0]); // no need to worry about coeff[i]==0 and coeff[i]==1 because of the parameterization pre-snapping
                    nv_pts.push_back(x);
                    nv_uv.push_back(u);
                }
                // new vertices must have only integer values of uv --- remove possible numerical imprecision
                for (int i : range( nv_pts.size())) {
                    for (int d : range(2)) {
                        nv_uv[i][d] = round(nv_uv[i][d]);
                    }
                }
            }




            // STEP 4: create new vertices and new faces
            int off = m.nverts();
            m.points.resize(off+nv_pts.size());
            for (int i : range( nv_pts.size())) {
                resp_facet[off+i] = orig_tri_fid[f];
                added_vertices[off+i] = true;
                m.points[off + i] = nv_pts[i];
            }
            for(int half : range( 2)) {
                std::vector <int> lv;
                std::vector <vec2> luv;

                // add original vertices
                int cir = cut[half];
                do {
                    lv.push_back(m.vert(f, cir));
                    luv.push_back(uv[m.facet_corner(f, cir)]);
                    cir = next_mod(cir, nbc);
                } while (cir != cut[(half + 1) % 2]);
                lv.push_back(m.vert(f, cir));
                luv.push_back(uv[m.facet_corner(f, cir)]);

                // add new vertices
                for( int i : range( nv_pts.size())) {
                    int ind = i;
                    if (half == 0) ind = nv_pts.size() - 1 - i;
                    lv.push_back(off + ind);
                    luv.push_back(nv_uv[ind]);
                }

//                int fid = m->facets.create_polygon(lv);
                int fid = m.create_facets(1, lv.size());
                for (int i : range(lv.size())) m.vert(fid, i) = lv[i];
//                m->facets.attributes().copy_item(fid, f);
                for (int fc : range(m.facet_size(fid))) {
                    uv[m.facet_corner(fid, fc)] = luv[fc];
                }
                to_kill.push_back(false);
            }
        }
        m.delete_facets(to_kill);








{ // no guarantees on the quality, simple fan-triangulate

        std::vector<bool> to_kill(m.nfacets(), true);
        int nb_facets = m.nfacets();
    int nb_triangles = 0;
        for (int f=0; f<nb_facets; f++)
        nb_triangles += (m.facet_size(f) - 2);

        to_kill.resize(nb_triangles+nb_facets, false);
  //      std::cerr << "deb1: " << m.nfacets() << " " << nb_triangles << std::endl;

    int off = m.create_facets(nb_triangles, 3);

    int cnt = 0;
    for (int f=0; f<nb_facets; f++) {
    int fs = m.facet_size(f);
        int bestv = 0;
        double bestval = 0;

        for (int lv : range(fs)) {
            vec3 e1 = (m.points[m.vert(f, lv)] - m.points[m.vert(f, (lv-1+fs)%fs)]).normalize();
            vec3 e2 = (-m.points[m.vert(f, lv)] + m.points[m.vert(f, (lv+1+fs)%fs)]).normalize();
            if (std::abs(e1*e2)>bestval) {bestval = std::abs(e1*e2); bestv = lv;}
        }
        std::cerr << bestv << std::endl;








        int v0 = m.vert(f, bestv);
        for (int lv=1; lv+1<fs; lv++) {
//          assert(cnt<=nb_triangles);
//        std::cerr << "deb2: " << f << " " << cnt << " " << m.nfacets() << std::endl;
              m.vert(off+cnt, 0) = v0;
            m.vert(off+cnt, 1) = m.vert(f, (bestv+lv)%fs);
            m.vert(off+cnt, 2) = m.vert(f, (bestv+lv+1)%fs);
            uv[m.facet_corner(off+cnt, 0)] = uv[m.facet_corner(f, bestv)];
            uv[m.facet_corner(off+cnt, 1)] = uv[m.facet_corner(f, (bestv+lv)%fs)];
            uv[m.facet_corner(off+cnt, 2)] = uv[m.facet_corner(f, (bestv+lv+1)%fs)];
            cnt++;
        }
    }
    m.delete_facets(to_kill);

}

{
    SurfaceConnectivity fec(m);
    FacetAttribute<int> chart(m);
    std::get<1>(attributes).push_back(std::make_pair("chart", chart.ptr));
    DisjointSet ds(m.nfacets());
    for (int c: range(m.ncorners())) {
        int opp  = fec.opposite(c);
        if (opp<0) continue;
        if (uv[c][0]==uv[fec.next(c)][0] || uv[c][1]==uv[fec.next(c)][1]) continue;
        ds.merge(fec.c2f[c], fec.c2f[opp]);
    }

    std::vector<int> id2set;
    ds.get_sets_id(id2set);

    for (int f : range(m.nfacets())) chart[f] = id2set[f];





    PointAttribute<int> verts_to_remove(m.points);
    for (int v: range(m.nverts())) verts_to_remove[v] = 1;
    for (int c: range(m.ncorners())) {
        if (uv[c][0]==int(uv[c][0]) && uv[c][1]==int(uv[c][1])) verts_to_remove[fec.from(c)] = 0;
    }
    std::get<0>(attributes).push_back(std::make_pair("verts_to_remove", verts_to_remove.ptr));
}

    write_geogram("read_test.geogram", m, attributes);

    return 0;
}

