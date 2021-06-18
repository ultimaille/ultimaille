#include <catch2/catch.hpp>

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <ultimaille/all.h>

using namespace UM;
static const double ftol = 1e-13;

TEST_CASE("PointSet + PointSetAttributes IO test", "[geogram]") {
    static const std::string filename = "ultimaille-test-pointset.geogram";
    PointSet m;
    *m.data = {{0,.9,0}, {0,.8,0}, {1.,.3,0.}};
    PointAttribute<bool> vbool(m);
    for (int v : range(m.size()))
        vbool[v] = rand()&1;
    write_by_extension(filename, m, {{{"vbool", vbool.ptr}}});

    PointSet m2;
    PointSetAttributes attrs = read_by_extension(filename, m2);

    REQUIRE( m.size()==m2.size() );

    PointAttribute<bool> vbool2("vbool", attrs, m2);
    for (int v=0; v<m.size(); v++) {
        REQUIRE( (m[v]-m2[v]).norm()<ftol );
        REQUIRE( vbool[v]==vbool2[v] );
    }
}

TEST_CASE("PolyLine + PolyLineAttributes IO test", "[geogram]") {
    static const std::string filename = "ultimaille-test-polyline.geogram";
    PolyLine m;
    *m.points.data = {{0,.9,0}, {0,.8,0}, {1.,.3,0.}};
    m.segments = {0,1,1,2};
    PointAttribute<bool> vbool(m.points);
    SegmentAttribute<int> sint(m);
    for (int v : range(m.nverts()))
        vbool[v] = rand()&1;
    for (int cc : range(m.nsegments()))
        sint[cc] = rand();
    write_by_extension(filename, m, {{{"vbool", vbool.ptr}}, {{"sint", sint.ptr}}});

    PolyLine m2;
    PolyLineAttributes attrs = read_by_extension(filename, m2);

    REQUIRE( m.nverts()==m2.nverts() );
    REQUIRE( m.nsegments()==m2.nsegments() );

    PointAttribute<bool> vbool2("vbool", attrs, m2);
    for (int v=0; v<m.nverts(); v++) {
        REQUIRE( (m.points[v]-m2.points[v]).norm()<ftol );
        REQUIRE( vbool[v]==vbool2[v] );
    }

    SegmentAttribute<int> sint2("sint", attrs, m2);
    for (int c=0; c<m.nsegments(); c++) {
        for (int lv=0; lv<2; lv++) {
            REQUIRE( m.vert(c, lv)==m2.vert(c, lv) );
        }
        REQUIRE( std::fabs(sint[c]-sint2[c])<ftol );
    }
}

TEST_CASE("Triangles + SurfaceAttributes IO test", "[geogram]") {
    static const std::string filename = "ultimaille-test-triangles.geogram";
    Triangles m;
    *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}};
    m.facets = {0,1,2,2,1,3};

    PointAttribute<bool> vbool(m.points);
    FacetAttribute<vec3> fvec3(m);
    CornerAttribute<int> cint(m);

    for (int v : range(m.nverts()))
        vbool[v] = rand()&1;

    for (int f : range(m.nfacets()))
        fvec3[f] = m.util.bary_verts(f);

    for (int c : range(m.ncorners()))
        cint[c] = rand();

    write_by_extension(filename, m, {{{"vbool", vbool.ptr}}, {{"fvec3", fvec3.ptr}}, {{"cint", cint.ptr}}});

    Triangles m2;
    SurfaceAttributes attrs = read_by_extension(filename, m2);

    REQUIRE( m.nverts()==m2.nverts() );
    REQUIRE( m.nfacets()==m2.nfacets() );

    PointAttribute<bool> vbool2("vbool", attrs, m2);
    for (int v=0; v<m.nverts(); v++) {
        REQUIRE( (m.points[v]-m2.points[v]).norm()<ftol );
        REQUIRE( vbool[v]==vbool2[v] );
    }

    FacetAttribute<vec3> fvec32("fvec3", attrs, m2);
    for (int f=0; f<m.nfacets(); f++) {
        for (int lv=0; lv<m.facet_size(f); lv++) {
            REQUIRE( m.vert(f, lv)==m2.vert(f, lv) );
        }
        REQUIRE( (fvec3[f]-fvec32[f]).norm()<ftol );
    }
    CornerAttribute<int> cint2("cint", attrs, m2);
    for (int c=0; c<m.ncorners(); c++) {
        REQUIRE( cint[c]==cint2[c] );
    }
}

TEST_CASE("Quads + SurfaceAttributes IO test", "[geogram]") {
    static const std::string filename = "ultimaille-test-quads.geogram";
    Quads m;
    *m.points.data = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {1,1,1}, {0,1,1}};
    m.facets = {0,1,2,3,3,2,4,5};

    PointAttribute<bool> vbool(m.points);
    FacetAttribute<vec3> fvec3(m);
    CornerAttribute<int> cint(m);

    for (int v : range(m.nverts()))
        vbool[v] = rand()&1;

    for (int f : range(m.nfacets()))
        fvec3[f] = m.util.bary_verts(f);

    for (int c : range(m.ncorners()))
        cint[c] = rand();

    write_by_extension(filename, m, {{{"vbool", vbool.ptr}}, {{"fvec3", fvec3.ptr}}, {{"cint", cint.ptr}}});

    Quads m2;
    SurfaceAttributes attrs = read_by_extension(filename, m2);

    REQUIRE( m.nverts()==m2.nverts() );
    REQUIRE( m.nfacets()==m2.nfacets() );

    PointAttribute<bool> vbool2("vbool", attrs, m2);
    for (int v=0; v<m.nverts(); v++) {
        REQUIRE( (m.points[v]-m2.points[v]).norm()<ftol );
        REQUIRE( vbool[v]==vbool2[v] );
    }

    FacetAttribute<vec3> fvec32("fvec3", attrs, m2);
    for (int f=0; f<m.nfacets(); f++) {
        for (int lv=0; lv<m.facet_size(f); lv++) {
            REQUIRE( m.vert(f, lv)==m2.vert(f, lv) );
        }
        REQUIRE( (fvec3[f]-fvec32[f]).norm()<ftol );
    }
    CornerAttribute<int> cint2("cint", attrs, m2);
    for (int c=0; c<m.ncorners(); c++) {
        REQUIRE( cint[c]==cint2[c] );
    }
}

TEST_CASE("Polygons + SurfaceAttributes IO test", "[geogram]") {
    static const std::string filename = "ultimaille-test-polygons.geogram";
    Polygons m;
    *m.points.data = {{ 1 , 0 ,0}, { 0.309017, 0.951057,0}, {-0.809017, 0.587785,0}, {-0.809017,-0.587785,0}, { 0.309017,-0.951057,0}};
    m.facets = {0,1,2,3,4,0,2};
    m.offset = {0,3,7};

    PointAttribute<bool> vbool(m.points);
    FacetAttribute<vec3> fvec3(m);
    CornerAttribute<int> cint(m);

    for (int v : range(m.nverts()))
        vbool[v] = rand()&1;

    for (int f : range(m.nfacets()))
        fvec3[f] = m.util.bary_verts(f);

    for (int c : range(m.ncorners()))
        cint[c] = rand();

    write_by_extension(filename, m, {{{"vbool", vbool.ptr}}, {{"fvec3", fvec3.ptr}}, {{"cint", cint.ptr}}});

    Polygons m2;
    SurfaceAttributes attrs = read_by_extension(filename, m2);

    REQUIRE( m.nverts()==m2.nverts() );
    REQUIRE( m.nfacets()==m2.nfacets() );

    PointAttribute<bool> vbool2("vbool", attrs, m2);
    for (int v=0; v<m.nverts(); v++) {
        REQUIRE( (m.points[v]-m2.points[v]).norm()<ftol );
        REQUIRE( vbool[v]==vbool2[v] );
    }

    FacetAttribute<vec3> fvec32("fvec3", attrs, m2);
    for (int f=0; f<m.nfacets(); f++) {
        REQUIRE(m.facet_size(f)==m2.facet_size(f));
        for (int lv=0; lv<m.facet_size(f); lv++) {
            REQUIRE( m.vert(f, lv)==m2.vert(f, lv) );
        }
        REQUIRE( (fvec3[f]-fvec32[f]).norm()<ftol );
    }
    CornerAttribute<int> cint2("cint", attrs, m2);
    for (int c=0; c<m.ncorners(); c++) {
        REQUIRE( cint[c]==cint2[c] );
    }
}

TEST_CASE("Tetrahedra + VolumeAttributes IO test", "[geogram]") {
    static const std::string filename = "ultimaille-test-tetrahedra.geogram";
    Tetrahedra m;
    *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,1,1}};
    m.cells = {0,1,2,3,4,3,2,1};

    PointAttribute<bool> vbool(m.points);
    CellAttribute<double> cdouble(m);
    CellFacetAttribute<vec3> cfvec3(m);
    CellCornerAttribute<int> ccint(m);

    for (int v : range(m.nverts()))
        vbool[v] = rand()&1;

    for (int c : range(m.ncells()))
        cdouble[c] = m.util.cell_volume(c);

    for (int cf : range(m.nfacets()))
        cfvec3[cf] = m.util.bary_facet(cf/m.nfacets_per_cell(), cf%m.nfacets_per_cell());

    for (int cc : range(m.ncorners()))
        ccint[cc] = rand();

    write_by_extension(filename, m, {{{"vbool", vbool.ptr}}, {{"cdouble", cdouble.ptr}}, {{"cfvec3", cfvec3.ptr}}, {{"ccint", ccint.ptr}}});

    Tetrahedra m2;
    VolumeAttributes attrs = read_by_extension(filename, m2);

    REQUIRE( m.nverts()==m2.nverts() );
    REQUIRE( m.ncells()==m2.ncells() );

    PointAttribute<bool> vbool2("vbool", attrs, m2);
    for (int v=0; v<m.nverts(); v++) {
        REQUIRE( (m.points[v]-m2.points[v]).norm()<ftol );
        REQUIRE( vbool[v]==vbool2[v] );
    }

    CellAttribute<double> cdouble2("cdouble", attrs, m2);
    for (int c=0; c<m.ncells(); c++) {
        for (int lv=0; lv<m.nverts_per_cell(); lv++) {
            REQUIRE( m.vert(c, lv)==m2.vert(c, lv) );
        }
        REQUIRE( std::fabs(cdouble[c]-cdouble2[c])<ftol );
    }

    CellFacetAttribute<vec3> cfvec32("cfvec3", attrs, m2);
    for (int cf=0; cf<m.nfacets(); cf++) {
        REQUIRE( (cfvec3[cf]-cfvec32[cf]).norm()<ftol );
    }

    CellCornerAttribute<int> ccint2("ccint", attrs, m2);
    for (int cc=0; cc<m.ncorners(); cc++) {
        REQUIRE( ccint[cc]==ccint2[cc] );
    }
}

TEST_CASE("Hexahedra + VolumeAttributes IO test", "[geogram]") {
    static const std::string filename = "ultimaille-test-hexahedra.geogram";
    Hexahedra m;
    *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}, {0,0,2}, {1,0,2}, {0,1,2}, {1,1,2}};
    m.cells = {0,1,2,3,4,5,6,7,4,5,6,7,8,9,10,11};

    PointAttribute<bool> vbool(m.points);
    CellAttribute<double> cdouble(m);
    CellFacetAttribute<vec3> cfvec3(m);
    CellCornerAttribute<int> ccint(m);


    for (int v : range(m.nverts()))
        vbool[v] = rand()&1;

    for (int c : range(m.ncells()))
        cdouble[c] = m.util.cell_volume(c);

    for (int cf : range(m.nfacets()))
        cfvec3[cf] = m.util.bary_facet(cf/m.nfacets_per_cell(), cf%m.nfacets_per_cell());

    for (int cc : range(m.ncorners()))
        ccint[cc] = rand();

    write_by_extension(filename, m, {{{"vbool", vbool.ptr}}, {{"cdouble", cdouble.ptr}}, {{"cfvec3", cfvec3.ptr}}, {{"ccint", ccint.ptr}}});

    Hexahedra m2;
    VolumeAttributes attrs = read_by_extension(filename, m2);

    REQUIRE( m.nverts()==m2.nverts() );
    REQUIRE( m.ncells()==m2.ncells() );

    PointAttribute<bool> vbool2("vbool", attrs, m2);
    for (int v=0; v<m.nverts(); v++) {
        REQUIRE( (m.points[v]-m2.points[v]).norm()<ftol );
        REQUIRE( vbool[v]==vbool2[v] );
    }

    CellAttribute<double> cdouble2("cdouble", attrs, m2);
    for (int c=0; c<m.ncells(); c++) {
        for (int lv=0; lv<m.nverts_per_cell(); lv++) {
            REQUIRE( m.vert(c, lv)==m2.vert(c, lv) );
        }
        REQUIRE( std::fabs(cdouble[c]-cdouble2[c])<ftol );
    }

    CellFacetAttribute<vec3> cfvec32("cfvec3", attrs, m2);
    for (int cf=0; cf<m.nfacets(); cf++) {
        REQUIRE( (cfvec3[cf]-cfvec32[cf]).norm()<ftol );
    }

    CellCornerAttribute<int> ccint2("ccint", attrs, m2);
    for (int cc=0; cc<m.ncorners(); cc++) {
        REQUIRE( ccint[cc]==ccint2[cc] );
    }
}

TEST_CASE("Wedges + VolumeAttributes IO test", "[geogram]") {
    static const std::string filename = "ultimaille-test-wedges.geogram";
    Wedges m;
    *m.points.data = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {0,0,2}, {1,0,2}, {0,1,2}};
    m.cells = {0,1,2,3,4,5,3,4,5,6,7,8};

    PointAttribute<bool> vbool(m.points);
    CellAttribute<double> cdouble(m);
    CellFacetAttribute<vec3> cfvec3(m);
    CellCornerAttribute<int> ccint(m);

    for (int v : range(m.nverts()))
        vbool[v] = rand()&1;

    for (int c : range(m.ncells()))
        cdouble[c] = m.util.cell_volume(c);

    for (int cf : range(m.nfacets()))
        cfvec3[cf] = m.util.bary_facet(cf/m.nfacets_per_cell(), cf%m.nfacets_per_cell());

    for (int cc : range(m.ncorners()))
        ccint[cc] = rand();

    write_by_extension(filename, m, {{{"vbool", vbool.ptr}}, {{"cdouble", cdouble.ptr}}, {{"cfvec3", cfvec3.ptr}}, {{"ccint", ccint.ptr}}});

    Wedges m2;
    VolumeAttributes attrs = read_by_extension(filename, m2);

    REQUIRE( m.nverts()==m2.nverts() );
    REQUIRE( m.ncells()==m2.ncells() );

    PointAttribute<bool> vbool2("vbool", attrs, m2);
    for (int v=0; v<m.nverts(); v++) {
        REQUIRE( (m.points[v]-m2.points[v]).norm()<ftol );
        REQUIRE( vbool[v]==vbool2[v] );
    }

    CellAttribute<double> cdouble2("cdouble", attrs, m2);
    for (int c=0; c<m.ncells(); c++) {
        for (int lv=0; lv<m.nverts_per_cell(); lv++) {
            REQUIRE( m.vert(c, lv)==m2.vert(c, lv) );
        }
        REQUIRE( std::fabs(cdouble[c]-cdouble2[c])<ftol );
    }

    CellFacetAttribute<vec3> cfvec32("cfvec3", attrs, m2);
    for (int cf=0; cf<m.nfacets(); cf++) {
        REQUIRE( (cfvec3[cf]-cfvec32[cf]).norm()<ftol );
    }

    CellCornerAttribute<int> ccint2("ccint", attrs, m2);
    for (int cc=0; cc<m.ncorners(); cc++) {
        REQUIRE( ccint[cc]==ccint2[cc] );
    }
}

TEST_CASE("Pyramids + VolumeAttributes IO test", "[geogram]") {
    static const std::string filename = "ultimaille-test-pyramids.geogram";
    Pyramids m;
    *m.points.data = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {.5,.5,.5}, {.5,.5,-.5}};
    m.cells = {0,1,2,3,4,3,2,1,0,5};

    PointAttribute<bool> vbool(m.points);
    CellAttribute<double> cdouble(m);
    CellFacetAttribute<vec3> cfvec3(m);
    CellCornerAttribute<int> ccint(m);

    for (int v : range(m.nverts()))
        vbool[v] = rand()&1;

    for (int c : range(m.ncells()))
        cdouble[c] = m.util.cell_volume(c);

    for (int cf : range(m.nfacets()))
        cfvec3[cf] = m.util.bary_facet(cf/m.nfacets_per_cell(), cf%m.nfacets_per_cell());

    for (int cc : range(m.ncorners()))
        ccint[cc] = rand();

    write_by_extension(filename, m, {{{"vbool", vbool.ptr}}, {{"cdouble", cdouble.ptr}}, {{"cfvec3", cfvec3.ptr}}, {{"ccint", ccint.ptr}}});

    Pyramids m2;
    VolumeAttributes attrs = read_by_extension(filename, m2);

    REQUIRE( m.nverts()==m2.nverts() );
    REQUIRE( m.ncells()==m2.ncells() );

    PointAttribute<bool> vbool2("vbool", attrs, m2);
    for (int v=0; v<m.nverts(); v++) {
        REQUIRE( (m.points[v]-m2.points[v]).norm()<ftol );
        REQUIRE( vbool[v]==vbool2[v] );
    }

    CellAttribute<double> cdouble2("cdouble", attrs, m2);
    for (int c=0; c<m.ncells(); c++) {
        for (int lv=0; lv<m.nverts_per_cell(); lv++) {
            REQUIRE( m.vert(c, lv)==m2.vert(c, lv) );
        }
        REQUIRE( std::fabs(cdouble[c]-cdouble2[c])<ftol );
    }

    CellFacetAttribute<vec3> cfvec32("cfvec3", attrs, m2);
    for (int cf=0; cf<m.nfacets(); cf++) {
        REQUIRE( (cfvec3[cf]-cfvec32[cf]).norm()<ftol );
    }

    CellCornerAttribute<int> ccint2("ccint", attrs, m2);
    for (int cc=0; cc<m.ncorners(); cc++) {
        REQUIRE( ccint[cc]==ccint2[cc] );
    }
}

