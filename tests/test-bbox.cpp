#define _USE_MATH_DEFINES
#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <iostream>
#include <ultimaille/all.h>

using namespace UM;

TEST_CASE("Test 1D bbox", "[bb]") {

    BBox1 bbox1({0.}, {1.});
    //BBox<2> b;
    CHECK_FALSE(bbox1.empty());

    // Check point is contained in bbox
    CHECK(bbox1.contains({.5}));

    // Check that bbox2 overlap bbox1
    BBox1 b2({0.8}, {1.2});
    CHECK(bbox1.intersect(b2));

    // Check that center is approximatively equal to {0.5}
    CHECK(std::abs((bbox1.center()-vec<1>{0.5}).norm()) < 1e-4);
    // Check that size is approximatively equal to {1.}
    CHECK(std::abs((bbox1.size() - vec<1>{1.}).norm()) < 1e-4);

    // Size before dilate
    double p_min = bbox1.min.data[0];
    double p_max = bbox1.max.data[0];
    // Dilate BBox
    bbox1.dilate(1.);
    // Check new size of bbox is correct
    CHECK(bbox1.min[0] <= p_min - 1.);
    CHECK(bbox1.max[0] >= p_max + 1.);
    

    BBox1 empty_bbox;
    CHECK(empty_bbox.empty());

}

TEST_CASE("Test 2D bbox", "[bb]") {

    BBox2 bbox1({0.,0.}, {1.,1.});
    //BBox<2> b;
    CHECK_FALSE(bbox1.empty());

    // Check point is contained in bbox
    CHECK(bbox1.contains({.5,.5}));

    // Check that bbox2 overlap bbox1
    BBox2 b2({0.8,0.8}, {1.2,1.2});
    CHECK(bbox1.intersect(b2));

    // Check that center is approximatively equal to {0.5, 0.5}
    CHECK(std::abs((bbox1.center()- vec2{0.5, 0.5}).norm()) < 1e-4);
    // Check that size is approximatively equal to {1., 1.}
    CHECK(std::abs((bbox1.size() - vec2{1., 1.}).norm()) < 1e-4);

    vec2 v;
//  std::cout <<v.x << std::endl;

    // Dilate BBox
    bbox1.dilate(1.);
    // Check new size of bbox is correct
    CHECK(bbox1.min[0] <= -1.);
    CHECK(bbox1.max[0] >= 2.);
    

    BBox2 empty_bbox;
    CHECK(empty_bbox.empty());

}

TEST_CASE("Test 1D hbbox", "[hb]") {

    BBox1 bbox1({0.}, {1.});
    BBox1 bbox2({1.1}, {2.});
    BBox1 bbox3({2.1}, {3.});
    std::vector<BBox1> bboxes;
    bboxes.push_back(bbox1);
    bboxes.push_back(bbox2);
    bboxes.push_back(bbox3);

    HBoxes1 hbbox(bboxes);
    std::vector<int> primitives;

    BBox1 bbox4({1.9}, {2.5});

    hbbox.intersect(bbox4, primitives);
    CHECK(primitives[1]==1);
    

}

TEST_CASE("Test 2D hbbox", "[hb]") {

    BBox2 bbox1({0., 0.}, {1., 1.});
    BBox2 bbox2({1.1, 1.1}, {2., 2.});
    BBox2 bbox3({2.1, 2.1}, {3., 3.});
    std::vector<BBox2> bboxes;
    bboxes.push_back(bbox1);
    bboxes.push_back(bbox2);
    bboxes.push_back(bbox3);

    HBoxes2 hbbox(bboxes);
    std::vector<int> primitives;

    BBox2 bbox4({1.9, 1.9}, {2.5, 2.5});

    hbbox.intersect(bbox4, primitives);
    CHECK(primitives[1]==1);

}
