# UltiMaille: the ultimate mesh processing library
There are lots of mesh processing libraries in the wild, excellent specimens are:
* [geogram](http://alice.loria.fr/software/geogram/doc/html/index.html)
* [libigl](https://github.com/libigl/libigl)
* [pmp](http://www.pmp-library.org/)
* [CGAL](https://www.cgal.org/)

I am, however, not satisfied with either of those. At the moment of this writing, Geogram, for instance, has 847 thousands (sic!) of lines of code.
I'll try to make a library under 10K loc able to do common geometry processing tasks for surfaces and volumes (polygonal and tetrahedral/hexahedral meshes).
Another reason to create yet-another-mesh-library is that I like explicit types and avoid `auto` as long as it reasonable.
In practice it means that I cannot use libigl for the simple reason that I do not know what this data represents:
```C++
Eigen::MatrixXd V;
Eigen::MatrixXi F;
```
Is it a polygonal surface or a tetrahedral mesh? If surface, is it triangulated or is it a generic polygonal mesh? I simply can not tell...


# Common principles
* This library is meant to have a *reasonable* performance. It means that I strive to make it as rapid as possible as long as it does not deteriorate the readability of the source code.
* All the library is built around STL containers (mainly `std::vector<int>`), normally it will have no `malloc/free/new/delete` instructions.
* Despite that, there will be no `size_t` and `iterator::` in the code. An int is an int. Period.
* There will be as little templates as it is reasonable, the default data types are `int` and `double`.

# Compile and run:
```sh
git clone https://github.com/ssloy/ultimaille.git
cd ultimaille
mkdir build
cd build
cmake ..
make
bin/test
```

