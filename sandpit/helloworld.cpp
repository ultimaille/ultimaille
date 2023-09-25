#include <iostream>
#include <ultimaille/all.h>

using namespace UM;

int main(int argc, char** argv) {
    const int n = 32;
    Triangles m;
    m.points.create_points(n+1);
    m.points[0] = {0,0,0};
    for (int i=0; i<n; i++) {
        double a = i/(double)n * 2*M_PI;
        m.points[i+1] = vec2{std::cos(a), std::sin(a)}.xy0();
    }

#if 0
    for (double t=0; t<5; t+=1e-5) {
        for (int i=0; i<n/2; i++) {
        m.points[i*2+1].z = t;
        m.points[i*2+2].z = -t;
        }
        if (std::abs((m.points[0]-m.points[1]).norm()-(m.points[1]-m.points[2]).norm()) < 1e-5) break;
    }
#endif

    m.create_facets(n);
    for (int i=0; i<n; i++) {
        m.vert(i, 0) = 0;
        m.vert(i, 1) = i+1;
        m.vert(i, 2) = (i+1)%n+1;
    }

    PointAttribute<vec2> tex_coord(m);
    for (int i=0; i<=n; i++)
        tex_coord[i] = m.points[i].xy();

    tex_coord[2] = {.7,-.04};

    m.connect();



    FacetAttribute<mat<3,2>> reference(m);
    FacetAttribute<double> area(m);
    PointAttribute<bool> lock(m);

    for (auto v : m.iter_vertices())
        lock[v] = v.on_boundary();

    std::vector<double> X(2*m.nverts(), 0.); // optimization variables

    for (auto v : m.iter_vertices()) // set up initialization
        for (int d : {0, 1})
            X[2*v+d] = tex_coord[v][d];

    for (auto t : m.iter_facets()) {
        area[t] = m.util.unsigned_area(t);

        vec2 A,B,C;
        m.util.project(t, A, B, C);

        mat<2,2> ST = {{B-A, C-A}};
        ST = ST;
        reference[t] = mat<3,2>{{ {-1,-1},{1,0},{0,1} }}*ST.invert_transpose();
    }

    const auto getJ = [&reference, &m](const std::vector<double>& X, int t)->mat<2, 2> { // get Jacobian matrix for triangle t
        mat<2, 2> J = {};
        for (int i : {0, 1, 2})
            for (int d : {0, 1})
                J[d] += reference[t][i] * X[m.vert(t, i)*2 + d];
        return J;
    };

    std::vector<SpinLock> spin_locks(X.size());

    const STLBFGS::func_grad_eval energy = [&](const std::vector<double>& X, double& F, std::vector<double>& G) {
        F = 0;
        std::fill(G.begin(), G.end(), 0);

#pragma omp parallel for reduction(+:F)
        for (int t=0; t<m.nfacets(); t++) {
            const mat<2, 2> J = getJ(X, t);

            const double f = J.sumsqr()/2.;

            F += f * area[t];

            for (int d : {0, 1}) {
                const vec2& dfda = J[d];
                for (int i : {0, 1, 2}) {
                    const int v = m.vert(t,i);

                    if (lock[v]) continue;
                    spin_locks[v*2+d].lock();

                    G[v*2+d] += (dfda*reference[t][i])*area[t];
                    spin_locks[v*2+d].unlock();
                }
            }
        }
        std::cerr << F << std::endl;
    };

    STLBFGS::Optimizer opt{energy};
    opt.run(X);

    for (int v : m.iter_vertices())
        tex_coord[v] = {X[2*v+0], X[2*v+1]};


     double F;
     std::vector<double> G(X.size());
     energy(X, F, G);
     std::cerr << F << std::endl;

    write_by_extension("l.obj", m, {{{"tex_coord", tex_coord.ptr}}, {}, {}});


    return 0;
}

