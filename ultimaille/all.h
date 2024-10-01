#include <ultimaille/algebra/vec.h>
#include <ultimaille/algebra/mat.h>
#include <ultimaille/algebra/quaternion.h>
#include <ultimaille/algebra/eigen.h>
#include <ultimaille/algebra/svd.h>
#include <ultimaille/algebra/covariance.h>

#include <ultimaille/sparse/vector.h>
#include <ultimaille/sparse/matrix.h>
#include <ultimaille/sparse/nullspace.h>
#include <ultimaille/sparse/linexpr.h>
#include <ultimaille/sparse/least_squares.h>

#include <ultimaille/helpers/colocate.h>
#include <ultimaille/helpers/constraints.h>
#include <ultimaille/helpers/disjointset.h>
#include <ultimaille/helpers/permutation.h>
#include <ultimaille/helpers/hilbert_sort.h>
#include <ultimaille/helpers/hboxes.h>
#include <ultimaille/helpers/knn.h>
#include <ultimaille/helpers/bvh.h>

#include <ultimaille/meter.h>

#include <ultimaille/io/xyz.h>
#include <ultimaille/io/obj.h>
#include <ultimaille/io/geogram.h>
#include <ultimaille/io/medit.h>
#include <ultimaille/io/vtk.h>
#include <ultimaille/io/by_extension.h>

#include <ultimaille/syntactic-sugar/range.h>
#include <ultimaille/syntactic-sugar/parallel.h>
#include <ultimaille/syntactic-sugar/assert.h>

#include <ultimaille/attributes.h>
#include <ultimaille/pointset.h>
#include <ultimaille/polyline.h>
#include <ultimaille/surface.h>
#include <ultimaille/volume.h>
#include <ultimaille/volume_connectivity.h>
#include <ultimaille/primitive_geometry.h>
#include <ultimaille/attr_binding.h>

#include <stlbfgs.h>
#include <OpenNL_psm/OpenNL_psm.h>

