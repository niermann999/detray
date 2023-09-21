/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray core include(s).
#include "detray/definitions/algebra.hpp"

namespace detray::test {

using scalar = detray::scalar;
template <typename T = test::scalar>
using algebra = ALGEBRA_PLUGIN;
using transform3 = dtransform3D<algebra<scalar>>;
using point2 = dpoint2D<algebra<scalar>>;
using point3 = dpoint3D<algebra<scalar>>;
using vector2 = dvector2D<algebra<scalar>>;
using vector3 = dvector3D<algebra<scalar>>;

#if DETRAY_ALGEBRA_ARRAY

// The std::array based algebra plugin is always available in the tests
template <typename T = test::scalar>
using algebra_t = detray::cmath<T>;
static constexpr char filenames[] = "cmath-";

#elif DETRAY_ALGEBRA_EIGEN

template <typename T = test::scalar>
using algebra_t = detray::eigen<T>;
static constexpr char filenames[] = "eigen-";

#elif DETRAY_ALGEBRA_SMATRIX

template <typename T = test::scalar>
using algebra_t = detray::smatrix<T>;
static constexpr char filenames[] = "smatrix-";

#elif DETRAY_ALGEBRA_VC

template <typename T = test::scalar>
using algebra_t = detray::vc_cmath<T>;
static constexpr char filenames[] = "vc_cmath-";
#endif

}  // namespace detray::test
