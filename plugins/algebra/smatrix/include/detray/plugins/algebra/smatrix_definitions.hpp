/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// TODO: Remove this when gcc fixes their false positives.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#endif

// Algebra-Plugins include
#include "algebra/smatrix_smatrix.hpp"

#define ALGEBRA_PLUGIN detray::smatrix

namespace detray {

// Define scalar type
using scalar = DETRAY_CUSTOM_SCALARTYPE;

/// Define affine transformation types
/// @{
template <typename V = DETRAY_CUSTOM_SCALARTYPE>
struct smatrix {
    /// Define scalar type
    using value_type = V;

    template <typename T>
    using simd = T;

    using boolean = bool;
    using scalar = value_type;
    using size_type = algebra::smatrix::size_type;
    using transform3D = algebra::smatrix::transform3<value_type>;
    using point2D = algebra::smatrix::point2<value_type>;
    using point3D = algebra::smatrix::point3<value_type>;
    using vector3D = algebra::smatrix::vector3<value_type>;

    template <std::size_t ROWS, std::size_t COLS>
    using matrix = algebra::smatrix::matrix_type<value_type, ROWS, COLS>;
};
/// @}

namespace getter {

using algebra::smatrix::storage::block;
using algebra::smatrix::storage::element;
using algebra::smatrix::storage::set_block;
using algebra::smatrix::storage::vector;
}  // namespace getter

namespace vector {

using algebra::smatrix::math::cross;
using algebra::smatrix::math::dot;
using algebra::smatrix::math::eta;
using algebra::smatrix::math::norm;
using algebra::smatrix::math::normalize;
using algebra::smatrix::math::perp;
using algebra::smatrix::math::phi;
using algebra::smatrix::math::theta;

}  // namespace vector

namespace matrix {

using algebra::smatrix::math::determinant;
using algebra::smatrix::math::identity;
using algebra::smatrix::math::inverse;
using algebra::smatrix::math::set_identity;
using algebra::smatrix::math::set_zero;
using algebra::smatrix::math::transpose;
using algebra::smatrix::math::zero;

}  // namespace matrix

}  // namespace detray
