/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/array_cmath.hpp"

#define ALGEBRA_PLUGIN detray::cmath

namespace detray {

// Define scalar type
using scalar = DETRAY_CUSTOM_SCALARTYPE;

/// Define affine transformation types
/// @{
template <typename V = DETRAY_CUSTOM_SCALARTYPE>
struct cmath {
    /// Define scalar type
    using value_type = V;

    template <typename T>
    using simd = T;

    using boolean = bool;
    using scalar = value_type;
    using size_type = algebra::array::size_type;
    using transform3D = algebra::array::transform3<value_type>;
    using point2D = algebra::array::point2<value_type>;
    using point3D = algebra::array::point3<value_type>;
    using vector3D = algebra::array::vector3<value_type>;

    template <std::size_t ROWS, std::size_t COLS>
    using matrix = algebra::array::matrix_type<value_type, ROWS, COLS>;
};
/// @}

namespace getter {

using algebra::cmath::storage::block;
using algebra::cmath::storage::element;
using algebra::cmath::storage::set_block;
using algebra::cmath::storage::vector;

}  // namespace getter

namespace vector {

// array specific implementations
using algebra::cmath::dot;
using algebra::cmath::normalize;

// generic implementations
using algebra::cmath::cross;
using algebra::cmath::eta;
using algebra::cmath::norm;
using algebra::cmath::perp;
using algebra::cmath::phi;
using algebra::cmath::theta;

}  // namespace vector

namespace matrix {

using algebra::cmath::identity;
using algebra::cmath::set_identity;
using algebra::cmath::set_zero;
using algebra::cmath::zero;

using algebra::cmath::determinant;
using algebra::cmath::inverse;
using algebra::cmath::transpose;

}  // namespace matrix

}  // namespace detray
