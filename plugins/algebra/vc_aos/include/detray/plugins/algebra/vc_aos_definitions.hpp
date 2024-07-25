/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/vc_aos.hpp"

#define ALGEBRA_PLUGIN detray::vc_aos

namespace detray {

/// Define scalar type
using scalar = DETRAY_CUSTOM_SCALARTYPE;

/// Define affine transformation types
/// @{
template <typename V = DETRAY_CUSTOM_SCALARTYPE>
struct vc_aos {
    /// Define scalar type
    using value_type = V;

    template <typename T>
    using simd = T;

    using boolean = bool;
    using scalar = value_type;
    using size_type = algebra::vc_aos::size_type;
    using transform3D = algebra::vc_aos::transform3<value_type>;
    using point2D = algebra::vc_aos::point2<value_type>;
    using point3D = algebra::vc_aos::point3<value_type>;
    using vector3D = algebra::vc_aos::vector3<value_type>;

    template <std::size_t ROWS, std::size_t COLS>
    using matrix = algebra::vc_aos::matrix_type<value_type, ROWS, COLS>;
};
/// @}

namespace getter {

using algebra::vc_aos::storage::block;
using algebra::vc_aos::storage::element;
using algebra::vc_aos::storage::set_block;
using algebra::vc_aos::storage::vector;

}  // namespace getter

namespace vector {

// Vc array specific
using algebra::vc_aos::math::cross;
using algebra::vc_aos::math::dot;
using algebra::vc_aos::math::eta;
using algebra::vc_aos::math::norm;
using algebra::vc_aos::math::normalize;

// No specific vectorized implementation needed
using algebra::vc_aos::math::perp;
using algebra::vc_aos::math::phi;
using algebra::vc_aos::math::theta;

}  // namespace vector

namespace matrix {

using algebra::vc_aos::math::identity;
using algebra::vc_aos::math::set_identity;
using algebra::vc_aos::math::set_zero;
using algebra::vc_aos::math::transpose;
using algebra::vc_aos::math::zero;

// Placeholder, until vectorization-friendly version is available
using algebra::vc_aos::math::determinant;
using algebra::vc_aos::math::inverse;

}  // namespace matrix

}  // namespace detray
