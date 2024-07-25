/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/eigen_eigen.hpp"

#define ALGEBRA_PLUGIN detray::eigen

namespace detray {

// Define scalar type
using scalar = DETRAY_CUSTOM_SCALARTYPE;

/// Define affine transformation types
/// @{
template <typename V = DETRAY_CUSTOM_SCALARTYPE>
struct eigen {
    /// Define scalar type
    using value_type = V;

    template <typename T>
    using simd = T;

    using boolean = bool;
    using scalar = value_type;
    using size_type = algebra::eigen::size_type;
    using transform3D = algebra::eigen::transform3<value_type>;
    using point2D = algebra::eigen::point2<value_type>;
    using point3D = algebra::eigen::point3<value_type>;
    using vector3D = algebra::eigen::vector3<value_type>;

    template <std::size_t ROWS, std::size_t COLS>
    using matrix = algebra::eigen::matrix_type<value_type, ROWS, COLS>;
};
/// @}

namespace getter {

using algebra::eigen::storage::block;
using algebra::eigen::storage::element;
using algebra::eigen::storage::set_block;
using algebra::eigen::storage::vector;

}  // namespace getter

namespace vector {

using algebra::eigen::math::cross;
using algebra::eigen::math::dot;
using algebra::eigen::math::eta;
using algebra::eigen::math::norm;
using algebra::eigen::math::normalize;

using algebra::eigen::math::perp;
using algebra::eigen::math::phi;
using algebra::eigen::math::theta;

}  // namespace vector

namespace matrix {

using algebra::eigen::math::determinant;
using algebra::eigen::math::identity;
using algebra::eigen::math::inverse;
using algebra::eigen::math::set_identity;
using algebra::eigen::math::set_zero;
using algebra::eigen::math::transpose;
using algebra::eigen::math::zero;

}  // namespace matrix

}  // namespace detray
