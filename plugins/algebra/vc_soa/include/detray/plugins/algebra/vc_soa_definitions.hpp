
/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Algebra-Plugins include
#include "algebra/vc_soa.hpp"

#define IS_SOA 1

namespace detray {

using algebra::storage::operator*;
using algebra::storage::operator/;
using algebra::storage::operator-;
using algebra::storage::operator+;

using scalar = DETRAY_CUSTOM_SCALARTYPE;

/// Define affine transformation types
/// @{
template <typename V = DETRAY_CUSTOM_SCALARTYPE>
struct vc_soa {
    /// Define scalar precision
    using value_type = V;

    template <typename T>
    using simd = Vc::Vector<T>;

    using boolean = Vc::Mask<V>;

    /// Linear Algebra type definitions
    /// @{
    using scalar = simd<value_type>;
    using size_type = algebra::vc_soa::size_type;
    using transform3D = algebra::vc_soa::transform3<value_type>;
    using point2D = algebra::vc_soa::point2<value_type>;
    using point3D = algebra::vc_soa::point3<value_type>;
    using vector3D = algebra::vc_soa::vector3<value_type>;

    template <std::size_t ROWS, std::size_t COLS>
    using matrix = algebra::vc_soa::matrix_type<value_type, ROWS, COLS>;
    /// @}
};
/// @}

namespace getter {

using algebra::vc_soa::storage::block;
using algebra::vc_soa::storage::element;
using algebra::vc_soa::storage::set_block;
using algebra::vc_soa::storage::vector;

}  // namespace getter

namespace vector {

using algebra::vc_soa::math::cross;
using algebra::vc_soa::math::dot;
using algebra::vc_soa::math::eta;
using algebra::vc_soa::math::norm;
using algebra::vc_soa::math::normalize;
using algebra::vc_soa::math::perp;
using algebra::vc_soa::math::phi;
using algebra::vc_soa::math::theta;

}  // namespace vector

namespace matrix {

using algebra::vc_soa::math::determinant;
using algebra::vc_soa::math::identity;
using algebra::vc_soa::math::inverse;
using algebra::vc_soa::math::set_identity;
using algebra::vc_soa::math::set_zero;
using algebra::vc_soa::math::transpose;
using algebra::vc_soa::math::zero;

}  // namespace matrix

}  // namespace detray
