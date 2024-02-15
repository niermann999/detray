/** Detray library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray {

/// Frame projection into a cartesian coordinate frame
template <typename algebra_t>
struct cartesian2D final : public local_frame<cartesian2D<algebra_t>> {

    /// @name Linear algebra types
    /// @{

    // Transform type
    using transform3_type = algebra_t;
    // Sclar type
    using scalar_type = typename algebra_t::scalar_type;
    // Point in 2D space
    using point2 = typename algebra_t::point2;
    // Point in 3D space
    using point3 = typename algebra_t::point3;
    // Vector in 3D space
    using vector3 = typename algebra_t::vector3;

    // Local point type in 2D cartesian coordinates
    using local_point = point2;

    /// @}

    /// Projection into the local frame of a reference surface
    ///
    /// @param trf3 the transformation of the reference surface
    /// @param p the point to be projected
    /// @param d an optional orientation
    ///
    /// @note the projected point is contrained to the local frame and not
    /// meaningful without access to the surface
    ///
    /// @returns the projected point
    DETRAY_HOST_DEVICE
    static constexpr loc_point project(const transform3 &/*trf*/,
                                       const point3 &p,
                                       const vector3 &/*d*/) {
        return {p[0], p[1]};
    }

    /// Projection into the local frame of a reference surface
    ///
    /// @param trf3 the transformation of the reference surface
    /// @param mask the mask of the reference surface
    /// @param p the constrained point
    /// @param d an optional orientation
    ///
    /// @returns the point in the local 3D cartesian frame
        template <typename mask_t>
    DETRAY_HOST_DEVICE
    static constexpr point3 project(
        const transform3_t &/*trf*/, const mask_t &/*mask*/, const loc_point &p,
        const vector3 &/*d*/) {
        return {p[0], p[1], 0.f};
    }

    /// @returns the normal vector in global coordinates
    template <typename mask_t>
    DETRAY_HOST_DEVICE static constexpr vector3 normal(const transform3 &trf3,
                                                       const loc_point & = {},
                                                       const mask_t & = {}) {
        return trf3.z();
    }

};  // struct cartesian2

}  // namespace detray
