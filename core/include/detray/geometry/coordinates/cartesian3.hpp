/** Detray library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray {

/// Frame projection into a cartesian coordinate frame
template <typename transform3_t>
struct cartesian3D final : public local_frame<cartesian3D<algebra_t>> {

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

    // Local point type in 3D cartesian coordinates
    using loc_point = point3;

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
        return p;
    }

    /// Projection from the local frame of a reference surface
    ///
    /// @param p the point to be projected
    ///
    /// @returns the projected point
        template <typename mask_t>
    DETRAY_HOST_DEVICE
    static constexpr point3 project(
        const transform3_t &/*trf*/, const mask_t &/*mask*/, const loc_point &p,
        const vector3 &/*d*/) {
        return p;
    }

    /// @returns the normal vector in global coordinates
    template <typename mask_t>
    DETRAY_HOST_DEVICE static constexpr vector3 normal(const transform3 &trf,
                                                       const point2 & = {},
                                                       const mask_t & = {}) {
        return trf.z();
    }
    /// @}

};  // struct cartesian3

}  // namespace detray
