/** Detray library, part of the ACTS project
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray {

template <typename transform3_t>
struct polar2D final : public local_frame<polar2D<algebra_t>> {

    /// @name Type definitions for the struct
    /// @{

    // Transform type
    using transform3_type = transform3_t;
    // Sclar type
    using scalar_type = typename transform3_type::scalar_type;
    // Point in 2D space
    using point2 = typename transform3_type::point2;
    // Point in 3D space
    using point3 = typename transform3_type::point3;
    // Vector in 3D space
    using vector3 = typename transform3_type::vector3;

    // Local point type in polar coordinates
    using loc_point = point2;


    /// @brief Tranform from the global frame to the local frame
    ///
    /// @param trf3 the transformation to the local frame
    /// @param p the point to be transformed
    ///
    /// @returns a point in the local 3D cartesian frame
    DETRAY_HOST_DEVICE
    static constexpr loc_point project(const transform3 &/*trf*/,
                                       const point3 &p,
                                       const vector3 &/*d*/) {
        return {getter::perp(p), getter::phi(p)};
    }

    /// Transforms from a local 2D cartesian point that is constrained to a
    /// reference surface to a point in the global 3D cartesian frame
    ///
    /// @param trf placement of the surface the point is contrained to
    /// @param p the contrained point
    ///
    /// @returns the point in global coordinates
    template <typename mask_t>
    DETRAY_HOST_DEVICE
    static constexpr point3 project(
        const transform3_t &/*trf*/, const mask_t &/*mask*/, const loc_point &p,
        const vector3 &/*d*/) {

        return {p[0] * math::cos(p[1]), p[0] * math::sin(p[1]), 0.f};
    }

    /// @returns the normal vector in global coordinates
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_t &trf,
                                                    const point2 & = {},
                                                    const mask_t & = {}) {
        return trf.z();
    }
};

}  // namespace detray
