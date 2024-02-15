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

template <typename algebra_t>
struct line2D final : public local_frame<line2D<algebra_t>> {

    /// @name Type definitions for the struct
    /// @{

    // Transform type
    using transform3_type = algebra_t;
    // Sclar type
    using scalar_t = typename algebra_t::scalar_type;
    // Point in 2D space
    using point2 = typename algebra_t::point2;
    // Point in 3D space
    using point3 = typename algebra_t::point3;
    // Vector in 3D space
    using vector3 = typename algebra_t::vector3;
    // Local point type in line coordinates
    using loc_point = point2;

    /// @}

    /// This method transforms a point from a global cartesian 3D frame to a
    /// local 3D line point
    DETRAY_HOST_DEVICE
    static constexpr loc_point project(const transform3 &trf,
                                       const point3 &loc_p3,
                                       const vector3 &d) {
        // Line direction
        const vector3 z = trf.z();

        // Line center
        const point3 t = trf.translation();

        // Radial vector
        const vector3 r = vector::cross(z, d);

        // Assign the sign depending on the position w.r.t line
        // Right: -1
        // Left: 1
        const scalar_t sign = vector::dot(r, t - p) > 0.f ? -1.f : 1.f;

        return {sign * getter::perp(loc_p3), loc_p3[2]};
    }

    /// This method transform from a local 2D line point to a point global
    /// cartesian 3D frame
    ///
    /// @note overwrites the base class implementation of @c local_to_global
    DETRAY_HOST_DEVICE
    static constexpr point3 local_to_global(
        const transform3_t &trf, const mask_t &/*mask*/, const loc_point &loc_p,
        const vector3 &d) {

        // Line direction
        const vector3 z = trf.z();

        // Radial vector
        const vector3 r = vector::cross(z, d);

        // Local Z poisition in global cartesian coordinate
        const point3 locZ_in_global =
            trf.point_to_global(point3{0.f, 0.f, loc_p[1]});

        return locZ_in_global + loc_p[0] * vector::normalize(r);
    }

    /// @returns the normal vector
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_t &trf3,
                                                    const point2 & = {},
                                                    const mask_t & = {}) {
        return trf.z();
    }
};

}  // namespace detray
