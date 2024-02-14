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
struct line2 {

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

    // Local point type in line coordinates
    using loc_point = point2;

    /// @}

    /** This method transforms a point from a global cartesian 3D frame to a
     * local 3D line point */
    DETRAY_HOST_DEVICE
    static inline point3 global_to_local(const transform3_t &trf,
                                         const point3 &p, const vector3 &d) {

        const auto local3 = trf.point_to_local(p);

        // Line direction
        const vector3 z = trf.z();

        // Line center
        const point3 t = trf.translation();

        // Radial vector
        const vector3 r = vector::cross(z, d);

        // Assign the sign depending on the position w.r.t line
        // Right: -1
        // Left: 1
        const scalar_type sign = vector::dot(r, t - p) > 0.f ? -1.f : 1.f;

        return {sign * getter::perp(local3), local3[2], getter::phi(local3)};
    }

    /** This method transform from a local 2D line point to a point global
     * cartesian 3D frame*/
    DETRAY_HOST_DEVICE static inline point3 local_to_global(
        const transform3_t &trf, const point3 &p) {
        const scalar_type R = math::abs(p[0]);
        const point3 local = {R * math::cos(p[2]), R * math::sin(p[2]), p[1]};

        return trf.point_to_global(local);
    }

    /** This method transform from a local 2D line point to a point global
     * cartesian 3D frame*/
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline point3 bound_local_to_global(
        const transform3_t &trf, const mask_t & /*mask*/, const point2 &p,
        const vector3 &d) {

        // Line direction
        const vector3 z = trf.z();

        // Radial vector
        const vector3 r = vector::cross(z, d);

        // Local Z poisition in global cartesian coordinate
        const point3 locZ_in_global =
            trf.point_to_global(point3{0.f, 0.f, p[1]});

        return locZ_in_global + p[0] * vector::normalize(r);
    }

    /// @returns the normal vector
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_t &trf3,
                                                    const point2 & = {},
                                                    const mask_t & = {}) {
        return trf3.z();
    }

    /// @returns the normal vector
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_t &trf3,
                                                    const point3 & = {}) {
        return trf3.z();
    }
};

}  // namespace detray
