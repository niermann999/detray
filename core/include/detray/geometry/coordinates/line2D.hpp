/** Detray library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/coordinates/local_frame.hpp"

namespace detray {

template <typename algebra_t>
struct line2D final : public local_frame<line2D, algebra_t> {

    /// @name Type definitions for the struct
    /// @{
    using transform3_type = algebra_t;
    using scalar_type = typename algebra_t::scalar_type;
    using point2 = typename algebra_t::point2;
    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;

    // Local point type in line coordinates
    using local_point = point2;
    /// @}

    /// This method transforms a point from a global cartesian 3D frame to a
    /// local 3D line point
    DETRAY_HOST_DEVICE
    static constexpr local_point global_to_local(const transform3_type &trf,
                                                 const point3 &p,
                                                 const vector3 &d) {

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

        return {sign * getter::perp(local3), local3[2]};
    }

    /// This method transform from a local 2D line point to a point global
    /// cartesian 3D frame
    ///
    /// @note overwrites the base class implementation of @c local_to_global
    template <typename mask_t>
    DETRAY_HOST_DEVICE static constexpr point3 local_to_global(
        const transform3_type &trf, const mask_t & /*mask*/,
        const local_point &loc_p, const vector3 &d) {

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
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_type &trf,
                                                    const point2 & = {},
                                                    const mask_t & = {}) {
        return trf.z();
    }
};

}  // namespace detray
