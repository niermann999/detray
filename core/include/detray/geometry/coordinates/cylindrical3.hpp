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

/// Frame projection into a 3D cylindrical coordinate frame
template <typename transform3_t>
struct cylindrical3 {

    /// @name Type definitions for the struct
    /// @{

    using transform3_type = transform3_t;
    /// Sclar type
    using scalar_type = typename transform3_type::scalar_type;
    /// Point in 2D space
    using point2 = typename transform3_type::point2;
    /// Point in 3D space
    using point3 = typename transform3_type::point3;
    /// Vector in 3D space
    using vector3 = typename transform3_type::vector3;

    // Local point type in 3D cylindrical coordinates
    using loc_point = point3;

    /// @}

    /** This method transforms a point from a global cartesian 3D frame to a
     * local 3D cylindrical point */
    DETRAY_HOST_DEVICE
    static inline point3 global_to_local(const transform3_t &trf,
                                         const point3 &p,
                                         const vector3 & /*d*/) {
        const auto local3 = trf.point_to_local(p);
        return {getter::perp(local3), getter::phi(local3), local3[2]};
    }

    /** This method transforms a point from a global cartesian 3D frame to a
     * bound 3D cylindrical point */
    DETRAY_HOST_DEVICE
    static inline loc_point project_to_axes(const transform3_t &trf,
                                            const point3 &p, const vector3 &d) {
        return cylindrical3<transform3_t>::global_to_local(trf, p, d);
    }

    /** This method transform from a local 3D cylindrical point to a point
     * global cartesian 3D frame*/
    DETRAY_HOST_DEVICE static inline point3 local_to_global(
        const transform3_t &trf, const point3 &p) {
        const scalar_type x{p[0] * math::cos(p[1])};
        const scalar_type y{p[0] * math::sin(p[1])};

        return trf.point_to_global(point3{x, y, p[2]});
    }

};  // struct cylindrical3

}  // namespace detray
