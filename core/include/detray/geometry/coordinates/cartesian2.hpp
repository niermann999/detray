/** Detray library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray {

/** Frame projection into a cartesian coordinate frame
 */
template <typename transform3_t>
struct cartesian2 {

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

    // Local point type in 2D cartesian coordinates
    using loc_point = point2;

    /// @}

    /// This method transforms a point from a global cartesian 3D frame to a
    /// local 3D cartesian point
    DETRAY_HOST_DEVICE
    static inline point3 global_to_local(const transform3_t &trf3,
                                         const point3 &p,
                                         const vector3 & /*d*/) {
        return trf3.point_to_local(p);
    }

    /** This method transforms a point from a global cartesian 3D frame to a
     * bound 2D cartesian point */
    DETRAY_HOST_DEVICE
    static inline loc_point project_to_axes(const transform3_t &trf3,
                                            const point3 &p,
                                            const vector3 & /*d*/) {
        auto loc_p = trf3.point_to_local(p);
        return {loc_p[0], loc_p[1]};
    }

    /** This method transform from a local 2D cartesian point to a point global
     * cartesian 3D frame*/
    DETRAY_HOST_DEVICE static inline point3 local_to_global(
        const transform3_t &trf3, const point3 &p) {
        return trf3.point_to_global(p);
    }

    /** This method transform from a local 2D cartesian point to a point global
     * cartesian 3D frame*/
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline point3 bound_local_to_global(
        const transform3_t &trf3, const mask_t & /*mask*/, const point2 &p,
        const vector3 & /*d*/) {

        return cartesian2<transform3_t>::local_to_global(trf3,
                                                         {p[0], p[1], 0.f});
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

};  // struct cartesian2

}  // namespace detray
