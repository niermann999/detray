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

/// Local frame projection into a 2D cylindrical coordinate frame
template <typename transform3_t>
struct cylindrical2 {

    /// @name Type definitions for the struct
    /// @{

    // Transform type
    using transform3_type = transform3_t;
    // Sclar type
    using scalar_type = typename transform3_t::scalar_type;
    // Point in 2D space
    using point2 = typename transform3_t::point2;
    // Point in 3D space
    using point3 = typename transform3_t::point3;
    // Vector in 3D space
    using vector3 = typename transform3_t::vector3;

    // Local point type in 2D cylindrical coordinates
    using loc_point = point2;

    /// @}

    /** This method transforms a point from a global cartesian 3D frame to a
     * local 2D cylindrical point */
    DETRAY_HOST_DEVICE
    static inline point3 global_to_local(const transform3_t &trf,
                                         const point3 &p,
                                         const vector3 & /*d*/) {
        const auto local3 = trf.point_to_local(p);

        return {getter::perp(local3) * getter::phi(local3), local3[2],
                getter::perp(local3)};
    }

    /// @returns the projection of a global position onto the grid axis
    /// @note Does not scale with r, due to projection effects in the navigation
    DETRAY_HOST_DEVICE
    static inline loc_point project_to_axes(const transform3_t &trf,
                                            const point3 &p,
                                            const vector3 & /*d*/) {
        const auto local3 = trf.point_to_local(p);

        return {getter::perp(local3) * getter::phi(local3), local3[2]};
    }

    /** This method transform from a local 2D cylindrical point to a point
     * global cartesian 3D frame*/
    DETRAY_HOST_DEVICE static inline point3 local_to_global(
        const transform3_t &trf, const point3 &p) {
        const scalar_type r{p[2]};
        const scalar_type phi{p[0] / r};
        const scalar_type x{r * math::cos(phi)};
        const scalar_type y{r * math::sin(phi)};
        const scalar_type z{p[1]};

        return trf.point_to_global(point3{x, y, z});
    }

    /** This method transform from a local 2D cylindrical point to a point
     * global cartesian 3D frame*/
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline point3 bound_local_to_global(
        const transform3_t &trf, const mask_t &mask, const point2 &p,
        const vector3 & /*dir*/) {

        return cylindrical2<transform3_t>::local_to_global(
            trf, {p[0], p[1], mask[mask_t::shape::e_r]});
    }

    /// @returns the normal vector
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_t &trf3,
                                                    const point2 &bound_pos,
                                                    const mask_t &mask) {
        const scalar_type phi{bound_pos[0] / mask[mask_t::shape::e_r]};
        const vector3 local_normal{math::cos(phi), math::sin(phi), 0.f};

        // normal vector in local coordinate
        return trf3.rotation() * local_normal;
    }

    /// @returns the normal vector given a local position @param loc_pos
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_t &trf3,
                                                    const point3 &loc_pos) {
        const scalar_type phi{loc_pos[0] / loc_pos[2]};
        const vector3 local_normal{math::cos(phi), math::sin(phi), 0.f};

        // normal vector in local coordinate
        return trf3.rotation() * local_normal;
    }

};  // struct cylindrical2

}  // namespace detray
