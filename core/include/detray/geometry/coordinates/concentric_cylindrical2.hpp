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

namespace detray {

/// Local frame projection into a 2D concentric cylindrical coordinate frame
/// (No rotation in coordinate transformation)
template <typename transform3_t>
struct concentric_cylindrical2D final : public local_frame<concentric_cylindrical2D<algebra_t>> {

    /// @name Type definitions for the struct
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

    // Local point type in 2D cylindrical coordinates
    using loc_point = point2;

    /// @}

    /// This method transforms a point from a global cartesian 3D frame to a
    /// local 2D cylindrical point
    DETRAY_HOST_DEVICE
    static constexpr loc_point project(const transform3 &/*trf*/,
                                       const point3 &p,
                                       const vector3 &/*d*/) {
        return {getter::phi(local3), local3[2]};
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
    static constexpr loc_point project(
        const transform3_t &/*trf*/, const mask_t &mask, const point2 &p,
        const vector3 &/*d*/) {

        const scalar r{mask[mask_t::shape::e_r]};

        return {r * math::cos(p[0]), r * math::sin(p[0]), p[1]};
    }

    /// @returns the normal vector in global coordinates
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_t &trf,
                                                    const loc_point &p,
                                                    const mask_t &mask) {
        const scalar_type phi{p[0]};
        const vector3 local_normal{math::cos(phi), math::sin(phi), 0.f};

        return trf.vector_to_global(local_normal);
    }

};  // struct concentric_cylindrical2

}  // namespace detray
