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
struct cylindrical2D final : public local_frame<cylindrical2D<algebra_t>> {

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

    // Local point type in 2D cylinderical coordinates
    using loc_point = point2;

    /// @}

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
        return {getter::perp(p) * getter::phi(p), p[2]};
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
        const transform3_t &/*trf*/, const mask_t &mask, const point2 &p,
        const vector3 &/*d*/) {

        const scalar_type r{mask[mask_t::shape::e_r]};
        const scalar_type phi{p[0] / r};

        return {r * math::cos(phi), r * math::sin(phi), p[1]};
    }

    /// @returns the normal vector in global coordinates
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_t &trf,
                                                    const point2 &p,
                                                    const mask_t &mask) {
        const scalar_type phi{p[0] / mask[mask_t::shape::e_r]};
        const vector3 local_normal{math::cos(phi), math::sin(phi), 0.f};

        return trf.vector_to_global(local_normal);
    }

};  // struct cylindrical2

}  // namespace detray
