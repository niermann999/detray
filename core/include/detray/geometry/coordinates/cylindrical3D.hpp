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

/// Frame projection into a 3D cylindrical coordinate frame
template <typename algebra_t>
struct cylindrical3D final : public local_frame<cylindrical3D, algebra_t> {

    /// @name Type definitions for the struct
    /// @{
    using transform3_type = algebra_t;
    using scalar_type = typename algebra_t::scalar_type;
    using point2 = typename algebra_t::point2;
    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;

    // Local point type in 3D cylinderical coordinates
    using local_point = point3;
    /// @}

    /// @returns the point projected onto the reference surface
    DETRAY_HOST_DEVICE
    static constexpr local_point to(const transform3_type & /*trf*/,
                                    const point3 &p, const vector3 & /*d*/) {
        return {getter::perp(p), getter::phi(p), p[2]};
    }

    /// @returns the point in the local 3D cartesian frame
    template <typename mask_t>
    DETRAY_HOST_DEVICE static constexpr point3 from(
        const transform3_type & /*trf*/, const mask_t & /*mask*/,
        const local_point &p, const vector3 & /*d*/) {

        return {p[0] * math::cos(p[1]), p[0] * math::sin(p[1]), p[2]};
    }

    /// @returns the normal vector in global coordinates
    template <typename mask_t>
    DETRAY_HOST_DEVICE static inline vector3 normal(const transform3_type &trf,
                                                    const point2 &p,
                                                    const mask_t & /*mask*/) {

        const vector3 local_normal{math::cos(p[1]), math::sin(p[1]), 0.f};

        return trf.vector_to_global(local_normal);
    }

};  // struct cylindrical3D

}  // namespace detray
