/** Detray library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray {

/// Tranformation to a local coordinate frame
template <class impl>
struct local_frame {

    /// @name Linear algebra types of the local frame implementation
    /// @{

    // Transform type
    using transform3_type = typename impl::transform3;
    // Scalar type
    using scalar_type = typename impl::scalar_type;
    // Point in 2D space
    using point2 = typename impl::point2;
    // Point in 3D space
    using point3 = typename impl::point3;
    // Vector in 3D space
    using vector3 = typename impl::vector3;
    // Local point type
    using local_point = typename impl::local_point;

    /// @}

    /// @brief Tranform from the global frame to the local frame
    ///
    /// @param trf3 the transformation to the local frame
    /// @param p the point to be transformed
    /// @param d additional direction needed by some local frames
    ///
    /// @returns a point in the local frame
    DETRAY_HOST_DEVICE
    static constexpr local_point global_to_local(const transform3 &trf,
                                            const point3 &p,
                                            const vector3 &d) {
        point3 loc_cartesian3D = trf.point_to_local(p);

        return impl::project(trf, loc_cartesian3D, d);
    }

    /// @brief Tranform from the local frame to the global frame
    ///
    /// @param trf3 the transformation to the local frame
    /// @param v the vector to be transformed (no translation)
    /// @param d additional direction needed by some local frames
    ///
    /// @returns a vector in the local frame
    DETRAY_HOST_DEVICE
    static constexpr point3 local_to_global(
        const transform3_t &trf, const mask_t &mask, const loc_point &loc_p,
        const vector3 &d) {

        point3 loc_cartesian3D = impl::project(trf, mask, loc_p, d);

        return trf.point_to_global(loc_cartesian3D);
    }
};

}  // namespace detray
