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
template <template <typename> class impl, typename algebra_t>
struct local_frame {

    /// @name Linear algebra types of the local frame implementation
    /// @{
    using transform3_type = algebra_t;
    using scalar_type = typename algebra_t::scalar_type;
    using point2 = typename algebra_t::point2;
    using point3 = typename algebra_t::point3;
    using vector3 = typename algebra_t::vector3;
    /// @}

    /// @brief Tranform from the global frame to the local frame
    ///
    /// @param trf3 the transformation to the local frame
    /// @param p the point to be transformed
    /// @param d additional direction needed by some local frames
    ///
    /// @returns a point in the local frame
    DETRAY_HOST_DEVICE
    static constexpr decltype(auto) global_to_local(const transform3_type &trf,
                                                    const point3 &p,
                                                    const vector3 &d) {

        return impl<algebra_t>::project(trf, trf.point_to_local(p), d);
    }

    /// @brief Tranform from the local frame to the global frame
    ///
    /// @param trf3 the transformation to the local frame
    /// @param v the vector to be transformed (no translation)
    /// @param d additional direction needed by some local frames
    ///
    /// @returns a vector in the local frame
    template <typename mask_t, typename loc_point_t>
    DETRAY_HOST_DEVICE static constexpr point3 local_to_global(
        const transform3_type &trf, const mask_t &mask,
        const loc_point_t &loc_p, const vector3 &d) {

        point3 loc_cartesian3D = impl<algebra_t>::project(trf, mask, loc_p, d);

        return trf.point_to_global(loc_cartesian3D);
    }
};

}  // namespace detray
