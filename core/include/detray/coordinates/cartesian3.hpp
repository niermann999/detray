/** Detray library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/coordinate_base.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/** Frame projection into a cartesian coordinate frame
 */
template <typename transform3_t>
struct cartesian3 final : public coordinate_base<cartesian3, transform3_t> {

    /// @name Type definitions for the struct
    /// @{

    /// Base type
    using base_type = coordinate_base<cartesian3, transform3_t>;
    /// Sclar type
    using scalar_type = typename base_type::scalar_type;
    /// Point in 2D space
    using point2 = typename base_type::point2;
    /// Point in 3D space
    using point3 = typename base_type::point3;
    /// Vector in 3D space
    using vector3 = typename base_type::vector3;

    // Local point type in 3D cartesian coordinates
    using loc_point = point3;

    /// @}

    /// @returns the point @param p in local coordinates
    DETRAY_HOST_DEVICE
    inline point3 global_to_local(const transform3_t &trf, const point3 &p,
                                  const vector3 & /*d*/) const {
        return trf.point_to_local(p);
    }

    /// @returns the local point @param p in global coordinates
    DETRAY_HOST_DEVICE inline point3 local_to_global(const transform3_t &trf,
                                                     const point3 &p) const {
        return trf.point_to_global(p);
    }

    /// @returns the vector @param v in local coordinates
    DETRAY_HOST_DEVICE
    inline vector3 vector_to_local(const transform3_t &trf, const vector3 &v,
                                   const vector3 & = {}) const {
        return trf.vector_to_local(v);
    }

    /// @returns the local vector @param v in global coordinates
    DETRAY_HOST_DEVICE inline vector3 vector_to_global(const transform3_t &trf,
                                                       const vector3 &v) const {
        return trf.vector_to_global(v);
    }

};  // struct cartesian3

}  // namespace detray
