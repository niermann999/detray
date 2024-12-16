/** Detray library, part of the ACTS project (R&D line)
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
#include "detray/definitions/units.hpp"
#include "detray/geometry/detail/concepts.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/intersection/intersection.hpp"

// System include(s)
#include <type_traits>

namespace detray {

/// A functor to find intersections between straight line and planar surface
template <algebra::concepts::aos algebra_t, bool do_debug>
struct ray_plane_intersector {
    /// linear algebra types
    /// @{
    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;
    /// @}

    template <typename surface_descr_t>
    using intersection_type =
        intersection2D<surface_descr_t, algebra_t, do_debug>;
    using ray_type = detail::ray<algebra_t>;
    template <typename other_algebra_t>
    using trajectory_type = detail::ray<other_algebra_t>;

    /// Operator function to find intersections between ray and planar mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam surface_descr_t is the type of surface handle
    ///
    /// @param ray is the input ray trajectory
    /// @param sf the surface handle the mask is associated with
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tol negative cutoff for the path
    ///
    /// @return the intersection
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE constexpr intersection_type<surface_descr_t> operator()(
        const ray_type &ray, const surface_descr_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const darray<scalar_type, 2u> mask_tolerance =
            {0.f, 1.f * unit<scalar_type>::mm},
        const scalar_type mask_tol_scalor = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        intersection_type<surface_descr_t> is;

        // Retrieve the surface normal & translation (context resolved)
        const vector3_type &sn = trf.z();
        const vector3_type &st = trf.translation();

        // Intersection code
        const point3_type &ro = ray.pos();
        const vector3_type &rd = ray.dir();
        const scalar_type denom = vector::dot(rd, sn);
        // this is dangerous
        if (denom != 0.f) {
            is.path = vector::dot(sn, st - ro) / denom;

            // Intersection is valid for navigation - continue
            if (is.path >= overstep_tol) {

                const point3_type p3 = ro + is.path * rd;
                const auto loc{mask_t::to_local_frame(trf, p3, rd)};
                if constexpr (intersection_type<surface_descr_t>::is_debug()) {
                    is.local = loc;
                }
                // Tolerance: per mille of the distance
                is.status = mask.is_inside(
                    loc, math::max(
                             mask_tolerance[0],
                             math::min(mask_tolerance[1],
                                       mask_tol_scalor * math::fabs(is.path))));
                is.sf_desc = sf;
                is.direction = !detail::signbit(is.path);
                is.volume_link = mask.volume_link();
            }
        } else {
            is.status = false;
        }

        return is;
    }
};

template <typename frame_t, typename algebra_t, bool do_debug>
struct ray_intersector_impl;

template <concepts::planar frame_t, algebra::concepts::aos algebra_t,
          bool do_debug>
struct ray_intersector_impl<frame_t, algebra_t, do_debug> {
    using type = ray_plane_intersector<algebra_t, do_debug>;
};

}  // namespace detray
