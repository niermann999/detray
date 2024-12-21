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

/// A functor to find intersections between trajectory and line mask
template <algebra::concepts::aos algebra_t, bool do_debug>
struct ray_line_intersector {
    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;

    template <typename surface_descr_t>
    using intersection_type =
        intersection2D<surface_descr_t, algebra_t, do_debug>;
    using ray_type = detail::ray<algebra_t>;
    template <typename other_algebra_t>
    using trajectory_type = detail::ray<other_algebra_t>;

    /// Operator function to find intersections between ray and line mask
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
    //
    /// @return the intersection
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const ray_type &ray, const surface_descr_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const darray<scalar_type, 2u> mask_tolerance =
            {0.f, 1.f * unit<scalar_type>::mm},
        const scalar_type mask_tol_scalor = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        intersection_type<surface_descr_t> is;

        // line direction
        const vector3_type &_z = trf.z();

        // line center
        const point3_type &_t = trf.translation();

        // track direction
        const vector3_type &_d = ray.dir();

        // track position
        const point3_type &_p = ray.pos();

        // Projection of line to track direction
        const scalar_type zd{vector::dot(_z, _d)};

        const scalar_type denom{1.f - (zd * zd)};

        // Case for wire is parallel to track
        if (denom < 1e-5f) {
            is.status = false;
            return is;
        }

        // vector from track position to line center
        const auto t2l = _t - _p;

        // t2l projection on line direction
        const scalar_type t2l_on_line{vector::dot(t2l, _z)};

        // t2l projection on track direction
        const scalar_type t2l_on_track{vector::dot(t2l, _d)};

        // path length to the point of closest approach on the track
        const scalar_type A{1.f / denom * (t2l_on_track - t2l_on_line * zd)};

        is.path = A;
        // Intersection is not valid for navigation - return early
        if (is.path >= overstep_tol) {

            // point of closest approach on the track
            const point3_type m = _p + _d * A;

            const auto loc{mask_t::to_local_frame(trf, m, _d)};
            if constexpr (intersection_type<surface_descr_t>::is_debug()) {
                is.local = loc;
            }
            // Tolerance: per mille of the distance
            is.status = mask.is_inside(
                loc,
                math::max(mask_tolerance[0],
                          math::min(mask_tolerance[1],
                                    mask_tol_scalor * math::fabs(is.path))));
            is.sf_desc = sf;
            is.direction = !detail::signbit(is.path);
            is.volume_link = mask.volume_link();
        }
        return is;
    }
};

template <typename frame_t, typename algebra_t, bool do_debug>
struct ray_intersector_impl;

template <concepts::linear shape_t, algebra::concepts::aos algebra_t,
          bool do_debug>
struct ray_intersector_impl<shape_t, algebra_t, do_debug> {
    using type = ray_line_intersector<algebra_t, do_debug>;
};

}  // namespace detray
