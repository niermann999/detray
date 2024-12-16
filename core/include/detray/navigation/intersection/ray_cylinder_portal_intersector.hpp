/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/concentric_cylindrical2D.hpp"
#include "detray/geometry/coordinates/cylindrical2D.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_cylinder_intersector.hpp"
#include "detray/utils/quadratic_equation.hpp"

// System include(s)
#include <type_traits>

namespace detray {

/// @brief A functor to find intersections between a straight line and a
/// cylindrical portal surface.
///
/// With the way the navigation works, only the closest one of the two possible
/// intersection points is needed in the case of a cylinderical portal surface.
template <algebra::concepts::aos algebra_t, bool do_debug>
struct ray_cylinder_portal_intersector
    : public ray_cylinder_intersector<algebra_t, do_debug> {

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

    /// Operator function to find intersections between ray and cylinder mask
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
    /// @return the closest intersection
    template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE inline intersection_type<surface_descr_t> operator()(
        const ray_type &ray, const surface_descr_t &sf, const mask_t &mask,
        const transform3_type &trf,
        const darray<scalar_type, 2u> mask_tolerance =
            {0.f, 1.f * unit<scalar_type>::mm},
        const scalar_type mask_tol_scalor = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        intersection_type<surface_descr_t> is;

        // Intersecting the cylinder from the inside yield one intersection
        // along the direction of the track and one behind it
        const auto qe = this->solve_intersection(ray, mask, trf);

        // Find the closest valid intersection
        if (qe.solutions() > 0 && qe.larger() > overstep_tol) {
            // Only the closest intersection that is outside the overstepping
            // tolerance is needed
            const scalar_type t{(qe.smaller() > overstep_tol) ? qe.smaller()
                                                              : qe.larger()};
            is = this->template build_candidate<surface_descr_t>(
                ray, mask, trf, t, mask_tolerance, mask_tol_scalor,
                overstep_tol);
            is.sf_desc = sf;
        } else {
            is.status = false;
        }

        return is;
    }
};

template <typename frame_t, typename algebra_t, bool do_debug>
struct ray_intersector_impl;

template <algebra::concepts::aos algebra_t, bool do_debug>
struct ray_intersector_impl<concentric_cylindrical2D<algebra_t>, algebra_t,
                            do_debug> {
    using type = ray_cylinder_portal_intersector<algebra_t, do_debug>;
};

}  // namespace detray
