/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"

// System include(s)
#include <array>

namespace detray {

/// Base type that defines the interface of all intersectors towards the
/// navigation
template <typename intersector_impl>
struct intersector_base : public intersector_impl {

    /// linear algebra types
    /// @{
    using algebra_type = typename intersector_impl::algebra_type;
    using T = typename algebra_type::value_type;
    using scalar_type = typename intersector_impl::scalar_type;
    using transform3_type = typename intersector_impl::transform3_type;
    /// @}

    template <typename surface_descr_t>
    using intersection_type =
        typename intersector_impl::template intersection_type<surface_descr_t>;

    /// Ray or helix
    template <typename other_algebra_t>
    using trajectory_type =
        typename intersector_impl::template trajectory_type<other_algebra_t>;

    /// Use the default intersector call interface
    using intersector_impl::operator();

    /// Interface to use a fixed mask tolerance @param mask_tolerance
    /*template <typename detector_t, typename other_algebra_t>
    DETRAY_HOST_DEVICE constexpr decltype(auto) operator()(
        const detector_t &det,
        const trajectory_type<other_algebra_t> &traj,
        const typename detector_t::surface_type &sf,
        const transform3_type &trf,
        const scalar_type mask_tolerance,
        const scalar_type overstep_tol = 0.f) const {
        return intersector_impl{}(traj, sf, mask, trf, {mask_tolerance, 0.f},
                                  0.f, overstep_tol);
    }*/

    /// Interface to use a fixed mask tolerance @param mask_tolerance
    template <typename other_algebra_t, typename surface_descr_t,
              typename mask_t>
    DETRAY_HOST_DEVICE constexpr decltype(auto) operator()(
        const trajectory_type<other_algebra_t> &traj, const surface_descr_t &sf,
        const mask_t &mask, const transform3_type &trf,
        const scalar_type mask_tolerance,
        const scalar_type overstep_tol = 0.f) const {
        return intersector_impl{}(traj, sf, mask, trf, {mask_tolerance, 0.f},
                                  0.f, overstep_tol);
    }

    /// Operator function to update an intersection
    ///
    /// @tparam mask_t is the input mask type
    ///
    /// @param traj is the input trajectory
    /// @param sfi the intersection to be updated
    /// @param mask is the input mask that defines the surface extent
    /// @param trf is the surface placement transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tol negative cutoff for the path
    template <typename other_algebra_t, typename surface_descr_t,
              typename mask_t>
    DETRAY_HOST_DEVICE constexpr void update(
        const trajectory_type<other_algebra_t> &traj,
        intersection_type<surface_descr_t> &sfi, const mask_t &mask,
        const transform3_type &trf,
        const std::array<scalar_type, 2u> &mask_tolerance = {0.f,
                                                             1.f * unit<T>::mm},
        const scalar_type mask_tol_scalor = 0.f,
        const scalar_type overstep_tol = 0.f) const {

        using result_t = std::invoke_result_t<
            intersector_impl, trajectory_type<other_algebra_t>, surface_descr_t,
            mask_t, transform3_type, std::array<scalar_type, 2u>, scalar_type,
            scalar_type>;

        if constexpr (std::same_as<
                          std::array<intersection_type<surface_descr_t>, 2>,
                          result_t>) {
            intersector_impl{}.update(traj, sfi, mask, trf, mask_tolerance,
                                      mask_tol_scalor, overstep_tol);
        } else {
            sfi =
                intersector_impl{}(traj, sfi.sf_desc, mask, trf, mask_tolerance,
                                   mask_tol_scalor, overstep_tol);
        }
    }

    /// From the intersection path, construct an intersection candidate and
    /// check it against the surface boundaries (mask).
    ///
    /// @returns the intersection candidate. Might be (partially) uninitialized
    /// if the overstepping tolerance is not met or the intersection lies
    /// outside of the mask.
    /*template <typename surface_descr_t, typename mask_t>
    DETRAY_HOST_DEVICE constexpr void
    mask_check(intersection_type<surface_descr_t>& is,
                    const ray_type &ray, mask_t &mask,
                    const transform3_type &trf,
                    const std::array<scalar_type, 2u> mask_tolerance,
                    const scalar_type mask_tol_scalor) const {

        // Construct the candidate only when needed
        const point3_type p3{ray.pos() + is.path * ray.dir()};

        const auto loc{mask_t::to_local_frame(trf, p3)};
        if constexpr (intersection_type<surface_descr_t>::is_debug()) {
            is.local = loc;
        }
        // Tolerance: per mille of the distance
        is.status = mask.is_inside(
            loc,
            math::max(mask_tolerance[0],
                        math::min(mask_tolerance[1],
                                mask_tol_scalor * math::fabs(is.path))));
        is.direction = !detail::signbit(is.path);
        is.volume_link = mask.volume_link();
    }*/
};

}  // namespace detray
