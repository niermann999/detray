/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/tracks/tracks.hpp"

namespace detray::detail {

/// Transform free track parameters to bound track parameter vector (no
/// covariance needed)
///
/// @param trf3 transform of the surface the bound vector should be defined on
/// @param free_param the free track parameters to be transformed
///
/// @returns the bound track parameter vector corresponding to the free track
///          parameters
template <typename local_frame_t>
DETRAY_HOST_DEVICE inline bound_track_vector<
    typename local_frame_t::algebra_type>
free_to_bound_vector(
    const dtransform3D<typename local_frame_t::algebra_type>& trf3,
    const free_track_parameters<typename local_frame_t::algebra_type>&
        free_param) {

    const auto dir = free_param.dir();

    const auto bound_local =
        local_frame_t::global_to_local(trf3, free_param.pos(), dir);

    // Construct the bound track parametrization without covariance
    return {bound_local, getter::phi(dir), getter::theta(dir), free_param.qop(),
            free_param.time()};
}

/// Convenience function to transform free track parameters to
/// bound track parameters
///
/// @param trf3 transform of the surface the bound vector should be defined on
/// @param free_param the free track parameters to be transformed
/// @param bcd barcode of the surface the bound vector should be defined on
/// @param cov the covariance matrix corresponding to the bound track vector
///
/// @returns the bound track parameters corresponding to the free track
///          parameters
template <typename local_frame_t>
DETRAY_HOST_DEVICE inline bound_track_parameters<
    typename local_frame_t::algebra_type>
free_to_bound_param(
    const dtransform3D<typename local_frame_t::algebra_type>& trf3,
    const free_track_parameters<typename local_frame_t::algebra_type>&
        free_param,
    const geometry::barcode bcd = {},
    const typename bound_track_parameters<
        typename local_frame_t::algebra_type>::covariance_type& cov = {}) {

    // Construct the bound track parametrization with covariance and barcode
    return {bcd, free_to_bound_vector<local_frame_t>(trf3, free_param), cov};
}

/// Transform a bound track parameter vector to free track parameters
///
/// @param trf3 transform of the surface the bound parameters are defined on
/// @param mask the mask of the surface the bound parameters are defined on
/// @param bound_vec the bound track vector to be transformed
///
/// @returns the free track parameters for the bound track parameter vector
template <typename mask_t>
DETRAY_HOST_DEVICE inline free_track_parameters<typename mask_t::algebra_type>
bound_to_free_param(
    const dtransform3D<typename mask_t::algebra_type>& trf3, const mask_t& mask,
    const bound_track_vector<typename mask_t::algebra_type>& bound_vec) {

    using local_frame_t = typename mask_t::local_frame_type;

    const auto dir = bound_vec.dir();

    const auto pos = local_frame_t::local_to_global(
        trf3, mask, bound_vec.bound_local(), dir);

    // Construct the free track parametrization
    return {pos, bound_vec.time(), dir, bound_vec.qop()};
}

}  // namespace detray::detail
