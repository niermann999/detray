/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"

namespace detray {

template <typename algebra_t>
struct parameter_transporter : actor {

    struct state {};

    /// Mask store visitor
    struct kernel {

        /// @name Type definitions for the struct
        /// @{

        // Transformation matching this struct
        using transform3_type = dtransform3D<algebra_t>;
        // scalar_type
        using scalar_type = dscalar<algebra_t>;

        /// @}

        template <typename mask_group_t, typename index_t,
                  typename propagator_state_t>
        DETRAY_HOST_DEVICE inline void operator()(
            const mask_group_t& /*mask_group*/, const index_t& /*index*/,
            const transform3_type& trf3, propagator_state_t& propagation) {

            using frame_t = typename mask_group_t::value_type::shape::
                template local_frame_type<algebra_t>;

            using jacobian_engine_t = detail::jacobian_engine<frame_t>;

            using bound_matrix_t = bound_matrix<algebra_t>;
            using bound_to_free_matrix_t =
                typename jacobian_engine_t::bound_to_free_matrix_type;

            using free_matrix_t = free_matrix<algebra_t>;
            using free_to_bound_matrix_t =
                typename jacobian_engine_t::free_to_bound_matrix_type;

            // Stepper and Navigator states
            auto& stepping = propagation._stepping;

            // Free vector
            const auto& free_vec = stepping().vector();

            // Convert free to bound vector
            stepping._bound_params.set_vector(
                detail::free_to_bound_vector<frame_t>(trf3, free_vec));

            // Free to bound jacobian at the destination surface
            const free_to_bound_matrix_t free_to_bound_jacobian =
                jacobian_engine_t::free_to_bound_jacobian(trf3, free_vec);

            // Transport jacobian in free coordinate
            const free_matrix_t& free_transport_jacobian =
                stepping._jac_transport;

            // Path correction factor
            free_matrix_t path_correction = jacobian_engine_t::path_correction(
                stepping().pos(), stepping().dir(), stepping.dtds(),
                stepping.dqopds(), trf3);

            const auto correction_term =
                matrix::identity<free_matrix_t>() + path_correction;

            auto new_cov = matrix::zero<bound_matrix_t>();

            if (propagation.param_type() == parameter_type::e_free) {

                const free_to_bound_matrix_t full_jacobian =
                    free_to_bound_jacobian * correction_term *
                    free_transport_jacobian;

                new_cov = full_jacobian * stepping().covariance() *
                          matrix::transpose(full_jacobian);

                propagation.set_param_type(parameter_type::e_bound);

            } else if (propagation.param_type() == parameter_type::e_bound) {
                // Bound to free jacobian at the departure surface
                const bound_to_free_matrix_t& bound_to_free_jacobian =
                    stepping._jac_to_global;

                stepping._full_jacobian =
                    free_to_bound_jacobian * correction_term *
                    free_transport_jacobian * bound_to_free_jacobian;

                new_cov = stepping._full_jacobian *
                          stepping._bound_params.covariance() *
                          matrix::transpose(stepping._full_jacobian);
            }

            // Calculate surface-to-surface covariance transport
            stepping._bound_params.set_covariance(new_cov);
        }
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& /*actor_state*/,
                                       propagator_state_t& propagation) const {
        const auto& navigation = propagation._navigation;

        // Do covariance transport when the track is on surface
        if (!(navigation.is_on_sensitive() ||
              navigation.encountered_sf_material())) {
            return;
        }

        using geo_cxt_t =
            typename propagator_state_t::detector_type::geometry_context;
        const geo_cxt_t ctx{};

        // Surface
        const auto sf = navigation.get_surface();

        sf.template visit_mask<kernel>(sf.transform(ctx), propagation);

        // Set surface link
        propagation._stepping._bound_params.set_surface_link(sf.barcode());
    }
};  // namespace detray

}  // namespace detray
