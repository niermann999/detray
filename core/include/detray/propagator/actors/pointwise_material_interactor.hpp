/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/materials/interaction.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/utils/ranges.hpp"

namespace detray {

template <typename T, template <typename> class algebra_t>
struct pointwise_material_interactor : actor {

    /// @name Linear algebra definitions
    /// @{
    using scalar_type = dscalar<algebra_t<T>>;
    using vector3 = dvector3D<algebra_t<T>>;
    using transform3_type = dtransform3D<algebra_t<T>>;
    using matrix_operator = typename transform3_type::matrix_actor;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    using bound_vector = matrix_type<e_bound_size, 1u>;
    using bound_matrix = matrix_type<e_bound_size, e_bound_size>;
    /// @}

    using interaction_type = interaction<scalar_type>;

    struct state {
        using vector3 = __plugin::vector3<scalar>;

        /// The particle mass
        scalar_type mass{105.7f * unit<scalar_type>::MeV};
        /// The particle pdg
        int pdg = 13;  // default muon

        /// Evaluated energy loss
        scalar_type e_loss{0.f};
        /// Evaluated projected scattering angle
        scalar_type projected_scattering_angle{0.f};
        /// Evaluated sigma of qoverp
        scalar_type sigma_qop{0.f};

        bool do_covariance_transport = true;
        bool do_energy_loss = true;
        bool do_multiple_scattering = true;

        DETRAY_HOST_DEVICE
        void reset() {
            e_loss = 0.f;
            projected_scattering_angle = 0.f;
            sigma_qop = 0.f;
        }
    };

    /// Material store visitor
    struct kernel {

        using scalar_type = typename interaction_type::scalar_type;
        using state = typename pointwise_material_interactor::state;

        template <typename material_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline bool operator()(
            const material_group_t &material_group,
            const index_t &material_range, state &s,
            const bound_track_parameters<transform3_type> &bound_params,
            const scalar_type cos_inc_angle, const scalar_type approach) const {

            const scalar qop = bound_params.qop();
            const scalar charge = bound_params.charge();

            bool success = false;

            for (const auto &mat :
                 detray::ranges::subrange(material_group, material_range)) {

                // return early in case of vacuum or zero thickness
                if (not mat) {
                    continue;
                }
                success |= true;

                const scalar_type path_segment{
                    mat.path_segment(cos_inc_angle, approach)};

                // Energy Loss
                if (s.do_energy_loss) {
                    s.e_loss = interaction_type().compute_energy_loss_bethe(
                        path_segment, mat, s.pdg, s.mass, qop, charge);
                }

                // @todo: include the radiative loss (Bremsstrahlung)
                if (s.do_energy_loss && s.do_covariance_transport) {
                    s.sigma_qop =
                        interaction_type()
                            .compute_energy_loss_landau_sigma_QOverP(
                                path_segment, mat, s.pdg, s.mass, qop, charge);
                }

                // Covariance update
                if (s.do_multiple_scattering) {
                    // @todo: use momentum before or after energy loss in
                    // backward mode?
                    s.projected_scattering_angle =
                        interaction_type().compute_multiple_scattering_theta0(
                            mat.path_segment_in_X0(cos_inc_angle, approach),
                            s.pdg, s.mass, qop, charge);
                }
            }

            return success;
        }
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void operator()(
        state &interactor_state, propagator_state_t &prop_state) const {

        interactor_state.reset();

        const auto &navigation = prop_state._navigation;

        // Do material interaction when the track is on surface
        if (navigation.is_on_module()) {

            auto &stepping = prop_state._stepping;

            this->update(stepping._bound_params, interactor_state,
                         static_cast<int>(navigation.direction()),
                         navigation.get_surface(),
                         navigation.current()->cos_incidence_angle);
        }
    }

    /// @brief Update the bound track parameter
    ///
    /// @param[out] bound_params bound track parameter
    /// @param[out] interactor_state actor state
    /// @param[in]  nav_dir navigation direction
    /// @param[in]  sf the surface
    template <typename surface_t>
    DETRAY_HOST_DEVICE inline void update(
        bound_track_parameters<transform3_type> &bound_params,
        state &interactor_state, const int nav_dir, const surface_t &sf,
        const scalar_type cos_inc_angle) const {

        // Closest approach of the track to a line surface. Otherwise this is
        // ignored.
        const auto approach{
            matrix_operator().element(bound_params.vector(), e_bound_loc0, 0)};

        const bool succeed = sf.template visit_material<kernel>(
            interactor_state, bound_params, cos_inc_angle, approach);

        if (succeed) {

            auto &covariance = bound_params.covariance();
            auto &vector = bound_params.vector();

            if (interactor_state.do_energy_loss) {

                update_qop(vector, bound_params.p(), bound_params.charge(),
                           interactor_state.mass, interactor_state.e_loss,
                           nav_dir);

                if (interactor_state.do_covariance_transport) {

                    update_qop_variance(covariance, interactor_state.sigma_qop,
                                        nav_dir);
                }
            }

            if (interactor_state.do_covariance_transport) {

                update_angle_variance(
                    covariance, bound_params.dir(),
                    interactor_state.projected_scattering_angle, nav_dir);
            }
        }
    }

    /// @brief Update the q over p of bound track parameter
    ///
    /// @param[out] vector vector of bound track parameter
    /// @param[in]  p momentum of the track
    /// @param[in]  q charge of the track
    /// @param[in]  m mass of the track
    /// @param[in]  e_loss energy loss
    /// @param[in]  sign navigation direction
    DETRAY_HOST_DEVICE
    inline void update_qop(bound_vector &vector, const scalar_type p,
                           const scalar_type q, const scalar_type m,
                           const scalar_type e_loss, const int sign) const {

        // Get new Energy
        const scalar_type nextE{
            std::sqrt(m * m + p * p) -
            std::copysign(e_loss, static_cast<scalar_type>(sign))};

        // Put particle at rest if energy loss is too large
        const scalar_type nextP{(m < nextE) ? std::sqrt(nextE * nextE - m * m)
                                            : 0.f};

        // For neutral particles, qoverp = 1/p
        getter::element(vector, e_bound_qoverp, 0) =
            (q != 0.f) ? q / nextP : 1.f / nextP;
    }

    /// @brief Update the variance of q over p of bound track parameter
    ///
    /// @param[out] covariance covariance matrix of bound track parameter
    /// @param[in]  sigma_qop variance of q over p
    /// @param[in]  sign navigation direction
    DETRAY_HOST_DEVICE inline void update_qop_variance(
        bound_matrix &covariance, const scalar_type sigma_qop,
        const int sign) const {

        const scalar_type variance_qop{sigma_qop * sigma_qop};

        matrix_operator().element(covariance, e_bound_qoverp, e_bound_qoverp) +=
            std::copysign(variance_qop, static_cast<scalar_type>(sign));
    }

    /// @brief Update the variance of phi and theta of bound track parameter
    ///
    /// @param[out] covariance covariance matrix of bound track parameter
    /// @param[in]  dir direction of track
    /// @param[in]  projected_scattering_angle projected scattering angle
    /// @param[in]  sign navigation direction
    DETRAY_HOST_DEVICE inline void update_angle_variance(
        bound_matrix &covariance, const vector3 &dir,
        const scalar_type projected_scattering_angle, const int sign) const {

        // variance of projected scattering angle
        const scalar_type var_scattering_angle{std::copysign(
            projected_scattering_angle * projected_scattering_angle,
            static_cast<scalar_type>(sign))};

        matrix_operator().element(covariance, e_bound_phi, e_bound_phi) +=
            var_scattering_angle / (1 - dir[2] * dir[2]);

        matrix_operator().element(covariance, e_bound_theta, e_bound_theta) +=
            var_scattering_angle;
    }
};

}  // namespace detray
