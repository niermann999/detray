/** Detray library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/coordinates/concentric_cylindrical2.hpp"
#include "detray/geometry/coordinates/cylindrical2.hpp"
#include "detray/propagator/detail/jacobian_kernel.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/tracks/detail/track_helper.hpp"
#include "detray/utils/invalid_values.hpp"

namespace detray::detail {

/// @brief Specialization for 2D cylindrical frames
template <typename algebra_t>
struct jacobian_kernel<cylindrical2<algebra_t>> {

    /// @name Type definitions for the struct
    /// @{
    using coordinate_frame = cylindrical2<algebra_t>;
    // Transform type
    using transform3_type = algebra_t;
    // Sclar type
    using scalar_type = typename algebra_t::scalar_type;
    // Point in 2D space
    using point2 = typename algebra_t::point2;
    // Point in 3D space
    using point3 = typename algebra_t::point3;
    // Vector in 3D space
    using vector3 = typename algebra_t::vector3;
    // Matrix operator
    using matrix_operator = typename algebra_t::matrix_actor;
    // Matrix size type
    using size_type = typename matrix_operator::size_ty;
    // 2D matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    // Rotation Matrix
    using rotation_matrix = matrix_type<3, 3>;
    // Shorthand vector/matrix types related to bound track parameters.
    using bound_vector = matrix_type<e_bound_size, 1>;
    using bound_matrix = matrix_type<e_bound_size, e_bound_size>;
    // Mapping from bound track parameters.
    using bound_to_free_matrix = matrix_type<e_free_size, e_bound_size>;
    // Shorthand vector/matrix types related to free track parameters.
    using free_vector = matrix_type<e_free_size, 1>;
    using free_matrix = matrix_type<e_free_size, e_free_size>;
    // Mapping from free track parameters.
    using free_to_bound_matrix = matrix_type<e_bound_size, e_free_size>;
    using free_to_path_matrix = matrix_type<1, e_free_size>;
    using path_to_free_matrix = matrix_type<e_free_size, 1>;
    // Track helper
    using track_helper = detail::track_helper<matrix_operator>;
    /// @}

    DETRAY_HOST_DEVICE
    static inline rotation_matrix reference_frame(const transform3_type &trf3,
                                                  const point3 &pos,
                                                  const vector3 &dir) {

        rotation_matrix rot = matrix_operator().template zero<3, 3>();

        // y axis of the new frame is the z axis of cylindrical coordinate
        const auto new_yaxis =
            matrix_operator().template block<3, 1>(trf3.matrix(), 0u, 2u);

        // z axis of the new frame is the vector normal to the cylinder surface
        const point3 local = coordinate_frame::global_to_local(trf3, pos, dir);
        const vector3 new_zaxis = coordinate_frame::normal(trf3, local);

        // x axis
        const vector3 new_xaxis = vector::cross(new_yaxis, new_zaxis);

        matrix_operator().element(rot, 0u, 0u) = new_xaxis[0];
        matrix_operator().element(rot, 1u, 0u) = new_xaxis[1];
        matrix_operator().element(rot, 2u, 0u) = new_xaxis[2];
        matrix_operator().template set_block<3, 1>(rot, new_yaxis, 0u, 1u);
        matrix_operator().element(rot, 0u, 2u) = new_zaxis[0];
        matrix_operator().element(rot, 1u, 2u) = new_zaxis[1];
        matrix_operator().element(rot, 2u, 2u) = new_zaxis[2];

        return rot;
    }

    DETRAY_HOST_DEVICE static inline free_to_path_matrix path_derivative(
        const transform3_type &trf3, const point3 &pos, const vector3 &dir,
        const vector3 & /*dtds*/) {

        free_to_path_matrix derivative =
            matrix_operator().template zero<1u, e_free_size>();

        const point3 local = coordinate_frame::global_to_local(trf3, pos, dir);
        const vector3 normal = coordinate_frame::normal(trf3, local);

        const vector3 pos_term = -1.f / vector::dot(normal, dir) * normal;

        matrix_operator().element(derivative, 0u, e_free_pos0) = pos_term[0];
        matrix_operator().element(derivative, 0u, e_free_pos1) = pos_term[1];
        matrix_operator().element(derivative, 0u, e_free_pos2) = pos_term[2];

        return derivative;
    }

    DETRAY_HOST_DEVICE
    static inline void set_bound_pos_to_free_pos_derivative(
        bound_to_free_matrix &bound_to_free_jacobian,
        const transform3_type &trf3, const point3 &pos, const vector3 &dir) {

        const auto frame = reference_frame(trf3, pos, dir);

        // Get d(x,y,z)/d(loc0, loc1)
        const auto bound_pos_to_free_pos_derivative =
            matrix_operator().template block<3, 2>(frame, 0u, 0u);

        matrix_operator().template set_block(bound_to_free_jacobian,
                                             bound_pos_to_free_pos_derivative,
                                             e_free_pos0, e_bound_loc0);
    }

    DETRAY_HOST_DEVICE
    static inline void set_free_pos_to_bound_pos_derivative(
        free_to_bound_matrix &free_to_bound_jacobian,
        const transform3_type &trf3, const point3 &pos, const vector3 &dir) {

        const auto frame = reference_frame(trf3, pos, dir);
        const auto frameT = matrix_operator().transpose(frame);

        // Get d(loc0, loc1)/d(x,y,z)
        const auto free_pos_to_bound_pos_derivative =
            matrix_operator().template block<2, 3>(frameT, 0u, 0u);

        matrix_operator().template set_block(free_to_bound_jacobian,
                                             free_pos_to_bound_pos_derivative,
                                             e_bound_loc0, e_free_pos0);
    }

    DETRAY_HOST_DEVICE
    static inline void set_bound_angle_to_free_pos_derivative(
        bound_to_free_matrix & /*bound_to_free_jacobian*/,
        const transform3_type & /*trf3*/, const point3 & /*pos*/,
        const vector3 & /*dir*/) {
        // Do nothing
    }
};

template <typename algebra_t>
struct jacobian_kernel<concentric_cylindrical2<algebra_t>>
    : public jacobian_kernel<cylindrical2<algebra_t>> {};

}  // namespace detray::detail
