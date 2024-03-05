/** Detray library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/coordinate_base.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/coordinates/cartesian2D.hpp"
#include "detray/tracks/bound_track_parameters.hpp"

namespace detray::detail {

/** Frame projection into a cartesian coordinate frame
 */
template <typename transform3_t>
struct coordinate_base<cartesian2D<transform3_t>> {

    /// @name Type definitions for the struct
    /// @{

    // Transform type
    using transform3_type = transform3_t;
    // The underlying geometric local frame
    using frame_type = cartesian2D<transform3_t>;
    // Sclar type
    using scalar_type = typename transform3_type::scalar_type;
    // Point in 2D space
    using point2 = typename transform3_type::point2;
    // Point in 3D space
    using point3 = typename transform3_type::point3;
    // Vector in 3D space
    using vector3 = typename transform3_type::vector3;
    // Matrix operator
    using matrix_operator = typename transform3_type::matrix_actor;
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
    // Matrix types
    using free_to_bound_matrix = matrix_type<e_bound_size, e_free_size>;
    using bound_to_free_matrix = matrix_type<e_free_size, e_bound_size>;
    using free_to_path_matrix = matrix_type<1, e_free_size>;

    /// @}

    /** This method transforms a point from a global cartesian 3D frame to a
     * local 3D cartesian point */
    DETRAY_HOST_DEVICE
    inline point3 global_to_local(const transform3_type &trf3, const point3 &p,
                                  const vector3 & /*d*/) const {
        return trf3.point_to_local(p);
    }

    /** This method transform from a local 2D cartesian point to a point global
     * cartesian 3D frame*/
    DETRAY_HOST_DEVICE inline point3 local_to_global(
        const transform3_type &trf3, const point3 &p) const {
        return trf3.point_to_global(p);
    }

    /** This method transform from a local 2D cartesian point to a point global
     * cartesian 3D frame*/
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline point3 bound_local_to_global(
        const transform3_type &trf3, const mask_t & /*mask*/, const point2 &p,
        const vector3 & /*d*/) const {

        return this->local_to_global(trf3, {p[0], p[1], 0.f});
    }

    /// @returns the normal vector
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline vector3 normal(const transform3_type &trf3,
                                             const point2 &p = {},
                                             const mask_t &m = {}) const {
        return frame_type::normal(trf3, p, m);
    }

    /// @returns the normal vector
    DETRAY_HOST_DEVICE inline vector3 normal(const transform3_type &trf3,
                                             const point3 & = {}) const {
        return trf3.z();
    }

    DETRAY_HOST_DEVICE inline rotation_matrix reference_frame(
        const transform3_type &trf3, const point3 & /*pos*/,
        const vector3 & /*dir*/) const {
        return trf3.rotation();
    }

    DETRAY_HOST_DEVICE inline free_to_path_matrix path_derivative(
        const transform3_type &trf3, const point3 & /*pos*/, const vector3 &dir,
        const vector3 & /*dtds*/) const {

        free_to_path_matrix derivative =
            matrix_operator().template zero<1u, e_free_size>();

        const vector3 normal = this->normal(trf3);

        const vector3 pos_term = -1.f / vector::dot(normal, dir) * normal;

        matrix_operator().element(derivative, 0u, e_free_pos0) = pos_term[0];
        matrix_operator().element(derivative, 0u, e_free_pos1) = pos_term[1];
        matrix_operator().element(derivative, 0u, e_free_pos2) = pos_term[2];

        return derivative;
    }

    DETRAY_HOST_DEVICE inline void set_bound_pos_to_free_pos_derivative(
        bound_to_free_matrix &bound_to_free_jacobian,
        const transform3_type &trf3, const point3 &pos,
        const vector3 &dir) const {

        const rotation_matrix frame = reference_frame(trf3, pos, dir);

        // Get d(x,y,z)/d(loc0, loc1)
        const matrix_type<3, 2> bound_pos_to_free_pos_derivative =
            matrix_operator().template block<3, 2>(frame, 0u, 0u);

        matrix_operator().template set_block(bound_to_free_jacobian,
                                             bound_pos_to_free_pos_derivative,
                                             e_free_pos0, e_bound_loc0);
    }

    DETRAY_HOST_DEVICE inline void set_free_pos_to_bound_pos_derivative(
        free_to_bound_matrix &free_to_bound_jacobian,
        const transform3_type &trf3, const point3 &pos,
        const vector3 &dir) const {

        const rotation_matrix frame = reference_frame(trf3, pos, dir);
        const rotation_matrix frameT = matrix_operator().transpose(frame);

        // Get d(loc0, loc1)/d(x,y,z)
        const matrix_type<2, 3> free_pos_to_bound_pos_derivative =
            matrix_operator().template block<2, 3>(frameT, 0, 0);

        matrix_operator().template set_block(free_to_bound_jacobian,
                                             free_pos_to_bound_pos_derivative,
                                             e_bound_loc0, e_free_pos0);
    }

    DETRAY_HOST_DEVICE inline void set_bound_angle_to_free_pos_derivative(
        bound_to_free_matrix & /*bound_to_free_jacobian*/,
        const transform3_type & /*trf3*/, const point3 & /*pos*/,
        const vector3 & /*dir*/) const {
        // Do nothing
    }
};  // struct cartesian2

}  // namespace detray::detail
