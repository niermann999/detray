/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"

namespace detray::detail {

template <typename algebra_t>
struct track_helper {

    using scalar_type = dscalar<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using point2_type = dpoint2D<algebra_t>;

    /// Track vector types
    using bound_vector_type = bound_vector<algebra_t>;
    using free_vector_type = free_vector<algebra_t>;

    DETRAY_HOST_DEVICE
    inline point3_type pos(const free_vector_type& free_vec) const {
        return {getter::element(free_vec, e_free_pos0, 0u),
                getter::element(free_vec, e_free_pos1, 0u),
                getter::element(free_vec, e_free_pos2, 0u)};
    }

    DETRAY_HOST_DEVICE
    inline void set_pos(free_vector_type& free_vec, const point3_type& pos) {
        getter::element(free_vec, e_free_pos0, 0u) = pos[0];
        getter::element(free_vec, e_free_pos1, 0u) = pos[1];
        getter::element(free_vec, e_free_pos2, 0u) = pos[2];
    }

    DETRAY_HOST_DEVICE
    inline vector3_type dir(const free_vector_type& free_vec) const {
        return {getter::element(free_vec, e_free_dir0, 0u),
                getter::element(free_vec, e_free_dir1, 0u),
                getter::element(free_vec, e_free_dir2, 0u)};
    }

    DETRAY_HOST_DEVICE
    inline void set_dir(free_vector_type& free_vec, const vector3_type& dir) {
        getter::element(free_vec, e_free_dir0, 0u) = dir[0];
        getter::element(free_vec, e_free_dir1, 0u) = dir[1];
        getter::element(free_vec, e_free_dir2, 0u) = dir[2];
    }

    DETRAY_HOST_DEVICE
    inline point2_type bound_local(const bound_vector_type& bound_vec) const {
        return {getter::element(bound_vec, e_bound_loc0, 0u),
                getter::element(bound_vec, e_bound_loc1, 0u)};
    }

    DETRAY_HOST_DEVICE
    inline vector3_type dir(const bound_vector_type& bound_vec) const {
        const scalar_type phi{getter::element(bound_vec, e_bound_phi, 0u)};
        const scalar_type theta{getter::element(bound_vec, e_bound_theta, 0u)};
        const scalar_type sinTheta{math::sin(theta)};

        return {math::cos(phi) * sinTheta, math::sin(phi) * sinTheta,
                math::cos(theta)};
    }

    DETRAY_HOST_DEVICE
    inline scalar_type p(const free_vector_type& free_vec,
                         const scalar_type q) const {
        assert(qop(free_vec) != 0.f);
        assert(q * qop(free_vec) > 0.f);
        return q / qop(free_vec);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type p(const bound_vector_type& bound_vec,
                         const scalar_type q) const {
        assert(qop(bound_vec) != 0.f);
        assert(q * qop(bound_vec) > 0.f);
        return q / qop(bound_vec);
    }

    DETRAY_HOST_DEVICE
    inline vector3_type mom(const free_vector_type& free_vec,
                            const scalar_type q) const {
        return p(free_vec, q) * dir(free_vec);
    }

    DETRAY_HOST_DEVICE
    inline vector3_type mom(const bound_vector_type& bound_vec,
                            const scalar_type q) const {
        return p(bound_vec, q) * dir(bound_vec);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qop(const free_vector_type& free_vec) const {
        return getter::element(free_vec, e_free_qoverp, 0u);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qop(const bound_vector_type& bound_vec) const {
        return getter::element(bound_vec, e_bound_qoverp, 0u);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qopT(const free_vector_type& free_vec) const {
        const auto dir = this->dir(free_vec);
        assert(vector::perp(dir) != 0.f);
        return getter::element(free_vec, e_free_qoverp, 0u) / vector::perp(dir);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qopT(const bound_vector_type& bound_vec) const {
        const scalar_type theta{getter::element(bound_vec, e_bound_theta, 0u)};
        const scalar_type sinTheta{math::sin(theta)};
        assert(sinTheta != 0.f);
        return getter::element(bound_vec, e_bound_qoverp, 0u) / sinTheta;
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qopz(const free_vector_type& free_vec) const {
        const auto dir = this->dir(free_vec);
        return getter::element(free_vec, e_free_qoverp, 0u) / dir[2];
    }

    DETRAY_HOST_DEVICE
    inline scalar_type qopz(const bound_vector_type& bound_vec) const {
        const scalar_type theta{getter::element(bound_vec, e_bound_theta, 0u)};
        const scalar_type cosTheta{math::cos(theta)};
        assert(cosTheta != 0.f);
        return getter::element(bound_vec, e_bound_qoverp, 0u) / cosTheta;
    }

    DETRAY_HOST_DEVICE
    inline scalar_type time(const free_vector_type& free_vec) const {
        return getter::element(free_vec, e_free_time, 0u);
    }

    DETRAY_HOST_DEVICE
    inline scalar_type time(const bound_vector_type& bound_vec) const {
        return getter::element(bound_vec, e_bound_time, 0u);
    }
};

}  // namespace detray::detail
