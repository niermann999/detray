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

namespace detray {

template <typename algebra_t>
struct free_track_parameters {

    /// @name Type definitions for the struct
    /// @{
    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;
    using matrix_operator = dmatrix_operator<algebra_t>;
    /// @}

    /// Default constructor
    constexpr free_track_parameters() = default;

    /// Construct from a 6-dim vector of parameters
    DETRAY_HOST_DEVICE
    explicit free_track_parameters(const free_vector<algebra_t>& vec) {
        m_pos[0] = matrix_operator().element(vec, e_free_pos0, 0u);
        m_pos[1] = matrix_operator().element(vec, e_free_pos1, 0u);
        m_pos[2] = matrix_operator().element(vec, e_free_pos2, 0u);
        m_dir[0] = matrix_operator().element(vec, e_free_dir0, 0u);
        m_dir[1] = matrix_operator().element(vec, e_free_dir1, 0u);
        m_dir[2] = matrix_operator().element(vec, e_free_dir2, 0u);
        m_qop = matrix_operator().element(vec, e_free_qoverp, 0u);
        m_time = matrix_operator().element(vec, e_free_time, 0u);
    }

    /// Construct from single parameters
    ///
    /// @param pos the global position
    /// @param time the time
    /// @param mom the global track momentum 3-vector
    /// @param q the particle charge
    DETRAY_HOST_DEVICE
    constexpr free_track_parameters(const point3_type& pos,
                                    const scalar_type time,
                                    const vector3_type& mom,
                                    const scalar_type q)
        : m_pos{pos},
          m_dir{vector::normalize(mom)},
          m_qop{q / getter::norm(mom)},
          m_time{time} {}

    /// @param rhs is the left hand side params for comparison
    DETRAY_HOST_DEVICE
    constexpr bool operator==(const free_track_parameters& rhs) const {
        for (unsigned int i = 0u; i < e_free_size; i++) {
            if (math::fabs((*this)[i] - rhs[i]) >
                std::numeric_limits<scalar_type>::epsilon()) {
                return false;
            }
        }

        return true;
    }

    /// Convenience access to the track parameters - const
    DETRAY_HOST_DEVICE
    constexpr scalar_type operator[](std::size_t i) const {
        switch (static_cast<free_indices>(i)) {
            case e_free_pos0:
                return m_pos[0];
            case e_free_pos1:
                return m_pos[1];
            case e_free_pos2:
                return m_pos[2];
            case e_free_dir0:
                return m_dir[0];
            case e_free_dir1:
                return m_dir[1];
            case e_free_dir2:
                return m_dir[2];
            case e_free_qoverp:
                return m_qop;
            case e_free_time:
                return m_time;
            default:
                return 0.f;
        }
    }

    /// Convenience access to the track parameters - non-const
    DETRAY_HOST_DEVICE
    constexpr scalar_type& operator[](std::size_t i) {
        switch (static_cast<free_indices>(i)) {
            case e_free_pos0:
                return m_pos[0];
            case e_free_pos1:
                return m_pos[1];
            case e_free_pos2:
                return m_pos[2];
            case e_free_dir0:
                return m_dir[0];
            case e_free_dir1:
                return m_dir[1];
            case e_free_dir2:
                return m_dir[2];
            case e_free_qoverp:
                return m_qop;
            case e_free_time:
                return m_time;
        }
    }

    /// @returns the global track position - const
    DETRAY_HOST_DEVICE
    constexpr const point3_type& pos() const { return m_pos; }

    /// @returns the global track position - non-const
    DETRAY_HOST_DEVICE
    constexpr point3_type& pos() { return m_pos; }

    /// Set the global track position
    DETRAY_HOST_DEVICE
    constexpr void set_pos(const vector3_type& pos) { m_pos = pos; }

    /// @returns the normalized, global track direction - const
    DETRAY_HOST_DEVICE
    constexpr const vector3_type& dir() const { return m_dir; }

    /// @returns the normalized, global track direction - non-const
    DETRAY_HOST_DEVICE
    constexpr vector3_type& dir() { return m_dir; }

    /// Set the global track direction
    /// @note Must be normalized!
    DETRAY_HOST_DEVICE
    constexpr void set_dir(const vector3_type& dir) { m_dir = dir; }

    /// @returns the time
    DETRAY_HOST_DEVICE
    constexpr scalar_type time() const { return m_time; }

    /// Set the time
    DETRAY_HOST_DEVICE
    constexpr void set_time(const scalar_type t) {
        assert(0.f <= t);
        m_time = t;
    }

    /// @returns the q/p value
    DETRAY_HOST_DEVICE
    constexpr scalar_type qop() const { return m_qop; }

    /// Set the q/p value
    DETRAY_HOST_DEVICE
    constexpr void set_qop(const scalar_type qop) { m_qop = qop; }

    /// @returns the q/p_T value
    DETRAY_HOST_DEVICE
    constexpr scalar_type qopT() const {
        const auto& dir = this->dir();
        assert(getter::perp(dir) != 0.f);
        return m_qop / getter::perp(dir);
    }

    /// @returns the q/p_z value
    DETRAY_HOST_DEVICE
    constexpr scalar_type qopz() const {
        const auto& dir = this->dir();
        return m_qop / dir[2];
    }

    /// @returns the absolute momentum
    DETRAY_HOST_DEVICE
    constexpr scalar_type p(const scalar_type q) const {
        assert(qop() != 0.f);
        assert(q * qop() > 0.f);
        return q / qop();
    }

    /// @returns the global momentum 3-vector
    DETRAY_HOST_DEVICE
    constexpr vector3_type mom(const scalar_type q) const {
        return p(q) * dir();
    }

    /// @returns the transverse momentum
    DETRAY_HOST_DEVICE
    constexpr scalar_type pT(const scalar_type q) const {
        assert(this->qop() != 0.f);
        return math::fabs(q / this->qop() * getter::perp(this->dir()));
    }

    /// @returns the absolute momentum z-component
    DETRAY_HOST_DEVICE
    constexpr scalar_type pz(const scalar_type q) const {
        assert(this->qop() != 0.f);
        return math::fabs(q / this->qop() * this->dir()[2]);
    }

    private:
    point3_type m_pos{0.f, 0.f, 0.f};
    vector3_type m_dir{0.f, 0.f, 0.f};
    scalar_type m_qop{0.f};
    scalar_type m_time{0.f};
};

}  // namespace detray
