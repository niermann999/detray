/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/barcode.hpp"

namespace detray {

template <typename algebra_t>
struct bound_track_vector {

    /// @name Type definitions for the struct
    /// @{
    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using point2_type = dpoint2D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using matrix_operator = dmatrix_operator<algebra_t>;
    /// @}

    /// Default constructor
    constexpr bound_track_vector() = default;

    /// Construct from a 6-dim vector of parameters
    DETRAY_HOST_DEVICE
    explicit bound_track_vector(const bound_vector<algebra_t>& vec) {
        m_bound_local[e_bound_loc0] =
            matrix_operator().element(vec, e_bound_loc0, 0u);
        m_bound_local[e_bound_loc1] =
            matrix_operator().element(vec, e_bound_loc1, 0u);
        m_phi = matrix_operator().element(vec, e_bound_phi, 0u);
        m_theta = matrix_operator().element(vec, e_bound_theta, 0u);
        m_qop = matrix_operator().element(vec, e_bound_qoverp, 0u);
        m_time = matrix_operator().element(vec, e_bound_time, 0u);
    }

    /// Construct from single parameters
    ///
    /// @param loc_p the bound local position
    /// @param phi the global phi angle of the track direction
    /// @param theta the global theta angle of the track direction
    /// @param qop the q/p value
    /// @param t the time
    DETRAY_HOST_DEVICE
    constexpr bound_track_vector(const point2_type& loc_p,
                                 const scalar_type phi, const scalar_type theta,
                                 const scalar_type qop, const scalar_type t)
        : m_bound_local{loc_p},
          m_phi{phi},
          m_theta{theta},
          m_qop{qop},
          m_time{t} {}

    /// @param rhs is the left hand side params for comparison
    DETRAY_HOST_DEVICE
    bool operator==(const bound_track_vector& rhs) const {

        constexpr auto e{std::numeric_limits<scalar_type>::epsilon()};

        bool ret{math::fabs(m_bound_local[0] - rhs.m_bound_local[0]) <= e};
        ret &= (math::fabs(m_bound_local[1] - rhs.m_bound_local[1]) <= e);
        ret &= (math::fabs(m_phi - rhs.m_phi) <= e);
        ret &= (math::fabs(m_theta - rhs.m_theta) <= e);
        ret &= (math::fabs(qop() - rhs.qop()) <= e);
        ret &= (math::fabs(time() - rhs.time()) <= e);

        return ret;
    }

    /// Convenience access to the track parameters - const
    DETRAY_HOST_DEVICE
    constexpr scalar_type operator[](const std::size_t i) const {
        switch (static_cast<bound_indices>(i)) {
            case e_bound_loc0:
                return m_bound_local[e_bound_loc0];
            case e_bound_loc1:
                return m_bound_local[e_bound_loc1];
            case e_bound_phi:
                return m_phi;
            case e_bound_theta:
                return m_theta;
            case e_bound_qoverp:
                return m_qop;
            case e_bound_time:
                return m_time;
            default:
                return 0.f;
        }
    }

    /// Convenience access to the track parameters - non-const
    DETRAY_HOST_DEVICE
    constexpr scalar_type& operator[](const std::size_t i) {
        switch (static_cast<bound_indices>(i)) {
            case e_bound_loc0:
                return m_bound_local[e_bound_loc0];
            case e_bound_loc1:
                return m_bound_local[e_bound_loc1];
            case e_bound_phi:
                return m_phi;
            case e_bound_theta:
                return m_theta;
            case e_bound_qoverp:
                return m_qop;
            case e_bound_time:
                return m_time;
        }
    }

    /// Access the track parameters as a 6-dim vector - const
    DETRAY_HOST_DEVICE
    bound_vector<algebra_t> vector() const {
        bound_vector<algebra_t> vector{};

        matrix_operator().set_block(vector, m_bound_local, e_bound_loc0, 0u);
        matrix_operator().element(vector, e_bound_phi, 0u) = m_phi;
        matrix_operator().element(vector, e_bound_theta, 0u) = m_theta;
        matrix_operator().element(vector, e_bound_qoverp, 0u) = m_qop;
        matrix_operator().element(vector, e_bound_time, 0u) = m_time;

        return vector;
    }

    /// @returns the bound local position
    DETRAY_HOST_DEVICE
    constexpr const point2_type& bound_local() const { return m_bound_local; }

    /// Set the bound local position
    DETRAY_HOST_DEVICE
    constexpr void set_bound_local(const point2_type& pos) {
        m_bound_local = pos;
    }

    /// @returns the global phi angle
    DETRAY_HOST_DEVICE
    constexpr scalar_type phi() const { return m_phi; }

    /// Set the global phi angle
    DETRAY_HOST_DEVICE
    constexpr void set_phi(const scalar_type phi) {
        assert(math::abs(phi) <= constant<scalar_type>::pi);
        m_phi = phi;
    }

    /// @returns the global theta angle
    DETRAY_HOST_DEVICE
    constexpr scalar_type theta() const { return m_theta; }

    /// Set the global theta angle
    DETRAY_HOST_DEVICE
    constexpr void set_theta(const scalar_type theta) {
        assert(0.f < theta);
        assert(theta <= constant<scalar_type>::pi);
        m_theta = theta;
    }

    /// @returns the global track direction
    DETRAY_HOST_DEVICE
    constexpr vector3_type dir() const {
        const scalar_type sinTheta{math::sin(m_theta)};

        return {math::cos(m_phi) * sinTheta, math::sin(m_phi) * sinTheta,
                math::cos(m_theta)};
    }

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
        const scalar_type sinTheta{math::sin(m_theta)};
        assert(sinTheta != 0.f);
        return m_qop / sinTheta;
    }

    /// @returns the q/p_z value
    DETRAY_HOST_DEVICE
    constexpr scalar_type qopz() const {
        const scalar_type cosTheta{math::cos(m_theta)};
        assert(cosTheta != 0.f);
        return m_qop / cosTheta;
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
        assert(qop() != 0.f);
        return math::fabs(q / qop() * getter::perp(dir()));
    }

    /// @returns the absolute momentum z-component
    DETRAY_HOST_DEVICE
    constexpr scalar_type pz(const scalar_type q) const {
        assert(qop() != 0.f);
        return math::fabs(q / qop() * dir()[2]);
    }

    protected:
    point2_type m_bound_local{0.f, 0.f};
    scalar_type m_phi{0.f};
    scalar_type m_theta{0.f};
    scalar_type m_qop{0.f};
    scalar_type m_time{0.f};
};

/// Combine the bound track parameter vector with the covariance and associated
/// surface
template <typename algebra_t>
struct bound_track_parameters : public bound_track_vector<algebra_t> {

    using base_type = bound_track_vector<algebra_t>;

    /// @name Type definitions for the struct
    /// @{
    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using matrix_operator = dmatrix_operator<algebra_t>;

    // Shorthand vector/matrix types related to bound track parameters.
    using track_vector_type = base_type;
    using covariance_type = bound_matrix<algebra_t>;

    /// @}

    /// Default constructor sets the covaraicne to zero
    bound_track_parameters() = default;

    /// Construct from the track parameters and set the covariance to zero
    DETRAY_HOST_DEVICE
    explicit bound_track_parameters(const base_type& vec)
        : base_type(vec),
          m_covariance(
              matrix_operator().template zero<e_bound_size, e_bound_size>()),
          m_barcode() {}

    /// Construct with associated surface barcode and covariance
    DETRAY_HOST_DEVICE
    bound_track_parameters(const geometry::barcode bcd, const base_type& vec,
                           const covariance_type& cov)
        : base_type(vec), m_covariance(cov), m_barcode(bcd) {}

    /// Equality operator
    DETRAY_HOST_DEVICE
    bool operator==(const bound_track_parameters& rhs) const {
        if (m_barcode != rhs.surface_link()) {
            return false;
        }

        if (!base_type::operator==(rhs)) {
            return false;
        }

        for (unsigned int i = 0u; i < e_bound_size; i++) {
            for (unsigned int j = 0u; j < e_bound_size; j++) {
                const auto lhs_val =
                    matrix_operator().element(m_covariance, i, j);
                const auto rhs_val =
                    matrix_operator().element(rhs.covariance(), i, j);

                if (math::fabs(lhs_val - rhs_val) >
                    std::numeric_limits<scalar_type>::epsilon()) {
                    return false;
                }
            }
        }
        return true;
    }

    /// Set the track parameter vector
    // DETRAY_HOST_DEVICE
    // void set_vector(const track_vector_type& v) { this->m_vector =
    // v.m_vector; }

    /// @returns the barcode of the associated surface
    DETRAY_HOST_DEVICE
    geometry::barcode surface_link() const { return m_barcode; }

    /// Set the barcode of the associated surface
    DETRAY_HOST_DEVICE
    void set_surface_link(geometry::barcode link) { m_barcode = link; }

    /// @returns the track parameter covariance - non-const
    DETRAY_HOST_DEVICE
    covariance_type& covariance() { return m_covariance; }

    /// @returns the track parameter covariance - const
    DETRAY_HOST_DEVICE
    const covariance_type& covariance() const { return m_covariance; }

    /// Set the track parameter covariance
    DETRAY_HOST_DEVICE
    void set_covariance(const covariance_type& c) { m_covariance = c; }

    private:
    covariance_type m_covariance =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    geometry::barcode m_barcode{};
};

}  // namespace detray
