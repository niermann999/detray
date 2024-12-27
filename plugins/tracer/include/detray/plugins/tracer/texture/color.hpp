/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/concepts.hpp"

// System include(s)
#include <array>
#include <cstdint>
#include <limits>
#include <ostream>
#include <type_traits>

namespace detray::texture {

namespace detail {

/// Determine the value of the maximal color intensity for each color channel
/// type
/// @{
template <typename T>
struct channel_max
    : public std::integral_constant<T, std::numeric_limits<T>::max()> {};

template <>
struct channel_max<float> {
    static constexpr float value{1.f};
};

template <>
struct channel_max<double> {
    static constexpr double value{1.};
};

template <typename T>
inline constexpr auto channel_max_v = channel_max<T>::value;
/// @}

}  // namespace detail

/// @brief holds rgb and alpha values
///
/// @tparam depth_t how to encode the color data (e.g. # bits, AoS vs SoA)
template <concepts::arithmetic depth_t = std::uint8_t>
struct color {

    using T = std::conditional_t<std::floating_point<depth_t>, depth_t, float>;
    using depth = depth_t;

    /// Maximal intensity per color channel
    static constexpr auto channel_max{detail::channel_max_v<depth>};

    /// Default constructor
    consteval color() = default;

    /// Construct from colors @param r (red), @param g (green), @param b (blue)
    /// and @param alpha values
    DETRAY_HOST_DEVICE
    constexpr color(const depth r, const depth g, const depth b,
                    const depth alpha = channel_max)
        : m_data{r, g, b, alpha} {}

    /// Broadcast constructor
    DETRAY_HOST_DEVICE
    explicit constexpr color(const depth_t value)
        : color(value, value, value) {}

    /// Conversion to a color with a different color depth
    template <typename other_depth_t>
    requires std::is_convertible_v<depth_t, other_depth_t>
        DETRAY_HOST_DEVICE constexpr operator color<other_depth_t>() const {
        // Conversion factor to scale maximal color intensity
        constexpr T s{static_cast<T>(detail::channel_max_v<other_depth_t>) /
                      static_cast<T>(channel_max)};

        std::array<T, 4> tmp{
            static_cast<T>(m_data[0]), static_cast<T>(m_data[1]),
            static_cast<T>(m_data[2]), static_cast<T>(m_data[3])};

        for (std::size_t i = 0u; i < 4u; ++i) {
            tmp[i] *= s;
        }

        return color<other_depth_t>{static_cast<other_depth_t>(tmp[0]),
                                    static_cast<other_depth_t>(tmp[1]),
                                    static_cast<other_depth_t>(tmp[2]),
                                    static_cast<other_depth_t>(tmp[3])};
    }

    /// Mix two colors with alpha values: This color over @param other
    DETRAY_HOST_DEVICE
    constexpr color alpha_over(const color& right) {
        const auto r{static_cast<color<T>>(right)};
        const auto l{static_cast<color<T>>(*this)};

        constexpr auto n{1.f / color<T>::channel_max};

        const auto alpha_fraction{r[3] * (1.f - l[3] * n)};
        const auto alpha{l[3] + alpha_fraction};
        const auto inv_alpha{1.f / alpha};

        color<T> res;
        std::array<T, 4> tmp1;
        std::array<T, 4> tmp2;
        for (std::size_t i = 0u; i < 4u; ++i) {
            tmp1[i] = l[i] * l[3];
            tmp2[i] = r[i] * alpha_fraction;
            res[i] = (tmp1[i] + tmp2[i]) * inv_alpha;
        }
        res[3] = alpha;

        return static_cast<color<depth_t>>(res);
    }

    /// Equality operator: Only considers exact match
    DETRAY_HOST_DEVICE
    constexpr depth operator==(const color& other) {
        return m_data == other.m_data;
    }

    /// Subscript operator @returns a color data point - const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) operator[](const std::size_t i) const {
        return m_data[i];
    }

    /// Subscript operator @returns a color data point - non-const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) operator[](const std::size_t i) {
        return m_data[i];
    }

    /// Scale the color by a factor @param factor
    DETRAY_HOST_DEVICE constexpr color& operator*=(const T factor) {
        *this = *this * factor;
        return *this;
    }

    /// Multiply with color @param left
    DETRAY_HOST_DEVICE constexpr color& operator*=(const color& right) {
        *this = *this * right;
        return *this;
    }

    /// Add a color @param left
    DETRAY_HOST_DEVICE constexpr color& operator+=(const color& right) {
        *this = *this + right;
        return *this;
    }

    /// Mixes two colors @param left and @param right by addition
    DETRAY_HOST_DEVICE
    template <typename>
    friend constexpr color operator+(const color& left, const color& right);

    /// Scale the color by a factor @param factor
    DETRAY_HOST_DEVICE
    template <typename, typename scalar_t>
    friend constexpr color operator*(const scalar_t factor, const color& right);

    /// Print the color data to stdout
    DETRAY_HOST
    template <typename>
    friend std::ostream& operator<<(std::ostream& os, const color& c);

    std::array<depth_t, 4> m_data{};
};

template <concepts::arithmetic depth_t>
std::ostream& operator<<(std::ostream& os, const color<depth_t>& c) {
    return os << "rgba: (" << c[0] << ", " << c[1] << ", " << c[2] << ", "
              << c[3] << ")";
}

template <>
std::ostream& operator<<(std::ostream& os, const color<std::uint8_t>& c) {
    return os << "rgba: (" << static_cast<std::size_t>(c[0]) << ", "
              << static_cast<std::size_t>(c[1]) << ", "
              << static_cast<std::size_t>(c[2]) << ", "
              << static_cast<std::size_t>(c[3]) << ")";
}

template <concepts::arithmetic depth_t>
constexpr color<depth_t> operator+(const color<depth_t>& left,
                                   const color<depth_t>& right) {
    color<depth_t> res;
    for (std::size_t i = 0u; i < 4u; ++i) {
        res[i] = left[i] + right[i];
    }

    return res;
}

template <concepts::arithmetic depth_t>
constexpr color<depth_t> operator*(const color<depth_t>& col,
                                   const typename color<depth_t>::T& factor) {
    using T = typename color<depth_t>::T;
    const auto c{static_cast<color<T>>(col)};

    color<T> res;
    for (std::size_t i = 0u; i < 4u; ++i) {
        res[i] = c[i] * factor;
    }

    return static_cast<color<depth_t>>(res);
}

template <concepts::arithmetic depth_t>
constexpr color<depth_t> operator*(const color<depth_t>& left,
                                   const color<depth_t>& right) {
    using T = typename color<depth_t>::T;
    const auto l{static_cast<color<T>>(left)};
    const auto r{static_cast<color<T>>(right)};

    constexpr auto n{1.f / color<T>::channel_max};

    color<T> res;
    for (std::size_t i = 0u; i < 4u; ++i) {
        res[i] = l[i] * r[i] * n;
    }

    return static_cast<color<depth_t>>(res);
}

}  // namespace detray::texture
