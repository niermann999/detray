/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <array>
#include <cstdint>
#include <limits>
#include <ostream>
#include <type_traits>

namespace detray::texture {

namespace detail {

/// Determine the value of the maximal color intensity per color depth
template <typename T>
struct max_intensity
    : public std::integral_constant<T, std::numeric_limits<T>::max()> {};

template <>
struct max_intensity<float> {
    static constexpr float value{1.f};
};

template <>
struct max_intensity<double> {
    static constexpr double value{1.};
};

}  // namespace detail

/// @brief holds rgb and alpha values for color shading
///
/// @tparam data_t how to store the rgb data: single color or soa
template <typename data_t = std::uint8_t>
struct color {

    using color_depth = data_t;

    /// Maximal intensity per color channel
    static constexpr auto max_I{detail::max_intensity<data_t>::value};

    /// Default constructor
    constexpr color() = default;

    /// Construct from colors @param r (red), @param g (green), @param b (blue)
    /// and @param alpha values
    DETRAY_HOST_DEVICE
    constexpr color(const data_t r, const data_t g, const data_t b,
                    const data_t alpha = {max_I - 1u})
        : m_data{r, g, b, alpha} {}

    /// Broadcast constructor
    explicit constexpr color(const data_t value) : color(value, value, value) {}

    /// Conversion to surface interface around constant detector type
    template <typename other_data_t,
              std::enable_if_t<std::is_convertible_v<data_t, other_data_t>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr operator color<other_data_t>() const {
        // Conversion factor to scale maximal color intensity
        constexpr float s{
            static_cast<float>(detail::max_intensity<other_data_t>::value) /
            max_I};

        return color<other_data_t>{
            static_cast<other_data_t>(this->m_data[0] * s),
            static_cast<other_data_t>(this->m_data[1] * s),
            static_cast<other_data_t>(this->m_data[2] * s),
            static_cast<other_data_t>(this->m_data[3] * s)};
    }

    /// Equality operator: Only considers exact match
    DETRAY_HOST_DEVICE
    constexpr data_t operator==(const color& other) {
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

    /// Scale the color by a normalization factor @param scalor
    DETRAY_HOST_DEVICE
    constexpr color& operator+=(const color& left) {
        m_data[0] += left.m_data[0];
        m_data[1] += left.m_data[1];
        m_data[2] += left.m_data[2];
        m_data[3] += left.m_data[3];

        return *this;
    }

    /// Scale the color by a normalization factor @param scalor
    template <typename scalar_t,
              std::enable_if_t<std::is_arithmetic_v<scalar_t>, bool> = true>
    DETRAY_HOST_DEVICE constexpr color& operator*=(const scalar_t scalor) {
        m_data[0] *= scalor;
        m_data[1] *= scalor;
        m_data[2] *= scalor;
        m_data[3] *= scalor;

        return *this;
    }

    /// Scale the color by a normalization factor @param scalor
    DETRAY_HOST_DEVICE
    constexpr color& operator*=(const color& left) {
        m_data[0] *= left.m_data[0];
        m_data[1] *= left.m_data[1];
        m_data[2] *= left.m_data[2];
        m_data[3] *= left.m_data[3];

        return *this;
    }

    /// Mixes two colors @param left and @param right by addition
    DETRAY_HOST_DEVICE
    template <typename>
    friend constexpr color operator+(const color& left, const color& right);

    /// Scale the color by a normalization factor @param scalor
    DETRAY_HOST_DEVICE
    template <typename, typename scalar_t>
    friend constexpr color operator*(const scalar_t scalor, const color& right);

    /// Print the color data to stdout
    DETRAY_HOST
    template <typename>
    friend std::ostream& operator<<(std::ostream& os, const color& c);

    std::array<data_t, 4> m_data{};
};

template <typename data_t>
std::ostream& operator<<(std::ostream& os, const color<data_t>& c) {
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

template <typename data_t>
constexpr color<data_t> operator+(const color<data_t>& left,
                                  const color<data_t>& right) {
    color<data_t> new_color;

    new_color[0] = left[0] + right[0];
    new_color[1] = left[1] + right[1];
    new_color[2] = left[2] + right[2];
    new_color[3] = left[3] + right[3];

    return new_color;
}

template <typename data_t, typename scalar_t>
constexpr color<data_t> operator*(const scalar_t& scalor,
                                  const color<data_t>& right) {
    color<data_t> new_color;

    new_color[0] = static_cast<data_t>(right[0] * scalor);
    new_color[1] = static_cast<data_t>(right[1] * scalor);
    new_color[2] = static_cast<data_t>(right[2] * scalor);
    new_color[3] = static_cast<data_t>(right[3] * scalor);

    return new_color;
}

}  // namespace detray::texture
