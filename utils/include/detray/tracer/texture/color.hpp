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
#include <iostream>
#include <limits>

namespace detray::texture {

/// @brief holds rgb and alpha values for color shading
///
/// @tparam data_t how to store the rgb data: single color or soa
template <typename data_t = std::uint8_t>
struct color {

    static constexpr auto max_I{std::numeric_limits<data_t>::max()};

    /// Default constructor
    constexpr color() = default;

    /// Construct from colors @param r (red), @param g (green), @param b (blue)
    /// and @param alpha values
    DETRAY_HOST_DEVICE
    constexpr color(const data_t r, const data_t g, const data_t b,
                    const data_t alpha = {max_I - 1u})
        : m_data{r, g, b, alpha} {}

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

    /// Mixes two colors @param left and @param right by addition
    DETRAY_HOST_DEVICE
    template <typename>
    friend constexpr color operator+(const color& left, const color& right);

    /// Scale the color by a normalization factor @param scalor
    DETRAY_HOST_DEVICE
    template <typename, typename scalar_t>
    friend constexpr color operator*(const color& left, const scalar_t scalor);

    /// Print the color data to stdout
    DETRAY_HOST
    template <typename>
    friend std::ostream& operator<<(std::ostream& os, const color& c);

    std::array<data_t, 4> m_data{};
};

template <typename data_t>
std::ostream& operator<<(std::ostream& os, const color<data_t>& c) {
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
constexpr color<data_t> operator*(const color<data_t>& left,
                                  const scalar_t& scalor) {
    color<data_t> new_color;

    new_color[0] = static_cast<data_t>(left[0] * scalor);
    new_color[1] = static_cast<data_t>(left[1] * scalor);
    new_color[2] = static_cast<data_t>(left[2] * scalor);
    new_color[3] = static_cast<data_t>(left[3] * scalor);

    return new_color;
}

}  // namespace detray::texture
