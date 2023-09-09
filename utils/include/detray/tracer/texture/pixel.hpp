/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/tracer/texture/color.hpp"

// System include(s)
#include <array>
#include <ostream>
#include <type_traits>

namespace detray::texture {

namespace detail {

/// @brief holds pixel coordinates and its color.
///
/// @tparam data_t how to store the pixel data.
template <uint D, typename data_t = uint, typename depth = uint8_t>
struct pixelD {

    using color_t = texture::color<depth>;
    using color_depth = depth;

    static constexpr uint Dim{D};

    /// Default constructor
    constexpr pixelD() = default;

    /// Construct from ist coordinates @param coord
    DETRAY_HOST_DEVICE
    constexpr pixelD(const std::array<data_t, D>& coord)
        : m_color{}, m_coord{coord} {}

    /// Construct from its coordinates @param coord and a color @param c
    DETRAY_HOST_DEVICE
    constexpr pixelD(const std::array<data_t, D>& coord, const color_t& c)
        : m_color{c}, m_coord{coord} {}

    template <typename other_depth_t,
              std::enable_if_t<std::is_convertible_v<depth, other_depth_t>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr operator pixelD<D, data_t, other_depth_t>()
        const {
        return pixelD<D, data_t, other_depth_t>{
            m_coord, static_cast<texture::color<other_depth_t>>(m_color)};
    }

    /// Equality operator: Only considers exact match
    DETRAY_HOST_DEVICE
    constexpr data_t operator==(const pixelD& other) {
        return (m_coord == other.m_coord) and (m_color == other.m_color);
    }

    /// Subscript operator @returns the pixel coordinates - const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) operator[](const std::size_t i) const {
        return m_coord[i];
    }

    /// Subscript operator @returns the pixel coordinates - non-const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) operator[](const std::size_t i) {
        return m_coord[i];
    }

    /// @returns the color of the pixel
    DETRAY_HOST_DEVICE
    constexpr const color_t& color() const { return m_color; }

    /// @returns the color of the pixel
    DETRAY_HOST_DEVICE
    constexpr color_t& color() { return m_color; }

    /// Set the color of the pixel to @param c
    DETRAY_HOST_DEVICE
    constexpr void set_color(const color_t& c) { m_color = c; }

    /// Mixes by adding their colors
    DETRAY_HOST_DEVICE
    constexpr pixelD operator+=(const color_t& c) {
        m_color += c;
        return *this;
    }

    /// Scale the pixel color
    DETRAY_HOST_DEVICE
    template <typename scalar_t,
              std::enable_if_t<std::is_arithmetic_v<scalar_t>, bool> = true>
    constexpr pixelD operator*=(const scalar_t scalor) {
        m_color *= scalor;
        return *this;
    }

    /// Print the pixel data to stdout
    DETRAY_HOST
    template <uint, typename, typename>
    friend std::ostream& operator<<(std::ostream&, const pixelD&);

    color_t m_color{};
    std::array<data_t, D> m_coord{};
};

template <uint D, typename data_t, typename depth>
std::ostream& operator<<(std::ostream& os, const pixelD<D, data_t, depth>& px) {
    return os << "pix: " << static_cast<uint>(px[0]) << ", "
              << static_cast<uint>(px[1]) << ", " << px.color();
}

}  // namespace detail

template <typename data_t = uint, typename depth = uint8_t>
using pixel = detail::pixelD<2, data_t, depth>;

template <typename data_t = uint, typename depth = uint8_t>
using voxel = detail::pixelD<3, data_t, depth>;

}  // namespace detray::texture
