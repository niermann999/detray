/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/plugins/tracer/texture/color.hpp"
#include "detray/utils/concepts.hpp"

// System include(s)
#include <array>
#include <ostream>
#include <type_traits>

namespace detray::texture {

namespace detail {

/// @brief holds pixel coordinates and its color.
///
/// @tparam data_t pixel coordinate type
template <unsigned int D, concepts::arithmetic depth = std::uint8_t,
          std::integral data_t = unsigned int>
struct pixelD {

    using color_t = texture::color<depth>;
    using color_depth = depth;

    static constexpr unsigned int Dim{D};

    /// Default constructor
    constexpr pixelD() = default;

    /// Construct from an array of coordinates @param coord
    DETRAY_HOST_DEVICE
    constexpr pixelD(const std::array<data_t, D>& coord) : m_coord{coord} {}

    /// Construct from its coordinates @param coord
    template <std::integral... C>
    DETRAY_HOST_DEVICE constexpr pixelD(const C... coord) : m_coord{coord...} {}

    /// Construct from an array of coordinates @param coord and a color @param c
    DETRAY_HOST_DEVICE
    constexpr pixelD(const color_t& c, const std::array<data_t, D>& coord)
        : m_coord{coord}, m_color{c} {}

    /// Construct from its coordinates @param coord and a color @param c
    template <std::integral... C>
    constexpr pixelD(const color_t& c, const C... coord)
        : m_coord{coord...}, m_color{c} {}

    template <concepts::arithmetic other_depth_t>
    requires std::is_convertible_v<depth, other_depth_t>
        DETRAY_HOST_DEVICE constexpr operator pixelD<D, data_t, other_depth_t>()
            const {
        return pixelD<D, data_t, other_depth_t>{
            m_coord, static_cast<texture::color<other_depth_t>>(m_color)};
    }

    /// Equality operator: Only considers exact match
    DETRAY_HOST_DEVICE
    constexpr data_t operator==(const pixelD& other) {
        return (m_coord == other.m_coord) && (m_color == other.m_color);
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
    constexpr color_t color() const { return m_color; }

    /// @returns the color of the pixel
    DETRAY_HOST_DEVICE
    constexpr color_t& color() { return m_color; }

    /// Set the color of the pixel to @param c
    DETRAY_HOST_DEVICE
    constexpr void set_color(const color_t& c) { m_color = c; }

    /// Mixes the pixel color by addition
    DETRAY_HOST_DEVICE
    constexpr pixelD operator+=(const color_t& c) {
        m_color += c;
        return *this;
    }

    /// Scale the pixel color
    DETRAY_HOST_DEVICE
    template <concepts::arithmetic scalar_t>
    constexpr pixelD operator*=(const scalar_t factor) {
        m_color *= factor;
        return *this;
    }

    /// Mixes the pixel color by multiplication
    DETRAY_HOST_DEVICE
    constexpr pixelD operator*=(const color_t& c) {
        m_color *= c;
        return *this;
    }

    /// Print the pixel data to stdout
    DETRAY_HOST
    template <unsigned int, typename, typename>
    friend std::ostream& operator<<(std::ostream&, const pixelD&);

    std::array<data_t, D> m_coord{};
    color_t m_color{};
};

template <unsigned int D, typename depth, typename data_t>
std::ostream& operator<<(std::ostream& os, const pixelD<D, depth, data_t>& px) {
    return os << "pix: " << static_cast<unsigned int>(px[0]) << ", "
              << static_cast<unsigned int>(px[1]) << ", " << px.color();
}

}  // namespace detail

template <typename depth = std::uint8_t, typename data_t = unsigned int>
using pixel = detail::pixelD<2, depth, data_t>;

template <typename depth = std::uint8_t, typename data_t = unsigned int>
using voxel = detail::pixelD<3, depth, data_t>;

}  // namespace detray::texture
