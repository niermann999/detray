/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/plugins/tracer/renderer/raw_image.hpp"

// System include(s)
#include <limits>
#include <ratio>

namespace detray {

/// @brief camera class that issues one or more rays per image pixel
template <concepts::algebra algebra_t,
          typename aspect_ratio_t = std::ratio<16, 9>>
class camera {
    using T = dvalue<algebra_t>;
    using scalar_t = dscalar<algebra_t>;
    using point3_t = dpoint3D<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;
    using transform3_t = dtransform3D<algebra_t>;

    public:
    using aspect_ratio = aspect_ratio_t;

    /// Construct the camera using its @param viewport_height , the ray
    /// @param origin and the @param focal_length
    DETRAY_HOST_DEVICE
    constexpr camera(const scalar_t viewport_height = 2.f,
                     const point3_t origin = {0.f, 0.f, 0.f},
                     const scalar_t focal_length = 1.f)
        : m_origin(origin) {
        constexpr T a{static_cast<T>(aspect_ratio::num) /
                      static_cast<T>(aspect_ratio::den)};
        const scalar_t viewport_width{a * viewport_height};

        m_horizontal = {viewport_width, 0.f, 0.f};
        m_vertical = {0.f, viewport_height, 0.f};
        point3_t av = (m_horizontal + m_vertical);
        m_lower_left_corner =
            m_origin - 0.5f * av - vector3_t{0.f, 0.f, focal_length};
    }

    /// @brief Shoot a single ray per pixel
    ///
    /// @param x x coordinate of the pixel
    /// @param y y coordinate of the pixel
    /// @param image the raw image
    ///
    /// @returns a ray, that passes the pixel
    DETRAY_HOST_DEVICE
    template <concepts::arithmetic color_depth>
    constexpr detail::ray<algebra_t> generate_ray(
        const scalar_t x, const scalar_t y,
        const raw_image<color_depth, aspect_ratio> &image) const {

        // percentage of pixel position of the width/height of the image
        const scalar_t u = x * (1.f / static_cast<T>(image.width() - 1u));
        const scalar_t v = y * (1.f / static_cast<T>(image.height() - 1u));

        return {m_origin, m_lower_left_corner + u * m_horizontal +
                              v * m_vertical - m_origin};
    }

    /// @brief Shoot multiple rays to shade a single pixel
    ///
    /// @param x x coordinate of the pixel
    /// @param y y coordinate of the pixel
    /// @param rand_gen random number generator to create direction variation
    /// @param image the raw image
    ///
    /// @returns a range of rays, that pass the pixel in random places
    DETRAY_HOST_DEVICE
    template <std::size_t SAMPLES, concepts::arithmetic color_depth,
              typename generator_t>
    constexpr std::array<detail::ray<algebra_t>, SAMPLES> generate_rays(
        const scalar_t x, const scalar_t y, generator_t &rand_gen,
        const raw_image<color_depth, aspect_ratio> &image) const {

        vector3_t pitch_x = m_horizontal;
        vector3_t pitch_y = m_vertical;
        pitch_x[0] /= image.width();
        pitch_y[1] /= image.height();

        std::array<detail::ray<transform3_t>, SAMPLES> rays;
        for (std::size_t i = 0u; i < SAMPLES; ++i) {
            rays[i] = generate_ray(x, y, image);
            auto &ray = rays[i];

            // Random modification of the ray direction
            const scalar_t px{-0.5f + rand_gen(0.f, 1.f)};
            const scalar_t py{-0.5f + rand_gen(0.f, 1.f)};

            ray.set_dir(ray.dir() + px * pitch_x + py * pitch_y);
        }

        return rays;
    }

    private:
    /// Ray origin
    point3_t m_origin;
    /// Lower left corner of the image
    point3_t m_lower_left_corner;
    /// Horizontal axis of the viewport
    vector3_t m_horizontal;
    /// Vertical axis of the viewport
    vector3_t m_vertical;
};

}  // namespace detray
