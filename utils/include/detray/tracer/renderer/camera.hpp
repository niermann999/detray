/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/tracer/renderer/raw_image.hpp"

// System include(s)
#include <limits>
#include <ratio>

namespace detray {

template <typename T, template <typename> class algebra_t,
          typename aspect_ratio = std::ratio<16, 9>>
class camera {
    using scalar_t = dscalar<algebra_t<T>>;
    using point3D = dpoint3D<algebra_t<T>>;
    using vector3D = dvector3D<algebra_t<T>>;
    using transform3D = dtransform3D<algebra_t<T>>;

    public:
    DETRAY_HOST_DEVICE
    camera(const scalar_t viewport_height = 2.0f,
           const point3D origin = {0.0f, 0.0f, 0.0f},
           const scalar_t focal_length = 1.0f)
        : m_origin(origin) {
        constexpr T a{static_cast<T>(aspect_ratio::num) /
                      static_cast<T>(aspect_ratio::den)};
        const scalar_t viewport_width{a * viewport_height};

        m_horizontal = point3D{viewport_width, 0.f, 0.f};
        m_vertical = point3D{0.f, viewport_height, 0.f};
        point3D av = (m_horizontal + m_vertical);
        m_lower_left_corner =
            m_origin - 0.5f * av - vector3D{0.f, 0.f, focal_length};
    }

    /// @brief Shoot a single ray per pixel
    ///
    /// @param x x coordinate of the pixel
    /// @param y y coordinate of the pixel
    /// @param image the raw image
    ///
    /// @returns a ray, that passes the pixel
    DETRAY_HOST_DEVICE
    template <typename color_depth>
    constexpr detail::ray<transform3D> get_ray(
        const scalar_t x, const scalar_t y,
        const raw_image<color_depth> &image) const {

        // percentage of pixel position of the width/height of the image
        const scalar_t u = x / static_cast<T>(image.width() - 1u);
        const scalar_t v = y / static_cast<T>(image.height() - 1u);

        detail::ray<transform3D> ray{
            m_origin, 0.f,
            m_lower_left_corner + u * m_horizontal + v * m_vertical - m_origin,
            0.f};
        ray.set_overstep_tolerance(-std::numeric_limits<T>::max());

        return ray;
    }

    /// @brief Shoot multiple rays to shade a single pixel
    ///
    /// @param x x coordinate of the pixel
    /// @param y y coordinate of the pixel
    /// @param rand_x x coordinate of the pixel
    /// @param rand_y y coordinate of the pixel
    /// @param image the raw image
    ///
    /// @returns a range of rays, that pass the pixel and its neighbors in
    /// random places
    DETRAY_HOST_DEVICE
    template <std::size_t SAMPLES, typename color_depth, typename generator_t>
    constexpr std::array<detail::ray<transform3D>, SAMPLES> get_rays(
        const scalar_t x, const scalar_t y,
        [[maybe_unused]] generator_t &rand_gen,
        const raw_image<color_depth> &image) const {

        if constexpr (SAMPLES == 1ul) {
            return {get_ray(x, y, image)};
        } else {

            vector3D pitch_x = m_horizontal;
            vector3D pitch_y = m_vertical;
            pitch_x[0] /= image.width();
            pitch_y[1] /= image.height();

            std::array<detail::ray<transform3D>, SAMPLES> rays;
            for (std::size_t i_s = 0u; i_s < SAMPLES; ++i_s) {
                rays[i_s] = get_ray(x, y, image);
                auto &ray = rays[i_s];

                // Random modification of the ray direction
                const scalar_t px{-0.5f + rand_gen(0.f, 1.f)};
                const scalar_t py{-0.5f + rand_gen(0.f, 1.f)};

                ray.set_dir(ray.dir() + px * pitch_x + py * pitch_y);
            }

            return rays;
        }
    }

    private:
    point3D m_origin;
    point3D m_lower_left_corner;
    vector3D m_horizontal;
    vector3D m_vertical;
};

}  // namespace detray
