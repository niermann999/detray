/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/ranges.hpp"

// Detray tracer include(s)
#include "detray/plugins/tracer/definitions/colors.hpp"
#include "detray/plugins/tracer/texture/color.hpp"
#include "detray/plugins/tracer/texture/pixel.hpp"

namespace detray {

/// @brief Image background class tag
struct image_background {};

/// @brief Single color image background
template <concepts::algebra algebra_t>
struct plain_background : public image_background {

    using algebra_type = algebra_t;

    /// Calculate the pixel color
    template <typename color_depth = std::uint8_t>
    DETRAY_HOST_DEVICE inline static consteval texture::color<color_depth> get(
        const detail::ray<algebra_t> &) {
        return m_color<color_depth>;
    }

    template <typename color_depth = std::uint8_t>
    static constexpr auto m_color = texture::white<color_depth>;
};

/// @brief Gradient background as described in
template <concepts::algebra algebra_t>
struct gradient_background : public image_background {

    using algebra_type = algebra_t;
    using point3_t = dpoint3D<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    /// Calculate the pixel color
    template <typename color_depth = std::uint8_t>
    DETRAY_HOST_DEVICE inline static constexpr texture::color<color_depth> get(
        const detail::ray<algebra_t> &ray) {
        const vector3_t dir = ray.dir();
        point3_t p1{1.0f, 1.0f, 1.0f};
        point3_t p2{0.85f, 0.85f, 1.0f};
        const auto t = 0.5f * dir[1] + 1.0f;
        point3_t p4 = ((1.0f - t) * p1 + t * p2);
        point3_t p3 = p4;
        if constexpr (std::is_integral_v<color_depth>) {
            p3 = static_cast<float>(texture::color<color_depth>::channel_max) *
                 p4;
        }

        return {
            static_cast<color_depth>(p3[0]), static_cast<color_depth>(p3[1]),
            static_cast<color_depth>(p3[2]), static_cast<color_depth>(255u)};
    }
};

/// @brief Gradient background as described in
template <class image_background_t>
struct inf_plane : public image_background {

    using algebra_type = typename image_background_t::algebra_type;

    /// Calculate the pixel color
    template <typename color_depth = std::uint8_t>
    DETRAY_HOST_DEVICE inline static constexpr texture::color<color_depth> get(
        const detail::ray<algebra_type> &ray) {

        // Flip background color at y = 0
        if (!detail::any_of(
                math::signbit(ray.pos()[1] + 0.1f * ray.dir()[1]))) {
            return texture::green<color_depth>;
        } else {
            return image_background_t::template get<color_depth>(ray);
        }
    }
};

/// Calculates the color of a background pixel using different backgound types
template <class image_background_t>
struct background_shader : public detray::actor {

    /// Set the pixel color
    template <typename scene_handle_t, typename intersector_state_t>
    DETRAY_HOST_DEVICE inline void operator()(intersector_state_t &intr_state,
                                              scene_handle_t &) const {
        using color_depth = typename scene_handle_t::color_depth;

        // Set only pixels that are not part of an object in the scene
        for (const auto &[r_idx, ray] :
             detray::views::enumerate(intr_state.rays())) {
            if (!intr_state.is_finished(r_idx) && !intr_state.m_is_hit[r_idx]) {

                auto c = image_background_t::template get<color_depth>(ray);

                // First ray?
                if (intr_state.m_refections[r_idx] == 1u) {
                    intr_state.m_colors[r_idx] = c;
                }
                // Ray went to infinity
                intr_state.m_refections[r_idx] = intr_state.m_max_reflections;
            }
        }
    }
};

}  // namespace detray
