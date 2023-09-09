/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/propagator/actor_chain.hpp"
#include "detray/tracer/renderer/raw_image.hpp"
#include "detray/tracer/renderer/single_shape.hpp"
#include "detray/tracer/texture/pixel.hpp"

namespace detray {

/// Executes the rendering steps sequentially
template <typename intersector_t, typename background_shader_t,
          typename... shaders_t>
using rendering_pipeline =
    actor_chain<dtuple, composite_actor<dtuple, intersector_t,
                                        background_shader_t, shaders_t...>>;

/// State that is passed through the pipeline per ray
///
/// Contains a pointer to the geometry and the ray that this pipeline instance
/// renders.
struct scene_handle {

    struct config {};

    template <typename geometry_t, typename color_depth, std::size_t SAMPLES>
    struct state {

        using transform3D = typename geometry_t::transform3D;
        using ray_t = detail::ray<dtransform3D<ALGEBRA_PLUGIN<detray::scalar>>>;

        DETRAY_HOST_DEVICE
        state(const geometry_t &geo, const raw_image<color_depth> &im,
              const std::array<
                  detail::ray<dtransform3D<ALGEBRA_PLUGIN<detray::scalar>>>,
                  SAMPLES> &ray)
            : m_geo{&geo}, m_image{&im}, m_ray{&ray} {
            m_colors.fill(texture::color<detray::scalar>{1});
        }

        /// Threadsafe interface
        /// @{
        const std::array<ray_t, SAMPLES> &rays() const { return *m_ray; }
        const geometry_t &geometry() const { return *m_geo; }
        /// @}

        /// The geometry handle
        const geometry_t *m_geo;
        /// The image handle
        const raw_image<color_depth> *m_image;
        /// The ray handle
        const std::array<ray_t, SAMPLES> *m_ray;
        /// The color for this ray (save with enough memory for color mixing)
        std::array<texture::color<detray::scalar>, SAMPLES> m_colors;
    };

#if __clang__
    template <typename geometry_t, typename color_depth>
    DETRAY_HOST_DEVICE state(
        const geometry_t &geo, const raw_image<color_depth> &im,
        const detail::ray<dtransform3D<ALGEBRA_PLUGIN<detray::scalar>>> &ray)
        -> state<geometry_t, color_depth>;
#endif
};

}  // namespace detray
