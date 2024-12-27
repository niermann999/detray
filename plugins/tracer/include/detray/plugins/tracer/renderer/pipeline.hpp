/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/propagator/actor_chain.hpp"

// Detray tracer include(s)
#include "detray/plugins/tracer/renderer/raw_image.hpp"
#include "detray/plugins/tracer/texture/pixel.hpp"

namespace detray {

/// Executes the rendering steps sequentially
template <typename intersector_t, typename background_shader_t,
          typename... shaders_t>
using rendering_pipeline = actor_chain<
    composite_actor<intersector_t, background_shader_t, shaders_t...>>;

/// @brief Global state that is passed through the pipeline
///
/// Contains a pointer to the geometry and the ray that this pipeline instance
/// renders.
template <typename geometry_t, typename color_depth_t>
struct scene_handle {

    struct config {};

    using algebra_t = typename geometry_t::algebra_type;
    using color_depth = color_depth_t;

    DETRAY_HOST_DEVICE
    scene_handle(const geometry_t &geo, const raw_image<color_depth> &im)
        : m_geo{&geo}, m_image{&im} {}

    /// Threadsafe interface
    /// @{
    const geometry_t &geometry() const { return *m_geo; }
    const raw_image<color_depth> &image() const { return *m_image; }
    raw_image<color_depth> &image() { return *m_image; }
    /// @}

    /// The geometry handle
    const geometry_t *m_geo;
    /// The image handle
    const raw_image<color_depth> *m_image;
};

}  // namespace detray
