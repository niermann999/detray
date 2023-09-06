/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/tracer/definitions/colors.hpp"
#include "detray/tracer/texture/color.hpp"
#include "detray/tracer/texture/detail/material_color_helper.hpp"

namespace detray {

/// Calculates the color of a pixel according to surface material
template <typename T, template <typename> class algebra_t>
struct material_shader : public detray::actor {

    template <typename scene_handle_t, typename intersector_state_t>
    DETRAY_HOST_DEVICE void operator()(state &, intersector_state_t &intr_state,
                                       scene_handle_t &sc) const {
        using color_depth = typename decltype(sc.m_pixel)::color_depth;

        if (intr_state.m_is_inside) {
            auto c = texture::detail::material_color_helper<color_depth>(
                intr_state.material());

            sc.m_pixel.set_color(c);
        }
    }
};

}  // namespace detray
