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
#include "detray/utils/ranges.hpp"

namespace detray {

/// Calculates the color of a pixel according to surface material
template <typename T, template <typename> class algebra_t>
struct material_shader : public detray::actor {

    template <typename scene_handle_t, typename intersector_state_t>
    DETRAY_HOST_DEVICE void operator()(state &, intersector_state_t &intr_state,
                                       scene_handle_t &sc) const {
        using color_depth =
            typename decltype(sc.m_colors)::value_type::color_depth;

        for (const std::size_t ray_idx :
             detray::views::iota(0ul, sc.rays().size())) {
            if (not intr_state.m_missed[ray_idx]) {
                auto c = texture::detail::material_color_helper<color_depth>(
                    *(intr_state.material()[ray_idx]));

                sc.m_colors[ray_idx] *= c;
            }
        }
    }
};

}  // namespace detray
