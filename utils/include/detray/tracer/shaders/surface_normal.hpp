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
#include "detray/utils/ranges.hpp"

namespace detray {

/// Calculates the color of a pixel according to surface normal
template <typename T, template <typename> class algebra_t>
struct sf_normal_shader : public detray::actor {

    template <typename scene_handle_t, typename intersector_state_t>
    DETRAY_HOST_DEVICE void operator()(state &, intersector_state_t &intr_state,
                                       scene_handle_t &sc) const {
        using color_depth = typename decltype(sc.m_pixel)::color_depth;

        for (const auto &[r_idx, ray] : detray::views::enumerate(sc.rays())) {
            if (intr_state.m_is_inside[r_idx]) {

                const auto &intr = intr_state.m_intersection[r_idx];
                auto normal = sc.geometry().mask()[r_idx].local_frame().normal(
                    sc.geometry().transform()[r_idx], intr.local);

                normal = normal + decltype(normal){1.f, 1.f, 1.f};
                normal = 255.99f * 0.5f * normal;

                // Of the masks that were tested, get the closest one that was
                // hit
                const auto idx = intr_state.closest_solution(r_idx);
                sc.m_pixel += {static_cast<color_depth>(normal[0][idx]),
                               static_cast<color_depth>(normal[1][idx]),
                               static_cast<color_depth>(normal[2][idx])};
            }
        }
    }
};

}  // namespace detray
