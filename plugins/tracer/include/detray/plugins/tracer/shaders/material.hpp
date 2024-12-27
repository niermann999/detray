/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/ranges.hpp"

// Detray tracer include(s)
#include "detray/plugins/tracer/definitions/colors.hpp"
#include "detray/plugins/tracer/shaders/detail/random_scattering.hpp"
#include "detray/plugins/tracer/texture/color.hpp"
#include "detray/plugins/tracer/texture/detail/material_color_map.hpp"

namespace detray {

/// Calculates the color of a pixel according to surface material
template <concepts::algebra algebra_t, typename generator_t>
struct material_shader : public detray::actor {

    using T = dvalue<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    struct state {

        DETRAY_HOST_DEVICE
        state(generator_t &rand_gen, const T min = 0.f, const T max = 1.f)
            : m_gen{&rand_gen}, m_min{min}, m_max{max} {}

        generator_t &random_gen() { return *m_gen; }
        T min() const { return m_min; }
        T max() const { return m_max; }

        generator_t *m_gen{nullptr};
        T m_min{0.f}, m_max{1.f};
    };

    template <typename scene_handle_t, typename intersector_state_t>
    DETRAY_HOST_DEVICE inline void operator()(state &mat_state,
                                              intersector_state_t &intr_state,
                                              scene_handle_t &sc) const {

        using color_depth = typename scene_handle_t::color_depth;

        std::size_t r_idx{0u};
        for (auto &ray : intr_state.rays()) {

            if (!intr_state.is_finished(r_idx) && intr_state.m_is_hit[r_idx]) {

                // Add the material color
                auto c = texture::detail::material_color_map<color_depth>(
                    *(intr_state.material()[r_idx]));

                // First ray?
                if (intr_state.m_refections[r_idx] == 1u) {
                    intr_state.m_colors[r_idx] = c;
                } else {
                    // intr_state.m_colors[r_idx] *= c;
                    intr_state.m_colors[r_idx] += (c * 0.1f);
                    intr_state.m_colors[r_idx] *= 0.5f;
                }

                // Scatter the ray
                const auto &intr = intr_state.closest_intersection(r_idx);

                const auto &trf =
                    sc.geometry().transform().at(intr.sf_desc.transform());
                const auto &mask =
                    sc.geometry().mask()[intr.sf_desc.mask().index()];

                const auto normals =
                    mask.get_local_frame().normal(trf, intr.local);

                random_scattering(ray, normals, intr.path,
                                  mat_state.random_gen(), mat_state.min(),
                                  mat_state.max());
            }
            ++r_idx;
        }
    }
};

}  // namespace detray
