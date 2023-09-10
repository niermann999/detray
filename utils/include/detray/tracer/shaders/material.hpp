/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include <iostream>

#include "detray/definitions/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/tracer/definitions/colors.hpp"
#include "detray/tracer/shaders/detail/random_scattering.hpp"
#include "detray/tracer/texture/color.hpp"
#include "detray/tracer/texture/detail/material_color_helper.hpp"
#include "detray/utils/ranges.hpp"

namespace detray {

/// Calculates the color of a pixel according to surface material
template <typename T, template <typename> class algebra_t, typename generator_t>
struct material_shader : public detray::actor {

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
    DETRAY_HOST_DEVICE inline void operator()(
        state &mat_state, const intersector_state_t &intr_state,
        scene_handle_t &sc) const {
        using color_depth =
            typename decltype(sc.m_colors)::value_type::color_depth;

        std::size_t r_idx{0u};
        for (auto &ray : sc.rays()) {

            if (!intr_state.m_finished[r_idx] and !intr_state.m_missed[r_idx]) {

                auto c = texture::detail::material_color_helper<color_depth>(
                    *(intr_state.material()[r_idx]));

                random_scattering<T, ALGEBRA_PLUGIN>(
                    ray, intr_state.path_to_closest(r_idx),
                    mat_state.random_gen(), mat_state.min(), mat_state.max());

                sc.m_colors[r_idx] *= c;
            }
            ++r_idx;
        }
    }
};

}  // namespace detray
