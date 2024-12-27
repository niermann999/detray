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
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/ranges.hpp"

// Detray tracer include(s)
#include "detray/plugins/tracer/renderer/pipeline.hpp"

// System include(s)
#include <limits>
#include <memory>

namespace detray {

// template<typename T>
// struct vc_aos;

/// Calculates the color of a pixel by ray surface intersection with the scene.
/// Starting point of the shader pipeline
template <typename geometry_t, typename aos_algebra_t, std::size_t SAMPLES = 1u,
          typename color_depth = std::uint8_t>
struct intersector : detray::actor {

    using algebra_t = typename geometry_t::algebra_type;
    using T = dvalue<algebra_t>;
    using scalar_t = dscalar<algebra_t>;

    using surface_t = typename geometry_t::surface_type;
    using intersection_t = intersection2D<surface_t, algebra_t, true>;
    using material_t = material<T>;

    struct state {

        using ray_t = detail::ray<aos_algebra_t>;

        DETRAY_HOST_DEVICE
        state(std::array<ray_t, SAMPLES> &rays, std::uint16_t max_ref = 1u,
              const T min = 0.f, const T max = std::numeric_limits<T>::max())
            : m_rays{&rays}, m_max_reflections{max_ref}, m_interval{min, max} {}

        /// @returns the material encountered in the intersection per ray
        const std::array<const material_t *, SAMPLES> &material() const {
            return m_material;
        }

        /// @returns the rays used for intersection for this pixel
        /// @{
        const std::array<ray_t, SAMPLES> &rays() const { return *m_rays; }
        std::array<ray_t, SAMPLES> &rays() { return *m_rays; }
        /// @}

        /// From potentially multiple intersected surfaces in the intersection,
        /// get the index of the closest one
        constexpr const intersection_t &closest_intersection(
            std::size_t ray_idx) const {
            // AoS
            if constexpr (algebra::concepts::aos<algebra_t>) {
                return m_intersection[ray_idx];
            } else {
                // Should be handled by the surface links
                const auto abs_path = math::abs(m_intersection[ray_idx].path);
                const auto m = ((abs_path == abs_path.min()) &&
                                m_intersection[ray_idx].status);

                return m_intersection[m.firstOne()];
            }
        }

        /// From potentially multiple intersected surfaces in the intersection,
        /// get the index of the closest one
        constexpr std::size_t closest_solution_idx(std::size_t ray_idx) const {
            // AoS
            if constexpr (algebra::concepts::aos<algebra_t>) {
                return ray_idx;
            } else {
                // Should be handled by the surface links
                const auto abs_path = math::abs(m_intersection[ray_idx].path);
                const auto m = ((abs_path == abs_path.min()) &&
                                m_intersection[ray_idx].status);

                return m.firstOne();
            }
        }

        /// From potentialliy multiple intersected surfaces in the intersection,
        /// get the index of the closest one
        constexpr T path_to_closest(std::size_t ray_idx) const {
            // AoS
            if constexpr (algebra::concepts::aos<algebra_t>) {
                return m_intersection[ray_idx].path;
            } else {
                return m_intersection[ray_idx]
                    .path[closest_solution_idx(ray_idx)];
            }
        }

        /// @returns @c true if any ray has hit a surface
        constexpr bool has_hit() const {
            bool is_hit = false;
            for (const bool b : m_is_hit) {
                is_hit |= b;
            }
            return is_hit;
        }

        /// @returns @c true if a ray has reached the maximum number of
        /// reflections
        constexpr bool is_finished(std::size_t ray_idx) const {
            return m_refections[ray_idx] >= m_max_reflections;
        }

        /// @returns @c true if any ray has hit a surface
        constexpr bool is_finished() const {
            bool finished = true;
            for (std::size_t i = 0u; i < SAMPLES; ++i) {
                finished &= is_finished(i);
            }
            return finished;
        }

        /// The rays
        std::array<ray_t, SAMPLES> *m_rays;
        /// Resulting intersection
        std::array<intersection_t, SAMPLES> m_intersection{};
        /// The colors for the rays
        std::array<texture::color<color_depth>, SAMPLES> m_colors;
        /// Pointer to the material of the surface
        std::array<const material_t *, SAMPLES> m_material{nullptr};
        /// Flag to the obseving colorizer/shaders that the surface was hit
        std::array<bool, SAMPLES> m_is_hit{false};
        /// Flag to the obseving colorizer/shaders that this ray is dead
        std::array<std::uint16_t, SAMPLES> m_refections{0u};
        /// Distance interval at which to consider an intersection
        std::array<T, 2> m_interval;
        /// How many reflections are allowed per ray
        std::uint16_t m_max_reflections{1u};
    };

    /// Intersect the ray with the mask. The closest intersection is in front of
    /// the @c m_intersections container
    DETRAY_HOST_DEVICE void operator()(
        state &intrs, const scene_handle<geometry_t, color_depth> &sc) const {
        const geometry_t &geo = sc.geometry();

        // Find the intersection information for every ray
        for (const auto &[r_idx, ray] :
             detray::views::enumerate(intrs.rays())) {
            // This is a background ray: No further intersections needed
            if (intrs.is_finished(r_idx)) {
                continue;
            }

            // Reset to test intersection again
            intrs.m_is_hit[r_idx] = false;
            intrs.m_material[r_idx] = nullptr;

            // Perform the intersection on every surface in the geometry
            for (const auto &[m_idx, mask] :
                 detray::views::enumerate(geo.mask())) {

                const auto mask_idx{static_cast<dindex>(m_idx)};
                const surface_t sf_desc{mask_idx,
                                        {0u, mask_idx},
                                        {0u, 0u},
                                        0u,
                                        surface_id::e_sensitive};

                using shape_t = typename std::decay_t<decltype(mask)>::shape;
                if (place_in_collection(
                        ray_intersector<shape_t, algebra_t, true>{}(
                            ray, sf_desc, mask, geo.transform()[m_idx],
                            darray<scalar_t, 2u>{0.00001f, 0.00001f}),
                        intrs.m_intersection[r_idx])) {

                    intrs.m_is_hit[r_idx] = true;

                    // Prefix sum to get the material
                    constexpr std::size_t simd_size{sizeof(scalar_t) /
                                                    sizeof(T)};
                    const std::size_t intr_idx{
                        intrs.closest_solution_idx(r_idx)};
                    const auto mat_idx = static_cast<dindex>(m_idx * simd_size);
                    intrs.m_intersection[r_idx].sf_desc.update_material(
                        mat_idx);
                    intrs.m_material[r_idx] =
                        std::addressof(geo.material()[mat_idx]);
                }
            }
            // Count reflections
            ++intrs.m_refections[r_idx];
        }
    }

    private:
    /// @brief check if the intersection is the next candidate
    DETRAY_HOST_DEVICE
    constexpr bool is_candidate(const intersection_t &new_sfi,
                                const intersection_t &current_sfi,
                                const scalar_t min_path = 0.01f) const {
        return detail::any_of(
            (math::abs(new_sfi.path) > min_path) &&
            (math::abs(new_sfi.path) < math::abs(current_sfi.path)) &&
            new_sfi.status);
    }

    /// Places the single solution of a ray-surface intersection @param sfi
    /// in the given container @param intersections, if the surfaces was hit.
    ///
    /// @returns @c true if the intersection was is valid.
    DETRAY_HOST_DEVICE bool place_in_collection(
        intersection_t &&sfi, intersection_t &intersection) const {

        if (is_candidate(sfi, intersection)) {
            intersection = std::move(sfi);
            return true;
        } else {
            return false;
        }
    }

    /// Places all of those solutions of a ray-surface intersection @param sfi
    /// in the given container @param intersections, that hit the surface
    ///
    /// @returns @c true if at least one valid intersection solution was found.
    DETRAY_HOST_DEVICE bool place_in_collection(
        std::array<intersection_t, 2> &&solutions,
        intersection_t &intersection) const {

        bool is_valid = false;

        for (auto &sfi : solutions) {

            if (is_candidate(sfi, intersection)) {
                intersection = std::move(sfi);
                is_valid = true;
            }
        }
        return is_valid;
    }
};

}  // namespace detray
