/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <iostream>
#include <limits>
#include <memory>

namespace detray {

/// Calculates the color of a pixel. Starting point of the shader pipeline
template <typename T, template <typename> class algebra_t, std::size_t SAMPLES,
          typename mask_t, typename material_t = material_slab<T>>
struct single_shape : detray::actor {

    using scalar_t = dscalar<algebra_t<T>>;
    using link_t = dsimd<algebra_t, std::uint_least16_t>;
    using surface_t = surface_descriptor<dtyped_index<>, dtyped_index<>, link_t,
                                         link_t, link_t>;
    using intersection_t = intersection2D<surface_t, T, algebra_t>;

    struct global_state {
        using point3D = dpoint3D<algebra_t<T>>;
        using vector3D = dvector3D<algebra_t<T>>;
        using transform3D = dtransform3D<algebra_t<T>>;

        /// Construct from surface data:
        DETRAY_HOST_DEVICE
        global_state(const std::vector<transform3D> &trf,
                     std::vector<mask_t> &&mask,
                     const std::vector<material_t> &mat)
            : m_trf{std::move(trf)},
              m_mask{std::move(mask)},
              m_material{std::move(mat)} {}

        /// Threadsafe interface
        /// @{
        const std::vector<transform3D> &transform() const { return m_trf; }
        const std::vector<mask_t> &mask() const { return m_mask; }
        const std::vector<material_t> &material() const { return m_material; }
        /// @}

        /// The surfaces data
        std::vector<transform3D> m_trf;
        std::vector<mask_t> m_mask;
        std::vector<material_t> m_material;
    };

    struct state {

        DETRAY_HOST_DEVICE
        state(const T min = 0.f, const T max = std::numeric_limits<T>::max())
            : m_interval{min, max} {}

        const std::array<const material_t *, SAMPLES> &material() const {
            return m_material;
        }

        /// From potentialliy multiple intersected surfaces in the intersection,
        /// get the index of the closest one
        inline constexpr std::size_t closest_solution(
            std::size_t ray_idx) const {
            // AoS?
            if constexpr (std::is_same_v<
                              decltype(m_intersection[ray_idx].status), bool>) {
                return 0;
            } else {
                // Should be handled by the surface links
                const auto abs_path =
                    math_ns::abs(m_intersection[ray_idx].path);
                const auto m = ((abs_path == abs_path.min()) &&
                                m_intersection[ray_idx].status);

                return m.firstOne();
            }
        }

        /// From potentialliy multiple intersected surfaces in the intersection,
        /// get the index of the closest one
        inline constexpr T path_to_closest(std::size_t ray_idx) const {
            // AoS?
            if constexpr (std::is_same_v<
                              decltype(m_intersection[ray_idx].status), bool>) {
                return m_intersection[ray_idx].path;
            } else {
                return m_intersection[ray_idx].path[closest_solution(ray_idx)];
            }
        }

        /// @returns @c true if any ray has hit a surface
        inline constexpr bool has_hit() const {
            bool is_hit = false;
            for (const bool b : m_is_hit) {
                is_hit |= b;
            }
            return is_hit;
        }

        std::array<T, 2> m_interval;
        /// Resulting intersection
        std::array<intersection_t, SAMPLES> m_intersection{};
        /// Pointer to the material of the surface
        std::array<const material_t *, SAMPLES> m_material{nullptr};
        /// Flag to the obseving colorizer/shaders that the surface was hit
        std::array<bool, SAMPLES> m_is_hit{false};
        /// Flag to the obseving colorizer/shaders that the surface was hit
        std::array<bool, SAMPLES> m_finished{false};
    };

    /// Intersect the ray with the mask. The closest intersection is in front of
    /// the @c m_intersections container
    template <typename scene_handle_t>
    DETRAY_HOST_DEVICE void operator()(state &loc_st,
                                       const scene_handle_t &sc) const {
        // In this special case, the geometry will be this actor's global_state
        const global_state &geo = sc.geometry();

        // Find the intersection information for every ray
        for (const auto &[r_idx, ray] : detray::views::enumerate(sc.rays())) {

            // This is a background ray: No further intersections needed
            if (loc_st.m_finished[r_idx]) {
                continue;
            }

            // Reset to test intersection again
            loc_st.m_is_hit[r_idx] = false;
            loc_st.m_material[r_idx] = nullptr;

            // Perform the intersection on every mask in the geometry
            for (const auto &[m_idx, mask] :
                 detray::views::enumerate(geo.mask())) {

                if (place_in_collection(
                        mask.template intersector<intersection_t>()(
                            ray, surface_t{}, mask, geo.transform()[m_idx]),
                        loc_st.m_intersection[r_idx])) {

                    loc_st.m_is_hit[r_idx] = true;

                    // Prefix sum to get the material
                    constexpr std::size_t simd_size{sizeof(scalar_t) /
                                                    sizeof(T)};
                    const std::size_t intr_idx{loc_st.closest_solution(r_idx)};
                    const auto mat_idx = m_idx * simd_size + intr_idx;
                    loc_st.m_material[r_idx] =
                        std::addressof(geo.material()[mat_idx]);
                }
            }
        }
    }

    private:
    /// @brief check if the intersection is the next candidate
    DETRAY_HOST_DEVICE
    inline constexpr bool is_candidate(const intersection_t &new_sfi,
                                       const intersection_t &current_sfi,
                                       const scalar_t min_path = 0.01f) const {
        return detail::any_of(
            (math_ns::abs(new_sfi.path) > min_path) &&
            (math_ns::abs(new_sfi.path) < math_ns::abs(current_sfi.path)) &&
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
