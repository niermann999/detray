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

    using link_t = dsimd<algebra_t, std::uint_least16_t>;
    using surface_t = surface_descriptor<dtyped_index<>, dtyped_index<>, link_t,
                                         link_t, link_t>;
    using intersection_t = intersection2D<surface_t, T, algebra_t>;

    struct global_state {
        using scalar_t = dscalar<algebra_t<T>>;
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
        std::size_t closest_solution(std::size_t ray_idx) const {
            // AoS?
            if constexpr (std::is_same_v<
                              decltype(m_intersection[ray_idx].status), bool>) {
                return 0;
            } else {
                // Should be handled by the surface links
                return m_intersection[ray_idx].status.firstOne();
            }
        }

        /// @returns @c true if any ray has hit a surface
        bool has_hit() const {
            // AoS?
            bool is_hit = false;
            for (std::size_t i = 0u; i < SAMPLES; ++i) {
                is_hit &= m_missed[i];
            }
            return is_hit;
        }

        std::array<T, 2> m_interval;
        /// Resulting intersection
        std::array<intersection_t, SAMPLES> m_intersection{};
        /// Pointer to the material of the surface
        std::array<const material_t *, SAMPLES> m_material{nullptr};
        /// Flag to the obseving colorizer/shaders that the surface was hit
        std::array<bool, SAMPLES> m_missed{false};
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
            if (loc_st.m_missed[r_idx]) {
                return;
            }

            // Reset to test intersection again
            loc_st.m_missed[r_idx] = true;
            loc_st.m_material[r_idx] = nullptr;

            // Perform the intersection on every mask in the geometry
            for (const auto &[m_idx, mask] :
                 detray::views::enumerate(geo.mask())) {

                if (place_in_collection(
                        mask.template intersector<intersection_t>()(
                            ray, surface_t{}, mask, geo.transform()[m_idx]),
                        loc_st.m_intersection[r_idx])) {

                    loc_st.m_missed[r_idx] = false;

                    const auto mat_idx = m_idx + loc_st.closest_solution(r_idx);
                    loc_st.m_material[r_idx] =
                        std::addressof(geo.material()[mat_idx]);
                }
            }
        }
    }

    private:
    /// Places the single solution of a ray-surface intersection @param sfi
    /// in the given container @param intersections, if the surfaces was hit.
    ///
    /// @returns @c true if the intersection was is valid.
    DETRAY_HOST_DEVICE bool place_in_collection(
        intersection_t &&sfi, intersection_t &intersection) const {
        if (detail::any_of(
                (math_ns::abs(sfi.path) < math_ns::abs(intersection.path)) &&
                sfi.status)) {
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
            if (detail::any_of((math_ns::abs(sfi.path) <
                                math_ns::abs(intersection.path)) &&
                               sfi.status)) {
                intersection = std::move(sfi);
                is_valid = true;
            }
        }
        return is_valid;
    }
};

}  // namespace detray
