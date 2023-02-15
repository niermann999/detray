/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/masks/masks.hpp"

// System include(s)
#include <cassert>
#include <limits>
#include <type_traits>
#include <vector>

namespace detray {

/// An axis aligned bounding box of a given @tparam shape_t
template <typename shape_t, typename scalar_t = scalar>
class axis_aligned_bounding_volume {
    public:
    /// Define geometric properties
    /// @{
    using shape = shape_t;
    using boundaries = typename shape_t::boundaries;
    using axes = typename shape_t::template axes<>;
    template <typename algebra_t>
    using local_frame = typename shape_t::template local_frame_type<algebra_t>;

    static constexpr std::size_t Dim{axes::Dim};
    /// @}

    /// Default constructor builds an infinite box
    constexpr axis_aligned_bounding_volume() = default;

    /// Constructor from mask boundary values
    template <typename... Args>
    DETRAY_HOST_DEVICE explicit constexpr axis_aligned_bounding_volume(
        unsigned int box_id, Args&&... args)
        : m_mask(box_id, std::forward<Args>(args)...) {}

    /// Construct around an arbitrary surface @param mask
    template <
        typename mask_t, typename S = shape_t,
        typename std::enable_if_t<std::is_same_v<S, cuboid3D<>>, bool> = true>
    DETRAY_HOST_DEVICE constexpr axis_aligned_bounding_volume(
        const mask_t& mask, unsigned int box_id, const scalar_t envelope)
        : m_mask{mask.local_min_bounds(envelope).values(), box_id} {
        // Make sure the box is actually 'bounding'
        assert(envelope > (1.f + std::numeric_limits<scalar_t>::epsilon()));
    }

    /// Construct from mask boundary vector
    DETRAY_HOST axis_aligned_bounding_volume(
        const std::vector<scalar_t>& values, unsigned int box_id)
        : m_mask(values, box_id) {
        assert(values.size() == shape::boundaries::e_size &&
               " Given number of boundaries does not match mask shape.");
    }

    /// Construct a bounding box around a set of boxes
    constexpr axis_aligned_bounding_volume(
        const std::vector<axis_aligned_bounding_volume>& aabbs) {}

    /// Subscript operator @returns a single box boundary.
    DETRAY_HOST_DEVICE
    constexpr auto operator[](const std::size_t i) const -> scalar_t {
        return m_mask[i];
    }

    /// @returns the bounds of the box, depending on its shape
    DETRAY_HOST_DEVICE
    constexpr auto id() const -> unsigned int { return m_mask.volume_link(); }

    /// @returns the bounds of the box, depending on its shape
    DETRAY_HOST_DEVICE
    constexpr auto bounds() const -> const mask<shape, unsigned int>& {
        return m_mask;
    }

    /// @brief Lower and upper point for minimum axis aligned bounding box of
    /// cuboid shape.
    ///
    /// Computes the min and max vertices in global coordinates from the local
    /// minimum aabb. The global aabb is not necessarily an minimum aabb.
    ///
    /// @param trf affine transformation
    ///
    /// @returns a new, transformed aabb.
    template <
        typename transform3_t, typename S = shape_t,
        typename std::enable_if_t<std::is_same_v<S, cuboid3D<>>, bool> = true>
    DETRAY_HOST_DEVICE auto transform(const transform3_t& trf) const
        -> axis_aligned_bounding_volume {

        using point3_t = typename transform3_t::point3;

        // The new axis vectors, scaled to the aabb dimensions
        // (e.g. max_x - min_x)
        const auto new_box_x =
            (m_mask[cuboid3D<>::e_max_x] - m_mask[cuboid3D<>::e_min_x]) *
            trf.x();
        const auto new_box_y =
            (m_mask[cuboid3D<>::e_max_y] - m_mask[cuboid3D<>::e_min_y]) *
            trf.y();
        const auto new_box_z =
            (m_mask[cuboid3D<>::e_max_z] - m_mask[cuboid3D<>::e_min_z]) *
            trf.z();

        // Transform the old min and max points to the global frame ...
        point3_t loc_min{m_mask[cuboid3D<>::e_min_x],
                         m_mask[cuboid3D<>::e_min_y],
                         m_mask[cuboid3D<>::e_min_z]};
        point3_t loc_max{m_mask[cuboid3D<>::e_max_x],
                         m_mask[cuboid3D<>::e_max_y],
                         m_mask[cuboid3D<>::e_max_z]};

        // ... and construct all corner points of the local aabb in global
        // coordinates
        std::array<point3_t, 8> glob_c_points;
        glob_c_points[0] = trf.point_to_global(loc_min);
        glob_c_points[1] = trf.point_to_global(loc_max);
        glob_c_points[2] = glob_c_points[0] + new_box_x;
        glob_c_points[3] = glob_c_points[0] + new_box_y;
        glob_c_points[4] = glob_c_points[0] + new_box_z;
        glob_c_points[5] = glob_c_points[1] - new_box_x;
        glob_c_points[6] = glob_c_points[1] - new_box_y;
        glob_c_points[7] = glob_c_points[1] - new_box_z;

        // Find min/max extent of the local aabb in global coordinates
        constexpr scalar_t inf{std::numeric_limits<scalar_t>::infinity()};
        scalar_t min_x{inf}, min_y{inf}, min_z{inf}, max_x{-inf}, max_y{-inf},
            max_z{-inf};
        for (const auto p : glob_c_points) {
            // Check every coordinate of the point
            min_x = p[0] < min_x ? p[0] : min_x;
            max_x = p[0] > max_x ? p[0] : max_x;
            min_y = p[1] < min_y ? p[1] : min_y;
            max_y = p[1] > max_y ? p[1] : max_y;
            min_z = p[2] < min_z ? p[2] : min_z;
            max_z = p[2] > max_z ? p[2] : max_z;
        }

        // Construct the transformed aabb
        return axis_aligned_bounding_volume{
            m_mask.volume_link(), min_x, min_y, min_z, max_x, max_y, max_z};
    }

    /// Checks whether a point lies inside the box. The point has to be defined
    /// in the coordinate frame that is spanned by the box axes.
    DETRAY_HOST_DEVICE
    template <typename point_t>
    constexpr auto is_inside(
        const point_t& loc_p,
        const scalar_t t = std::numeric_limits<scalar_t>::epsilon()) const
        -> intersection::status {
        return m_mask.is_inside(loc_p, t);
    }

    /// Intersect the box with a ray
    DETRAY_HOST_DEVICE
    template <typename algebra_t>
    constexpr bool intersect(
        const detail::ray<algebra_t>& ray,
        const scalar_t t = std::numeric_limits<scalar_t>::epsilon()) const {
        static_assert(std::is_same_v<shape, cuboid3D<>>,
                      "aabbs are only implemented in cuboid shape for now");
        return m_mask.intersector()(ray, m_mask, t);
    }

    /*DETRAY_HOST_DEVICE
    template<typename algebra_t>
    constexpr auto intersect(
        const axis_aligned_bounding_volume &aabb,
        const scalar_t t = std::numeric_limits<scalar_t>::epsilon()) const
        -> intersection::status {
        return m_mask.is_overlap(aabb.bounds());
    }*/

    /*DETRAY_HOST_DEVICE
    template<typename algebra_t>
    constexpr auto intersect(
        const detail::frustum<algebra_t> frustum,
        const scalar_t t = std::numeric_limits<scalar_t>::epsilon()) const
        -> intersection::status {
        ....
    }*/

    private:
    /// Keeps the box boundary values and id
    mask<shape, unsigned int> m_mask;
};

}  // namespace detray