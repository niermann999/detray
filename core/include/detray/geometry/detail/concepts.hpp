/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/detail/type_traits.hpp"

namespace detray::concepts {

/// Geometric objects with a planar geometry (rectangles, disks etc)
template <typename T>
concept planar =
    detail::is_planar_v<T> || detail::is_planar_v<typename T::frame_type>;

/// Geometric objects with a cylindrical geometry (2D/3D cylinders)
template <typename T>
concept cylindrical = detail::is_cylindrical_v<T> ||
                      detail::is_cylindrical_v<typename T::frame_type>;

/// Geometric objects with a line geometry (straw tubes, wire chamber cells)
template <typename T>
concept linear =
    detail::is_linear_v<T> || detail::is_linear_v<typename T::frame_type>;

/// Geometric objects with a spherical geometry
template <typename T>
concept spherical =
    detail::is_spherical_v<T> || detail::is_spherical_v<typename T::frame_type>;

}  // namespace detray::concepts
