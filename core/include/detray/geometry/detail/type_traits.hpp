/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <type_traits>

namespace detray::detail {

template <typename = void>
struct is_planar : public std::false_type {};

template <typename T>
inline constexpr bool is_planar_v{is_planar<T>::value};

template <typename = void>
struct is_cylindrical : public std::false_type {};

template <typename T>
inline constexpr bool is_cylindrical_v{is_cylindrical<T>::value};

template <typename = void>
struct is_linear : public std::false_type {};

template <typename T>
inline constexpr bool is_linear_v{is_linear<T>::value};

template <typename = void>
struct is_spherical : public std::false_type {};

template <typename T>
inline constexpr bool is_spherical_v{is_spherical<T>::value};

}  // namespace detray::detail
