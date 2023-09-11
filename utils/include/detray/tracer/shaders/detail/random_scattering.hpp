/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/utils/random_numbers.hpp"

// System include(s)
#include <type_traits>

namespace detray {

namespace detail {

/// Wrapper for random number generatrion
template <typename T = detray::scalar, template <typename> class algebra_t,
          typename generator_t = random_numbers<>>
DETRAY_HOST_DEVICE vector3D random_direction(generator_t &rand_gen,
                                             const T min = 0.f,
                                             const T max = 1.f) {

    using vector3D = dvector3D<algebra_t<T>>;

    return detray::vector::normalize(
        vector3D{rand_gen(min, max), rand_gen(min, max), rand_gen(min, max)});
};

}  // namespace detail

/// @returns a ray that is randomly scattered with respect to a surface
template <typename T = detray::scalar, template <typename> class algebra_t,
          typename generator_t>
DETRAY_HOST_DEVICE void random_scattering(
    detail::ray<dtransform3D<algebra_t<T>>> &ray,
    const dscalar<algebra_t<T>> path, generator_t &rand_gen, const T min = 0.f,
    const T max = 1.f) {

    using scalar_t = dscalar<algebra_t<T>>;
    using vector3D = dvector3D<algebra_t<T>>;

    const vector3D rand_dir =
        detail::random_direction<T, algebra_t>(rand_gen, min, max);
    const scalar_t sign =
        math_ns::copysign(1.f, -vector::dot(rand_dir, ray.dir()));

    // std::cout << "In mat shader: " << path << std::endl;
    // std::cout << ray << std::endl;

    ray.set_pos(ray.pos() + path * ray.dir());
    ray.set_dir(sign * rand_dir);
};

}  // namespace detray
