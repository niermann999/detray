/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/navigation/detail/ray.hpp"

// Detray test include(s)
#include "detray/test/utils/simulation/event_generator/random_numbers.hpp"

// System include(s)
#include <type_traits>

namespace detray {

namespace detail {

/// Wrapper for random number generatrion
template <concepts::algebra algebra_t,
          typename generator_t = random_numbers<dvalue<algebra_t>>>
DETRAY_HOST_DEVICE dvector3D<algebra_t> random_direction(
    generator_t &rand_gen, const dvalue<algebra_t> min = 0.f,
    const dvalue<algebra_t> max = 1.f) {

    using vector3_t = dvector3D<algebra_t>;

    return detray::vector::normalize(vector3_t{
        rand_gen({min, max}), rand_gen({min, max}), rand_gen({min, max})});
};

}  // namespace detail

/// @returns a ray that is randomly scattered with respect to a surface
template <concepts::algebra algebra_t, typename generator_t>
DETRAY_HOST_DEVICE void random_scattering(detail::ray<algebra_t> &ray,
                                          const dvector3D<algebra_t> &sf_normal,
                                          const dscalar<algebra_t> path,
                                          generator_t &rand_gen,
                                          const dvalue<algebra_t> min = 0.f,
                                          const dvalue<algebra_t> max = 1.f) {

    using scalar_t = dscalar<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    const vector3_t rand_dir =
        detail::random_direction<algebra_t>(rand_gen, min, max);
    const scalar_t sign =
        math::copysign(1.f, -vector::dot(rand_dir, sf_normal));

    ray.set_pos(ray.pos() + path * ray.dir());
    ray.set_dir(sign * vector::normalize(rand_dir));
};

}  // namespace detray
