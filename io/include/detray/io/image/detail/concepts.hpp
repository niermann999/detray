/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <concepts>
#include <vector>

namespace detray::io::concepts {

template <typename I>
concept image = requires(const I im) {

    typename I::color;
    typename I::aspect_ratio;

    { im.width() }
    ->std::integral;

    { im.height() }
    ->std::integral;

    { im.n_pixels() }
    ->std::integral;

    { im.pixel_data() }
    ->std::same_as<const std::vector<typename I::color>&>;
};

}  // namespace detray::io::concepts
