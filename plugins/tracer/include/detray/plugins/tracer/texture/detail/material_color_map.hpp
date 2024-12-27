/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/predefined_materials.hpp"

// Detray tracer include(s)
#include "detray/plugins/tracer/definitions/colors.hpp"

namespace detray::texture::detail {

/// @brief holds rgb and alpha values for color shading of different materials
template <typename color_depth = std::uint8_t, typename scalar_t = float>
DETRAY_HOST_DEVICE inline constexpr detray::texture::color<color_depth>
material_color_map(const material<scalar_t> &mat) {
    // color based on material
    if (mat == detray::beryllium<scalar_t>{} ||
        mat == detray::beryllium_tml<scalar_t>{}) {
        return detray::texture::grey<color_depth>;
    } else if (mat == detray::aluminium<scalar_t>{}) {
        return detray::texture::light_grey<color_depth>;
    } else if (mat == detray::tungsten<scalar_t>{}) {
        return detray::texture::dim_grey<color_depth>;
    } else if (mat == detray::gold<scalar_t>{}) {
        return detray::texture::golden_rod<color_depth>;
    } else if (mat == detray::silicon<scalar_t>{} ||
               mat == detray::silicon_tml<scalar_t>{}) {
        return detray::texture::dark_grey<color_depth>;
    } else {
        // default for unknown material
        return detray::texture::dark_red<color_depth>;
    }
}

}  // namespace detray::texture::detail
