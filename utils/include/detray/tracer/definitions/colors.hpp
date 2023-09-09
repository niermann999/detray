/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/tracer/texture/color.hpp"

namespace detray::texture {

/// Color definitions
/// @{

// Macro for declaring rgb colors
#define DETRAY_DECLARE_COLOR(COLOR_NAME, R, G, B)                            \
    template <typename color_depth = uint8_t>                                \
    inline constexpr auto COLOR_NAME = texture::color<color_depth>{R, G, B}; \
    template <>                                                              \
    inline constexpr auto COLOR_NAME<float> =                                \
        texture::color<float>{R / 256.f, G / 256.f, B / 256.f, 1.f};         \
    template <>                                                              \
    inline constexpr auto COLOR_NAME<double> =                               \
        texture::color<double>{R / 256., G / 256., B / 256., 1.};

// https://www.w3schools.com/colors/color_tryit.asp?hex=F0F8FF
// basic
DETRAY_DECLARE_COLOR(black, 0, 0, 0);
DETRAY_DECLARE_COLOR(white, 255, 255, 255);
DETRAY_DECLARE_COLOR(red, 255, 0, 0);
DETRAY_DECLARE_COLOR(lime, 0, 255, 0);
DETRAY_DECLARE_COLOR(blue, 0, 0, 255);
DETRAY_DECLARE_COLOR(yellow, 255, 255, 0);

DETRAY_DECLARE_COLOR(alice_blue, 240, 248, 255)
DETRAY_DECLARE_COLOR(antique_white, 250, 235, 215)
DETRAY_DECLARE_COLOR(aqua, 0, 255, 255)
DETRAY_DECLARE_COLOR(aquamarine, 127, 255, 212)
DETRAY_DECLARE_COLOR(azure, 240, 255, 255)
DETRAY_DECLARE_COLOR(beige, 245, 245, 220)
DETRAY_DECLARE_COLOR(bisque, 255, 228, 196)
DETRAY_DECLARE_COLOR(blanched_almond, 255, 235, 205)
DETRAY_DECLARE_COLOR(blue_violet, 138, 43, 226)
DETRAY_DECLARE_COLOR(brown, 165, 42, 42)
DETRAY_DECLARE_COLOR(burly_wood, 222, 184, 135)
DETRAY_DECLARE_COLOR(cadet_blue, 95, 158, 160)
DETRAY_DECLARE_COLOR(chartreuse, 127, 255, 0)
DETRAY_DECLARE_COLOR(chocolate, 210, 105, 30)
DETRAY_DECLARE_COLOR(coral, 255, 127, 80)
DETRAY_DECLARE_COLOR(cornflower_blue, 100, 149, 237)
DETRAY_DECLARE_COLOR(cornsilk, 255, 248, 220)
DETRAY_DECLARE_COLOR(crimson, 220, 20, 60)
DETRAY_DECLARE_COLOR(cyan, 0, 255, 255)
DETRAY_DECLARE_COLOR(dark_blue, 0, 0, 139)
DETRAY_DECLARE_COLOR(dark_cyan, 0, 139, 139)
DETRAY_DECLARE_COLOR(dark_golden_rod, 184, 134, 11)
DETRAY_DECLARE_COLOR(dark_grey, 169, 169, 169)
DETRAY_DECLARE_COLOR(dark_green, 0, 100, 0)
DETRAY_DECLARE_COLOR(dark_khaki, 189, 183, 107)
DETRAY_DECLARE_COLOR(dark_magenta, 139, 0, 139)
DETRAY_DECLARE_COLOR(dark_olive_green, 85, 107, 47)
DETRAY_DECLARE_COLOR(dark_orange, 255, 140, 0)
DETRAY_DECLARE_COLOR(dark_orchid, 153, 50, 204)
DETRAY_DECLARE_COLOR(dark_red, 139, 0, 0)
DETRAY_DECLARE_COLOR(dark_salmon, 233, 150, 122)
DETRAY_DECLARE_COLOR(dark_sea_green, 143, 188, 143)
DETRAY_DECLARE_COLOR(dark_slate_blue, 72, 61, 139)
DETRAY_DECLARE_COLOR(dark_slate_grey, 47, 79, 79)
DETRAY_DECLARE_COLOR(dark_turquoise, 0, 206, 209)
DETRAY_DECLARE_COLOR(dark_violet, 148, 0, 211)
DETRAY_DECLARE_COLOR(deep_pink, 255, 20, 147)
DETRAY_DECLARE_COLOR(deep_sky_blue, 0, 191, 255)
DETRAY_DECLARE_COLOR(dim_grey, 105, 105, 105)
DETRAY_DECLARE_COLOR(dodger_blue, 30, 144, 255)
DETRAY_DECLARE_COLOR(fire_brick, 178, 34, 34)
DETRAY_DECLARE_COLOR(floral_white, 255, 250, 240)
DETRAY_DECLARE_COLOR(forest_green, 34, 139, 34)
DETRAY_DECLARE_COLOR(fuchsia, 255, 0, 255)
DETRAY_DECLARE_COLOR(gainsboro, 220, 220, 220)
DETRAY_DECLARE_COLOR(ghost_white, 248, 248, 255)
DETRAY_DECLARE_COLOR(gold, 255, 215, 0)
DETRAY_DECLARE_COLOR(golden_rod, 218, 165, 32)
DETRAY_DECLARE_COLOR(grey, 128, 128, 128)
DETRAY_DECLARE_COLOR(green, 0, 128, 0)
DETRAY_DECLARE_COLOR(green_yellow, 173, 255, 47)
DETRAY_DECLARE_COLOR(honey_dew, 240, 255, 240)
DETRAY_DECLARE_COLOR(hot_pink, 255, 105, 180)
DETRAY_DECLARE_COLOR(indian_red, 205, 92, 92)
DETRAY_DECLARE_COLOR(indigo, 75, 0, 130)
DETRAY_DECLARE_COLOR(ivory, 255, 255, 240)
DETRAY_DECLARE_COLOR(khaki, 240, 230, 140)
DETRAY_DECLARE_COLOR(lavender, 230, 230, 250)
DETRAY_DECLARE_COLOR(lavender_blush, 255, 240, 245)
DETRAY_DECLARE_COLOR(lawn_green, 124, 252, 0)
DETRAY_DECLARE_COLOR(lemon_chiffon, 255, 250, 205)
DETRAY_DECLARE_COLOR(light_blue, 173, 216, 230)
DETRAY_DECLARE_COLOR(light_coral, 240, 128, 128)
DETRAY_DECLARE_COLOR(light_cyan, 224, 255, 255)
DETRAY_DECLARE_COLOR(light_golden_rod_yellow, 250, 250, 210)
DETRAY_DECLARE_COLOR(light_grey, 211, 211, 211)
DETRAY_DECLARE_COLOR(light_green, 144, 238, 144)
DETRAY_DECLARE_COLOR(light_pink, 255, 182, 193)
DETRAY_DECLARE_COLOR(light_salmon, 255, 160, 122)
DETRAY_DECLARE_COLOR(light_sea_green, 32, 178, 170)
DETRAY_DECLARE_COLOR(light_sky_blue, 135, 206, 250)
DETRAY_DECLARE_COLOR(light_slate_grey, 119, 136, 153)
DETRAY_DECLARE_COLOR(light_steel_blue, 176, 196, 222)
DETRAY_DECLARE_COLOR(light_yellow, 255, 255, 224)
DETRAY_DECLARE_COLOR(lime_green, 50, 205, 50)
DETRAY_DECLARE_COLOR(linen, 250, 240, 230)
DETRAY_DECLARE_COLOR(magenta, 255, 0, 255)
DETRAY_DECLARE_COLOR(maroon, 128, 0, 0)
DETRAY_DECLARE_COLOR(medium_aqua_marine, 102, 205, 170)
DETRAY_DECLARE_COLOR(medium_blue, 0, 0, 205)
DETRAY_DECLARE_COLOR(medium_orchid, 186, 85, 211)
DETRAY_DECLARE_COLOR(medium_purple, 147, 112, 219)
DETRAY_DECLARE_COLOR(medium_sea_green, 60, 179, 113)
DETRAY_DECLARE_COLOR(medium_slate_blue, 123, 104, 238)
DETRAY_DECLARE_COLOR(medium_spring_green, 0, 250, 154)
DETRAY_DECLARE_COLOR(medium_turquoise, 72, 209, 204)
DETRAY_DECLARE_COLOR(medium_violet_red, 199, 21, 133)
DETRAY_DECLARE_COLOR(midnight_blue, 25, 25, 112)
DETRAY_DECLARE_COLOR(mint_cream, 245, 255, 250)
DETRAY_DECLARE_COLOR(misty_rose, 255, 228, 225)
DETRAY_DECLARE_COLOR(moccasin, 255, 228, 181)
DETRAY_DECLARE_COLOR(navajo_white, 255, 222, 173)
DETRAY_DECLARE_COLOR(navy, 0, 0, 128)
DETRAY_DECLARE_COLOR(old_lace, 253, 245, 230)
DETRAY_DECLARE_COLOR(olive, 128, 128, 0)
DETRAY_DECLARE_COLOR(olive_drab, 107, 142, 35)
DETRAY_DECLARE_COLOR(orange, 255, 165, 0)
DETRAY_DECLARE_COLOR(orange_red, 255, 69, 0)
DETRAY_DECLARE_COLOR(orchid, 218, 112, 214)
DETRAY_DECLARE_COLOR(pale_golden_rod, 238, 232, 170)
DETRAY_DECLARE_COLOR(pale_green, 152, 251, 152)
DETRAY_DECLARE_COLOR(pale_turquoise, 175, 238, 238)
DETRAY_DECLARE_COLOR(pale_violet_red, 219, 112, 147)
DETRAY_DECLARE_COLOR(papaya_whip, 255, 239, 213)
DETRAY_DECLARE_COLOR(peach_puff, 255, 218, 185)
DETRAY_DECLARE_COLOR(peru, 205, 133, 63)
DETRAY_DECLARE_COLOR(pink, 255, 192, 203)
DETRAY_DECLARE_COLOR(plum, 221, 160, 221)
DETRAY_DECLARE_COLOR(powder_blue, 176, 224, 230)
DETRAY_DECLARE_COLOR(purple, 128, 0, 128)
DETRAY_DECLARE_COLOR(rebecca_purple, 102, 51, 153)
DETRAY_DECLARE_COLOR(rosy_brown, 188, 143, 143)
DETRAY_DECLARE_COLOR(royal_blue, 65, 105, 225)
DETRAY_DECLARE_COLOR(saddle_brown, 139, 69, 19)
DETRAY_DECLARE_COLOR(salmon, 250, 128, 114)
DETRAY_DECLARE_COLOR(sandy_brown, 244, 164, 96)
DETRAY_DECLARE_COLOR(sea_green, 46, 139, 87)
DETRAY_DECLARE_COLOR(sea_shell, 255, 245, 238)
DETRAY_DECLARE_COLOR(sienna, 60, 82, 45)
DETRAY_DECLARE_COLOR(silver, 192, 192, 192)
DETRAY_DECLARE_COLOR(sky_blue, 135, 206, 235)
DETRAY_DECLARE_COLOR(slate_blue, 106, 90, 205)
DETRAY_DECLARE_COLOR(slate_grey, 112, 128, 144)
DETRAY_DECLARE_COLOR(snow, 255, 250, 250)
DETRAY_DECLARE_COLOR(spring_green, 0, 255, 127)
DETRAY_DECLARE_COLOR(steel_blue, 70, 130, 180)
DETRAY_DECLARE_COLOR(tan, 210, 180, 140)
DETRAY_DECLARE_COLOR(teal, 0, 128, 128)
DETRAY_DECLARE_COLOR(thistle, 216, 191, 216)
DETRAY_DECLARE_COLOR(tomato, 255, 99, 71)
DETRAY_DECLARE_COLOR(turquoise, 64, 224, 208)
DETRAY_DECLARE_COLOR(violet, 238, 130, 238)
DETRAY_DECLARE_COLOR(wheat, 245, 222, 179)
DETRAY_DECLARE_COLOR(white_smoke, 245, 245, 245)
DETRAY_DECLARE_COLOR(yellow_green, 154, 205, 50)
/// @}

}  // namespace detray::texture
