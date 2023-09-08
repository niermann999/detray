/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/io/image/ppm_writer.hpp"
#include "detray/tracer/renderer/raw_image.hpp"

using namespace detray;

namespace {

/// Simple color gradient
template <typename color_depth>
inline void write_test_image(raw_image<color_depth> &im) {
    // Iterate through pixel matrix
    for (float i_y{static_cast<float>(im.height() - 1u)}; i_y >= 0.f;
         i_y -= 1.f) {
        for (float i_x{0.f}; i_x < static_cast<float>(im.width()); i_x += 1.f) {
            const float r{i_x / static_cast<float>(im.width())};
            const float g{i_y / static_cast<float>(im.height())};
            const float b{0.2f};

            const texture::color<color_depth> c_grad{
                static_cast<color_depth>(255.99f * r),
                static_cast<color_depth>(255.99f * g),
                static_cast<color_depth>(255.99f * b), 0u};

            im.set_pixel(static_cast<uint>(i_x), static_cast<uint>(i_y),
                         c_grad);
        }
    }
}

}  // namespace

int main() {

    io::ppm_writer<unsigned int> ppm{};

    // write a test image
    raw_image<unsigned int> image{500u};
    write_test_image(image);
    ppm.write(image, "test");

    return EXIT_SUCCESS;
}
