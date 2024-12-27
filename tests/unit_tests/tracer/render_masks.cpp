/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algebra include(s).
//#include "detray/plugins/algebra/array_definitions.hpp"
#include "detray/plugins/algebra/vc_aos_definitions.hpp"
#include "detray/plugins/algebra/vc_soa_definitions.hpp"

// Project include(s).
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/intersection/soa/ray_cylinder_intersector.hpp"
#include "detray/navigation/intersection/soa/ray_cylinder_portal_intersector.hpp"
#include "detray/navigation/intersection/soa/ray_line_intersector.hpp"
#include "detray/navigation/intersection/soa/ray_plane_intersector.hpp"
//#include "detray/navigation/intersection/soa/sphere_intersector.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
//#include "detray/masks/sphere2D.hpp"
#include "detray/materials/predefined_materials.hpp"

// Detray tracer include(s)
#include "detray/plugins/tracer/renderer/camera.hpp"
#include "detray/plugins/tracer/renderer/geometry.hpp"
#include "detray/plugins/tracer/renderer/intersector.hpp"
#include "detray/plugins/tracer/renderer/pipeline.hpp"
#include "detray/plugins/tracer/renderer/raw_image.hpp"
#include "detray/plugins/tracer/shaders/background.hpp"
#include "detray/plugins/tracer/shaders/material.hpp"
#include "detray/plugins/tracer/texture/color.hpp"
#include "detray/plugins/tracer/texture/pixel.hpp"

// Detray io include(s)
#include "detray/io/image/ppm_writer.hpp"

// Detray test include(s)
#include "detray/test/utils/simulation/event_generator/random_numbers.hpp"

// System include(s)
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <tuple>

using namespace detray;

namespace {

/// Render a shape
template <concepts::algebra algebra_t, std::size_t SAMPLES,
          typename color_depth, typename aspect_ratio, typename mask_t,
          typename material_t,
          class im_background_t =
              inf_plane<gradient_background<detray::vc_aos<dvalue<algebra_t>>>>>
inline void render_mask(raw_image<color_depth, aspect_ratio> &im,
                        std::vector<mask_t> &&mask,
                        const std::vector<dtransform3D<algebra_t>> &trf,
                        const std::vector<material_t> &mat) {

    using T = dvalue<algebra_t>;
    using scalar_t = dscalar<algebra_t>;
    using generator_t = detail::random_numbers<T>;
    using geometry_t = simple_geometry<algebra_t, mask_t>;

    // Rendering steps
    using intersector_t =
        intersector<geometry_t, detray::vc_aos<T>, SAMPLES, color_depth>;
    using backgr_shader_t = background_shader<im_background_t>;
    using mat_shader_t = material_shader<algebra_t, generator_t>;
    // The rendering pipeline: The intersector finds the shape intersections
    using pipeline_t =
        rendering_pipeline<intersector_t, backgr_shader_t, mat_shader_t>;

    // Random numbers for ray generation
    auto rand_gen = generator_t{};
    const T viewport_height = 2.0f;
    const dpoint3D<detray::vc_aos<T>> origin{0.0f, 0.0f, 0.0f};

    camera<detray::vc_aos<T>, aspect_ratio> cam(viewport_height, origin);

    // For the single shape render, the scene is actually encoded directly in
    // the single shape intersector
    geometry_t geo{std::move(trf), std::move(mask), std::move(mat)};

// Iterate through pixel matrix
#pragma omp parallel for collapse(2)
    for (std::size_t i_y = 0u; i_y < im.height(); ++i_y) {
        for (std::size_t i_x = 0u; i_x < im.width(); ++i_x) {

            // Ray to render the pixel at (i_x, i_y)
            auto rays =
                cam.template generate_rays<SAMPLES>(i_x, i_y, rand_gen, im);

            // Strap the global geometry state and the image together
            scene_handle<geometry_t, color_depth> scene{geo, im};

            // Finds the intersections between the ray and the geometry
            typename intersector_t::state intrs{rays, 1000};
            typename mat_shader_t::state mat_state{rand_gen};

            auto pipeline_state = std::tie(intrs, mat_state);

            // Run while at leat one ray is hitting a surface
            do {
                pipeline_t{}(pipeline_state, scene);
            } while (!intrs.is_finished());

            // Average the pixels for this ray boundle (antialiasing)
            texture::color<T> col{};
            for (const auto &c : intrs.m_colors) {
                col += static_cast<texture::color<T>>(c);
            }
            // Normalize the pixel color
            constexpr T scalor{1.f / SAMPLES};
            col *= scalor;

            im.set_pixel(i_x, i_y,
                         static_cast<texture::color<color_depth>>(col));
        }
    }
}

}  // namespace

using T = float;

/// Linear algebra implementation using SoA memory layout
using algebra_v = detray::vc_soa<T>;

/// Linear algebra implementation using AoS memory layout
using algebra_s = detray::vc_aos<T>;

int main() {

    using image_t = raw_image<std::uint8_t>;

    io::ppm_writer<image_t> ppm{};

    image_t image{500u};

    //
    // Render single shape
    //

    using vector3D_s = dvector3D<algebra_s>;
    using vector3D_v = dvector3D<algebra_v>;

    constexpr std::size_t simd_size{dscalar<algebra_v>::size()};
    constexpr std::size_t n_samples{100ul};

    // Affine transform matrix to place the shapes

    // SoA
    const T im_height{static_cast<T>(image.height())};
    const T im_width{static_cast<T>(image.width())};
    vector3D_v x_v{1.0f, 0.0f, 0.0f};
    vector3D_v z_v{0.0f, 0.0f, 1.f};
    vector3D_v t_v{30.0f, -20.0f, 0.0f};
    t_v[0] = t_v[0].Random();
    t_v[0] = 0.1f * (im_width * t_v[0] - 0.5f * im_width);
    t_v[1] = t_v[1].Random();
    t_v[1] = 0.1f * (im_height * t_v[1] - 0.5f * im_height);
    t_v[2] = -120.1f * math::abs(t_v[1].Random());

    std::vector<dtransform3D<algebra_v>> trfs_v;
    trfs_v.emplace_back(t_v, z_v, x_v);

    // AoS
    std::vector<dtransform3D<algebra_s>> trfs_s;
    trfs_s.reserve(simd_size);
    for (std::size_t i = 0; i < simd_size; ++i) {
        vector3D_s x_s{x_v[0][i], x_v[1][i], x_v[2][i]};
        vector3D_s z_s{z_v[0][i], z_v[1][i], z_v[2][i]};
        vector3D_s t_s{t_v[0][i], t_v[1][i], t_v[2][i]};

        trfs_s.emplace_back(t_s, z_s, x_s);
    }

    // Different materials per surface
    std::vector<material<T>> mat{beryllium<T>{}, aluminium<T>{}, gold<T>{},
                                 silicon<T>{},   tungsten<T>{},  gold<T>{},
                                 aluminium<T>{}, silicon<T>{}};

    // render a rectangle mask

    // AoS
    mask<rectangle2D, algebra_s> rect2_s{0u, 0.01f * im_width,
                                         0.01f * im_height};
    std::vector<mask<rectangle2D, algebra_s>> rect2_vec(simd_size, rect2_s);

    auto start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_s, n_samples>(image, std::move(rect2_vec), trfs_s, mat);
    auto end = std::chrono::high_resolution_clock::now();

    auto time_aos{
        static_cast<double>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count()) /
        1000'000.};
    std::cout << "\nRectangle AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "rectangle_AoS");

    // SoA
    /*std::vector<mask<rectangle2D, algebra_v>> rect2_v;
    rect2_v.emplace_back(0u, 0.01f * im_width, 0.01f * im_height);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_v, n_samples>(image, std::move(rect2_v), trfs_v, mat);
    end = std::chrono::high_resolution_clock::now();

    auto time_soa{
        static_cast<double>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count()) /
        1000'000.};
    std::cout << "Rectangle SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "rectangle_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;*/

    // render a trapezoid mask

    // AoS
    const mask<trapezoid2D, algebra_s> trpz2_s{0u, 10.f, 30.f, 20.f,
                                               1.f / 40.f};
    std::vector<mask<trapezoid2D, algebra_s>> trpz2_vec(simd_size, trpz2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_s, n_samples>(image, std::move(trpz2_vec), trfs_s, mat);

    end = std::chrono::high_resolution_clock::now();

    time_aos =
        static_cast<double>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count()) /
        1000'000.;
    std::cout << "\nTrapezoid AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "trapezoid_AoS");

    // SoA
    /*std::vector<mask<trapezoid2D<soa::plane_intersector>, std::uint_least16_t,
                     algebra_v>>
        trap2_v;
    trap2_v.emplace_back(0u, 10.f, 30.f, 20.f, 1.f / 40.f);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_v, n_samples>(image, std::move(trap2_v), trfs_v,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Trapezoid SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "trapezoid_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;*/

    // render a ring mask

    // AoS
    const mask<ring2D, algebra_s> ring2_s{0u, 12.f, 20.f};
    std::vector<mask<ring2D, algebra_s>> ring2_vec(simd_size, ring2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_s, n_samples>(image, std::move(ring2_vec), trfs_s, mat);
    end = std::chrono::high_resolution_clock::now();

    time_aos =
        static_cast<double>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count()) /
        1000'000.;
    std::cout << "\nRing AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "ring_AoS");

    // SoA
    /*std::vector<mask<ring2D<soa::plane_intersector>, std::uint_least16_t,
                     algebra_v>>
        ring2_v;
    ring2_v.emplace_back(0u, 12.f, 20.f);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_v, n_samples>(image, std::move(ring2_v), trfs_v,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Ring SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "ring_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;*/

    // render an annulus mask

    // AoS
    const mask<annulus2D, algebra_s> ann2_s{0u,       5.f,  13.0f, 0.74195f,
                                            1.33970f, -2.f, 2.f,   0.f};
    std::vector<mask<annulus2D, algebra_s>> ann2_vec(simd_size, ann2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_s, n_samples>(image, std::move(ann2_vec), trfs_s, mat);
    end = std::chrono::high_resolution_clock::now();

    time_aos =
        static_cast<double>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count()) /
        1000'000.;
    std::cout << "\nAnnulus AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "annulus_AoS");

    // SoA
    /*std::vector<mask<annulus2D<soa::plane_intersector>, std::uint_least16_t,
                     algebra_v>>
        ann2_v;
    ann2_v.emplace_back(0u, 5.f, 13.0f, 0.74195f, 1.33970f, -2.f, 2.f, 0.f);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_v, n_samples>(image, std::move(ann2_v), trfs_v,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Annulus SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "annulus_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;

    // render a spherical mask

    // AoS
    const mask<sphere2D<>, std::uint_least16_t, algebra_s> sph2_s{0u,
                                                                          10.f};
    std::vector<mask<sphere2D<>>> sph2_vec(simd_size, sph2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_s, n_samples>(image, std::move(sph2_vec),
                                              trfs_s, mat);
    end = std::chrono::high_resolution_clock::now();

    time_aos = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "\nSphere AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "sphere_AoS");

    // SoA
    std::vector<mask<sphere2D<soa::sphere_intersector>, std::uint_least16_t,
                     algebra_v>>
        sph2_v;
    sph2_v.emplace_back(0u, 10.f);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_v, n_samples>(image, std::move(sph2_v), trfs_v,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Sphere SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "sphere_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;*/

    // render a line mask

    // AoS
    const mask<line<true>, algebra_s> ln2_s{0u, 10.f,
                                            std::numeric_limits<T>::max()};
    std::vector<mask<line<true>, algebra_s>> ln2_vec(simd_size, ln2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_s, n_samples>(image, std::move(ln2_vec), trfs_s, mat);
    end = std::chrono::high_resolution_clock::now();

    time_aos =
        static_cast<double>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count()) /
        1000'000.;
    std::cout << "\nLine AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "line_AoS");

    // SoA
    /*std::vector<mask<line<true, soa::line_intersector>, std::uint_least16_t,
                     algebra_v>>
        ln2_v;
    ln2_v.emplace_back(0u, 10.f, std::numeric_limits<scalar>::max());

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_v, n_samples>(image, std::move(ln2_v), trfs_v,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Line SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "line_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;*/

    // render a cylinder mask

    // AoS
    const mask<cylinder2D, algebra_s> cyl2_s{0u, 0.5f * im_height,
                                             0.5f * im_width, 0.7f * im_width};
    std::vector<mask<cylinder2D, algebra_s>> cyl2_vec(simd_size, cyl2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_s, n_samples>(image, std::move(cyl2_vec), trfs_s, mat);
    end = std::chrono::high_resolution_clock::now();

    time_aos =
        static_cast<double>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count()) /
        1000'000.;
    std::cout << "\nCylinder AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "cylinder_AoS");

    // SoA
    /*std::vector<mask<cylinder2D<false, soa::cylinder_intersector>,
                     std::uint_least16_t, algebra_v>>
        cyl2_v;
    cyl2_v.emplace_back(0u, 0.5f * im_height, 0.5f * im_width, 0.7f * im_width);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_v, n_samples>(image, std::move(cyl2_v), trfs_v,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Cylinder SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "cylinder_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;*/

    // render a portal cylinder mask

    // AoS
    const mask<concentric_cylinder2D, algebra_s> pt_cyl2_s{
        0u, 0.5f * im_height, 0.5f * im_width, 0.7f * im_width};
    std::vector<mask<concentric_cylinder2D, algebra_s>> pt_cyl2_vec(simd_size,
                                                                    pt_cyl2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_s, n_samples>(image, std::move(pt_cyl2_vec), trfs_s,
                                      mat);
    end = std::chrono::high_resolution_clock::now();

    time_aos =
        static_cast<double>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                .count()) /
        1000'000.;
    std::cout << "\nPortal Cylinder AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "portal_cylinder_AoS");

    // SoA
    /*std::vector<mask<cylinder2D<false, soa::cylinder_portal_intersector>,
                     std::uint_least16_t, algebra_v>>
        pt_cyl2_v;
    pt_cyl2_v.emplace_back(0u, 0.5f * im_height, 0.5f * im_width,
                           0.7f * im_width);

    start = std::chrono::high_resolution_clock::now();
    render_mask<algebra_v, n_samples>(image, std::move(pt_cyl2_v),
                                              trfs_v, mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Portal Cylinder SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "portal_cylinder_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl << std::endl;*/

    return EXIT_SUCCESS;
}
