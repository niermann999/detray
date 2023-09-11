/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algebra include(s).
#include "detray/plugins/algebra/vc_soa_definitions.hpp"

// Project include(s).
#include "detray/intersection/cylinder_portal_intersector.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/soa/cylinder_intersector.hpp"
#include "detray/intersection/soa/cylinder_portal_intersector.hpp"
#include "detray/intersection/soa/line_intersector.hpp"
#include "detray/intersection/soa/plane_intersector.hpp"
#include "detray/intersection/soa/sphere_intersector.hpp"
#include "detray/io/image/ppm_writer.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/sphere2D.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/tracer/renderer/camera.hpp"
#include "detray/tracer/renderer/pipeline.hpp"
#include "detray/tracer/renderer/raw_image.hpp"
#include "detray/tracer/shaders/background.hpp"
#include "detray/tracer/shaders/material.hpp"
#include "detray/tracer/texture/color.hpp"
#include "detray/tracer/texture/pixel.hpp"
#include "detray/utils/random_numbers.hpp"

// System include(s)
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <tuple>

using namespace detray;

namespace {

/// Render a shape
template <typename T, template <typename> class algebra_t, std::size_t SAMPLES,
          typename color_depth, typename aspect_ratio, typename mask_t,
          typename material_t,
          class im_background_t = gradient_background<T, ALGEBRA_PLUGIN>>
inline void render_mask(raw_image<color_depth, aspect_ratio> &im,
                        std::vector<mask_t> &&mask,
                        const std::vector<dtransform3D<algebra_t<T>>> &trf,
                        const std::vector<material_t> &mat) {
    using generator_t = detail::random_numbers<>;

    // Rendering steps
    using intersector_t =
        single_shape<T, algebra_t, SAMPLES, mask_t, material_t>;
    // using backgr_shader_t = background_shader<inf_plane<im_background_t>>;
    using backgr_shader_t = background_shader<im_background_t>;
    using mat_shader_t = material_shader<T, ALGEBRA_PLUGIN, generator_t>;
    // The rendering pipeline: The intersector finds the shape intersections
    using pipeline_t =
        rendering_pipeline<intersector_t, backgr_shader_t, mat_shader_t>;

    // Random numbers for ray generation
    auto rand_gen = generator_t{};
    const T viewport_height = 2.0f;
    const dpoint3D<ALGEBRA_PLUGIN<T>> origin{0.0f, 0.0f, 0.0f};

    camera<T, ALGEBRA_PLUGIN, aspect_ratio> cam(viewport_height, origin);

    // For the single shape render, the scene is actually encoded directly in
    // the single shape intersector
    typename intersector_t::global_state geo{std::move(trf), std::move(mask),
                                             std::move(mat)};

// Iterate through pixel matrix
#pragma omp parallel for collapse(2)
    for (std::size_t i_y = 0u; i_y < im.height(); ++i_y) {
        for (std::size_t i_x = 0u; i_x < im.width(); ++i_x) {

            // Ray to render the pixel at (i_x, i_y)
            auto rays = cam.template get_rays<SAMPLES>(i_x, i_y, rand_gen, im);

            // Strap the global geometry state and the thread-local ray together
            scene_handle::state scene{geo, im, rays};

            // Finds the intersections between the ray and the geometry
            typename intersector_t::state intrs{};
            typename mat_shader_t::state mat_state{rand_gen};

            actor::state empty{};  // Dummy for actors that don't define a state
            auto pipeline_state = std::tie(empty, intrs, mat_state);

            // Run while at leat one ray is hitting a surface
            std::size_t n_reflections{0};  // < prevent infinite reflections
            do {
                pipeline_t{}(pipeline_state, scene);
                ++n_reflections;

                /*if (intrs.has_hit()) {
                    std::cout << "REFLECTION: " << n_reflections << std::endl;
                    for (const auto&[i, r] : detray::views::enumerate(rays)) {
                        std::cout << intrs.m_intersection[i] << std::endl;
                        std::cout << intrs.path_to_closest(i) << std::endl;
                        //std::cout << intrs.m_intersection[i].path <<
                std::endl; std::cout << r << std::endl;
                    }
                }*/

            } while (intrs.has_hit() and n_reflections < 100);

            // Average the pixels for this ray boundle (antialiasing)
            texture::pixel<std::size_t, T> pix{{i_x, i_y}};
            for (const auto &c : scene.m_colors) {
                pix += c;
            }
            // Normalize the pixel color
            constexpr T scalor{1.f / SAMPLES};
            pix *= scalor;

            im.set_pixel(pix);
        }
    }
}

}  // namespace

/// Linear algebra implementation using SoA memory layout
template <typename T>
using algebra_v = detray::vc_soa<T>;

/// Linear algebra implementation using AoS memory layout
template <typename T>
using algebra_s = detray::cmath<T>;

int main() {

    using color_depth = std::uint8_t;

    io::ppm_writer<color_depth> ppm{};

    raw_image<color_depth> image{500u};

    //
    // Render single shape
    //

    using vector3D_s = dvector3D<algebra_s<scalar>>;
    using vector3D_v = dvector3D<algebra_v<scalar>>;

    constexpr std::size_t simd_size{dscalar<algebra_v<scalar>>::size()};
    constexpr std::size_t n_samples{100ul};

    // Affine transform matrix to place the shapes

    // SoA
    vector3D_v x_v{1.0f, 0.0f, 0.0f};
    vector3D_v z_v{0.0f, 0.0f, 1.f};
    vector3D_v t_v{30.0f, -20.0f, 0.0f};
    t_v[0] = t_v[0].Random();
    t_v[0] = 0.1f * (image.width() * t_v[0] - 0.5f * image.width());
    t_v[1] = t_v[1].Random();
    t_v[1] = 0.1f * (image.height() * t_v[1] - 0.5f * image.height());
    t_v[2] = -120.1f * math_ns::abs(t_v[1].Random());

    std::vector<dtransform3D<algebra_v<scalar>>> trfs_v;
    trfs_v.emplace_back(t_v, z_v, x_v);

    // AoS
    std::vector<dtransform3D<algebra_s<scalar>>> trfs_s;
    trfs_s.reserve(simd_size);
    for (std::size_t i = 0; i < simd_size; ++i) {
        vector3D_s x_s{x_v[0][i], x_v[1][i], x_v[2][i]};
        vector3D_s z_s{z_v[0][i], z_v[1][i], z_v[2][i]};
        vector3D_s t_s{t_v[0][i], t_v[1][i], t_v[2][i]};

        trfs_s.emplace_back(t_s, z_s, x_s);
    }

    // Different materials per surface
    std::vector<material<scalar>> mat{beryllium<scalar>{}, aluminium<scalar>{},
                                      gold<scalar>{},      silicon<scalar>{},
                                      tungsten<scalar>{},  gold<scalar>{},
                                      aluminium<scalar>{}, silicon<scalar>{}};

    // render a rectangle mask

    // AoS
    mask<rectangle2D<>> rect2_s{0u, 0.01f * image.width(),
                                0.01f * image.height()};
    std::vector<mask<rectangle2D<>>> rect2_vec(simd_size, rect2_s);

    auto start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_s, n_samples>(image, std::move(rect2_vec),
                                              trfs_s, mat);
    auto end = std::chrono::high_resolution_clock::now();

    auto time_aos{
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
            .count() /
        1000'000.};
    std::cout << "\nRectangle AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "rectangle_AoS");

    // SoA
    std::vector<mask<rectangle2D<soa::plane_intersector>, std::uint_least16_t,
                     algebra_v<scalar>>>
        rect2_v;
    rect2_v.emplace_back(0u, 0.01f * image.width(), 0.01f * image.height());

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_v, n_samples>(image, std::move(rect2_v), trfs_v,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    auto time_soa{
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
            .count() /
        1000'000.};
    std::cout << "Rectangle SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "rectangle_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;

    // render a trapezoid mask

    // AoS
    const mask<trapezoid2D<>> trpz2_s{0u, 10.f, 30.f, 20.f, 1.f / 40.f};
    std::vector<mask<trapezoid2D<>>> trpz2_vec(simd_size, trpz2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_s, n_samples>(image, std::move(trpz2_vec),
                                              trfs_s, mat);

    end = std::chrono::high_resolution_clock::now();

    time_aos = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "\nTrapezoid AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "trapezoid_AoS");

    // SoA
    std::vector<mask<trapezoid2D<soa::plane_intersector>, std::uint_least16_t,
                     algebra_v<scalar>>>
        trap2_v;
    trap2_v.emplace_back(0u, 10.f, 30.f, 20.f, 1.f / 40.f);

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_v, n_samples>(image, std::move(trap2_v), trfs_v,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Trapezoid SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "trapezoid_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;

    // render a ring mask

    // AoS
    const mask<ring2D<>> ring2_s{0u, 12.f, 20.f};
    std::vector<mask<ring2D<>>> ring2_vec(simd_size, ring2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_s, n_samples>(image, std::move(ring2_vec),
                                              trfs_s, mat);
    end = std::chrono::high_resolution_clock::now();

    time_aos = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "\nRing AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "ring_AoS");

    // SoA
    std::vector<mask<ring2D<soa::plane_intersector>, std::uint_least16_t,
                     algebra_v<scalar>>>
        ring2_v;
    ring2_v.emplace_back(0u, 12.f, 20.f);

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_v, n_samples>(image, std::move(ring2_v), trfs_v,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Ring SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "ring_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;

    // render an annulus mask

    // AoS
    const mask<annulus2D<>> ann2_s{0u,       5.f,  13.0f, 0.74195f,
                                   1.33970f, -2.f, 2.f,   0.f};
    std::vector<mask<annulus2D<>>> ann2_vec(simd_size, ann2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_s, n_samples>(image, std::move(ann2_vec),
                                              trfs_s, mat);
    end = std::chrono::high_resolution_clock::now();

    time_aos = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "\nAnnulus AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "annulus_AoS");

    // SoA
    std::vector<mask<annulus2D<soa::plane_intersector>, std::uint_least16_t,
                     algebra_v<scalar>>>
        ann2_v;
    ann2_v.emplace_back(0u, 5.f, 13.0f, 0.74195f, 1.33970f, -2.f, 2.f, 0.f);

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_v, n_samples>(image, std::move(ann2_v), trfs_v,
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
    const mask<sphere2D<>, std::uint_least16_t, algebra_s<scalar>> sph2_s{0u,
                                                                          10.f};
    std::vector<mask<sphere2D<>>> sph2_vec(simd_size, sph2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_s, n_samples>(image, std::move(sph2_vec),
                                              trfs_s, mat);
    end = std::chrono::high_resolution_clock::now();

    time_aos = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "\nSphere AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "sphere_AoS");

    // SoA
    std::vector<mask<sphere2D<soa::sphere_intersector>, std::uint_least16_t,
                     algebra_v<scalar>>>
        sph2_v;
    sph2_v.emplace_back(0u, 10.f);

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_v, n_samples>(image, std::move(sph2_v), trfs_v,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Sphere SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "sphere_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;

    // render a line mask

    // AoS
    const mask<line<true>, std::uint_least16_t, algebra_s<scalar>> ln2_s{
        0u, 10.f, std::numeric_limits<scalar>::max()};
    std::vector<mask<line<true>>> ln2_vec(simd_size, ln2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_s, n_samples>(image, std::move(ln2_vec), trfs_s,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    time_aos = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "\nLine AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "line_AoS");

    // SoA
    std::vector<mask<line<true, soa::line_intersector>, std::uint_least16_t,
                     algebra_v<scalar>>>
        ln2_v;
    ln2_v.emplace_back(0u, 10.f, std::numeric_limits<scalar>::max());

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_v, n_samples>(image, std::move(ln2_v), trfs_v,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Line SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "line_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;

    // render a cylinder mask

    // AoS
    const mask<cylinder2D<>, std::uint_least16_t, algebra_s<scalar>> cyl2_s{
        0u, 0.5f * image.height(), 0.5f * image.width(), 0.7f * image.width()};
    std::vector<mask<cylinder2D<>>> cyl2_vec(simd_size, cyl2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_s, n_samples>(image, std::move(cyl2_vec),
                                              trfs_s, mat);
    end = std::chrono::high_resolution_clock::now();

    time_aos = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "\nCylinder AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "cylinder_AoS");

    // SoA
    std::vector<mask<cylinder2D<false, soa::cylinder_intersector>,
                     std::uint_least16_t, algebra_v<scalar>>>
        cyl2_v;
    cyl2_v.emplace_back(0u, 0.5f * image.height(), 0.5f * image.width(),
                        0.7f * image.width());

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_v, n_samples>(image, std::move(cyl2_v), trfs_v,
                                              mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Cylinder SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "cylinder_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl;

    // render a portal cylinder mask

    // AoS
    const mask<cylinder2D<false, cylinder_portal_intersector>,
               std::uint_least16_t, algebra_s<scalar>>
        pt_cyl2_s{0u, 0.5f * image.height(), 0.5f * image.width(),
                  0.7f * image.width()};
    std::vector<mask<cylinder2D<false, cylinder_portal_intersector>>>
        pt_cyl2_vec(simd_size, pt_cyl2_s);

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_s, n_samples>(image, std::move(pt_cyl2_vec),
                                              trfs_s, mat);
    end = std::chrono::high_resolution_clock::now();

    time_aos = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "\nPortal Cylinder AoS: " << time_aos << " ms" << std::endl;

    ppm.write(image, "portal_cylinder_AoS");

    // SoA
    std::vector<mask<cylinder2D<false, soa::cylinder_portal_intersector>,
                     std::uint_least16_t, algebra_v<scalar>>>
        pt_cyl2_v;
    pt_cyl2_v.emplace_back(0u, 0.5f * image.height(), 0.5f * image.width(),
                           0.7f * image.width());

    start = std::chrono::high_resolution_clock::now();
    render_mask<scalar, algebra_v, n_samples>(image, std::move(pt_cyl2_v),
                                              trfs_v, mat);
    end = std::chrono::high_resolution_clock::now();

    time_soa = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
                   .count() /
               1000'000.;
    std::cout << "Portal Cylinder SoA: " << time_soa << " ms" << std::endl;

    ppm.write(image, "portal_cylinder_SoA");

    std::cout << "Speedup: " << time_aos / time_soa << std::endl << std::endl;

    return EXIT_SUCCESS;
}
