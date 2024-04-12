/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algebra include(s).
#include "detray/plugins/algebra/vc_aos_definitions.hpp"
#include "detray/plugins/algebra/vc_soa_definitions.hpp"

// Detray core include(s).
#include "detray/definitions/detail/containers.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"

// Detray test include(s).
#include "detray/test/common/types.hpp"
#include "detray/test/common/utils/planes_along_direction.hpp"

// Google Benchmark include(s)
#include <benchmark/benchmark.h>

using namespace detray;

static constexpr unsigned int theta_steps{2000u};
static constexpr unsigned int phi_steps{2000u};
static constexpr unsigned int n_surfaces{16u};

/// Linear algebra implementation using SoA memory layout
using algebra_v = detray::vc_soa<test::scalar>;

/// Linear algebra implementation using AoS memory layout
using algebra_s = detray::vc_aos<test::scalar>;
// using algebra_s = detray::cmath<test::scalar>;

// Size of an SoA batch
constexpr std::size_t simd_size{dscalar<algebra_v>::size()};

using ray_t = detail::ray<algebra_s>;

namespace {

enum mask_ids : unsigned int {
    e_rectangle2 = 0,
    e_cylinder2 = 1,
    e_conc_cylinder3 = 2,
};

enum material_ids : unsigned int {
    e_slab = 0,
};

// Helper type definitions.
using mask_link_t = dtyped_index<mask_ids, dindex>;
using material_link_t = dtyped_index<material_ids, dindex>;

template <typename transform3_t>
using surface_desc_t =
    surface_descriptor<mask_link_t, material_link_t, transform3_t>;

/// Generate a number of test rays
std::vector<ray_t> generate_rays() {

    using ray_generator_t = uniform_track_generator<ray_t>;

    // Iterate through uniformly distributed momentum directions
    auto ray_generator = ray_generator_t{};
    ray_generator.config().theta_steps(theta_steps).phi_steps(phi_steps);

    std::vector<ray_t> rays;
    std::copy(ray_generator.begin(), ray_generator.end(),
              std::back_inserter(rays));

    return rays;
}

/// Generate the translation distances to place the surfaces
template <typename algebra_t>
dvector<dscalar<algebra_t>> get_dists(std::size_t n) {

    using scalar_t = dscalar<algebra_t>;

    dvector<test::scalar> dists;

    for (std::size_t i = 1u; i <= n; ++i) {
        dists.push_back(static_cast<scalar_t>(i));
    }

    return dists;
}

/// Specialization for hthe SOA memory layout (need n/simd_size samples)
template <>
dvector<dscalar<algebra_v>> get_dists<algebra_v>(std::size_t n) {

    using scalar_t = dscalar<algebra_v>;
    using value_t = typename algebra_v::value_type;

    dvector<scalar_t> dists;
    dists.resize(static_cast<std::size_t>(std::ceil(n / simd_size)));
    for (std::size_t i = 0u; i < dists.size(); ++i) {
        dists[i] = scalar_t::IndexesFromZero() +
                   scalar_t(static_cast<value_t>(i)) * simd_size +
                   scalar_t(1.f);
    }

    return dists;
}

}  // namespace

/// This benchmark runs intersection with the planar intersector
void BM_INTERSECT_PLANES_AOS(benchmark::State& state) {

    using mask_t = mask<rectangle2D, std::uint_least16_t, algebra_s>;

    auto planes = test::planes_along_direction<algebra_s>(
        get_dists<algebra_s>(n_surfaces), dvector3D<algebra_s>{1.f, 1.f, 1.f});

    constexpr mask_t rect{0u, 100.f, 200.f};
    std::vector<mask_t> masks(planes.size(), rect);

    const auto rays = generate_rays();
    const auto pi = ray_intersector<rectangle2D, algebra_s>{};

#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::size_t hit{0u};
    std::size_t miss{0u};
#endif

    for (auto _ : state) {
#ifdef DETRAY_BENCHMARK_PRINTOUTS
        hit = 0u;
        miss = 0u;
#endif

        // Iterate through uniformly distributed momentum directions
        for (const auto& ray : rays) {

            for (std::size_t i = 0u; i < planes.size(); ++i) {
                const auto& plane = planes[i];
                auto is = pi(ray, plane, masks[i], plane.transform());

#ifdef DETRAY_BENCHMARK_PRINTOUTS
                if (is.status) {
                    ++hit;
                } else {
                    ++miss;
                }
#endif

                benchmark::DoNotOptimize(is);
            }
        }
    }
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::cout << mask_t::shape::name << " AoS: hit/miss ... " << hit << " / "
              << miss << " (total: " << rays.size() * masks.size() << ")"
              << std::endl;
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_INTERSECT_PLANES_AOS)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with the planar intersector
void BM_INTERSECT_PLANES_SOA(benchmark::State& state) {

    using mask_t = mask<rectangle2D, std::uint_least16_t, algebra_v>;
    using vector3_t = dvector3D<algebra_v>;

    auto planes = test::planes_along_direction<algebra_v>(
        get_dists<algebra_v>(n_surfaces), vector3_t{1.f, 1.f, 1.f});

    std::vector<mask_t> masks{};
    for (std::size_t i = 0u; i < planes.size(); ++i) {
        masks.emplace_back(0u, 100.f, 200.f);
    }

    const auto rays = generate_rays();
    const auto pi = ray_intersector<rectangle2D, algebra_v>{};

#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::size_t hit{0u};
    std::size_t miss{0u};
#endif

    for (auto _ : state) {

#ifdef DETRAY_BENCHMARK_PRINTOUTS
        hit = 0u;
        miss = 0u;
#endif

        // Iterate through uniformly distributed momentum directions
        for (const auto& ray : rays) {

            for (std::size_t i = 0u; i < planes.size(); ++i) {
                const auto& plane = planes[i];
                auto is = pi(ray, plane, masks[i], plane.transform());

                benchmark::DoNotOptimize(is);

#ifdef DETRAY_BENCHMARK_PRINTOUTS
                hit += is.status.count();
                miss += simd_size - is.status.count();
#endif
            }
        }
    }

#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::cout << mask_t::shape::name << " SoA: hit/miss ... " << hit << " / "
              << miss << " (total: " << rays.size() * masks.size() * simd_size
              << ")" << std::endl;
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_INTERSECT_PLANES_SOA)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with the cylinder intersector
void BM_INTERSECT_CYLINDERS_AOS(benchmark::State& state) {

    using transform3_t = dtransform3D<algebra_s>;
    using scalar_t = dscalar<algebra_s>;

    using mask_t = mask<cylinder2D, std::uint_least16_t, algebra_s>;

    std::vector<mask_t> masks;
    for (const scalar_t r : get_dists<algebra_s>(n_surfaces)) {
        masks.emplace_back(0u, r, -100.f, 100.f);
    }

    mask_link_t mask_link{mask_ids::e_conc_cylinder3, 0u};
    material_link_t material_link{material_ids::e_slab, 0u};
    surface_desc_t<transform3_t> cyl_desc(
        transform3_t{}, mask_link, material_link, 0u, surface_id::e_sensitive);

    // Iterate through uniformly distributed momentum directions
    const auto rays = generate_rays();
    const auto cci = ray_intersector<cylinder2D, algebra_s>{};
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::size_t hit{0u};
    std::size_t miss{0u};
#endif
    for (auto _ : state) {
#ifdef DETRAY_BENCHMARK_PRINTOUTS
        hit = 0u;
        miss = 0u;
#endif

        // Iterate through uniformly distributed momentum directions
        for (const auto& ray : rays) {

            for (const auto& cylinder : masks) {
                auto is = cci(ray, cyl_desc, cylinder, cyl_desc.transform());

                static_assert(is.size() == 2u, "Wrong number of solutions");
#ifdef DETRAY_BENCHMARK_PRINTOUTS
                for (const auto& i : is) {
                    if (i.status) {
                        ++hit;
                    } else {
                        ++miss;
                    }
                }
#endif
                benchmark::DoNotOptimize(is);
            }
        }
    }
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::cout << mask_t::shape::name << " AoS: hit/miss ... " << hit << " / "
              << miss << " (total: " << rays.size() * masks.size() << ")"
              << std::endl;
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_INTERSECT_CYLINDERS_AOS)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with the cylinder intersector
void BM_INTERSECT_CYLINDERS_SOA(benchmark::State& state) {

    using transform3_t = dtransform3D<algebra_v>;
    using scalar_t = dscalar<algebra_v>;

    using mask_t = mask<cylinder2D, std::uint_least16_t, algebra_v>;

    std::vector<mask_t> masks;
    for (const scalar_t r : get_dists<algebra_v>(n_surfaces)) {
        masks.emplace_back(0u, r, -100.f, 100.f);
    }

    mask_link_t mask_link{mask_ids::e_conc_cylinder3, 0u};
    material_link_t material_link{material_ids::e_slab, 0u};
    surface_desc_t<transform3_t> cyl_desc(
        transform3_t{}, mask_link, material_link, 0u, surface_id::e_sensitive);

    // Iterate through uniformly distributed momentum directions
    const auto rays = generate_rays();
    const auto cci = ray_intersector<cylinder2D, algebra_v>{};
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::size_t hit{0u};
    std::size_t miss{0u};
#endif
    for (auto _ : state) {
#ifdef DETRAY_BENCHMARK_PRINTOUTS
        hit = 0u;
        miss = 0u;
#endif

        // Iterate through uniformly distributed momentum directions
        for (const auto& ray : rays) {

            for (const auto& cylinder : masks) {
                auto is = cci(ray, cyl_desc, cylinder, cyl_desc.transform());

                static_assert(is.size() == 2u, "Wrong number of solutions");
#ifdef DETRAY_BENCHMARK_PRINTOUTS
                for (const auto& i : is) {
                    hit += i.status.count();
                    miss += simd_size - i.status.count();
                }
#endif
                benchmark::DoNotOptimize(is);
            }
        }
    }
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::cout << mask_t::shape::name << " SoA: hit/miss ... " << hit << " / "
              << miss << " (total: " << rays.size() * masks.size() * simd_size
              << ")" << std::endl;
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_INTERSECT_CYLINDERS_SOA)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with the concentric cylinder intersector
void BM_INTERSECT_CONCETRIC_CYLINDERS_AOS(benchmark::State& state) {

    using transform3_t = dtransform3D<algebra_s>;
    using scalar_t = dscalar<algebra_s>;

    using mask_t = mask<concentric_cylinder2D, std::uint_least16_t, algebra_s>;

    std::vector<mask_t> masks;
    for (const scalar_t r : get_dists<algebra_s>(n_surfaces)) {
        masks.emplace_back(0u, r, -100.f, 100.f);
    }

    mask_link_t mask_link{mask_ids::e_conc_cylinder3, 0u};
    material_link_t material_link{material_ids::e_slab, 0u};
    surface_desc_t<transform3_t> cyl_desc(
        transform3_t{}, mask_link, material_link, 0u, surface_id::e_sensitive);

    // Iterate through uniformly distributed momentum directions
    const auto rays = generate_rays();
    const auto cci = ray_intersector<concentric_cylinder2D, algebra_s>{};
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::size_t hit{0u};
    std::size_t miss{0u};
#endif

    for (auto _ : state) {
#ifdef DETRAY_BENCHMARK_PRINTOUTS
        hit = 0u;
        miss = 0u;
#endif
        // Iterate through uniformly distributed momentum directions
        for (const auto& ray : rays) {

            for (const auto& cylinder : masks) {
                auto is = cci(ray, cyl_desc, cylinder, cyl_desc.transform());
#ifdef DETRAY_BENCHMARK_PRINTOUTS
                if (is.status) {
                    ++hit;
                } else {
                    ++miss;
                }
#endif
                benchmark::DoNotOptimize(is);
            }
        }
    }
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::cout << mask_t::shape::name << " AoS: hit/miss ... " << hit << " / "
              << miss << " (total: " << rays.size() * masks.size() << ")"
              << std::endl;
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_INTERSECT_CONCETRIC_CYLINDERS_AOS)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with the concentric cylinder intersector
void BM_INTERSECT_CONCETRIC_CYLINDERS_SOA(benchmark::State& state) {

    using transform3_t = dtransform3D<algebra_v>;
    using scalar_t = dscalar<algebra_v>;

    using mask_t = mask<concentric_cylinder2D, std::uint_least16_t, algebra_v>;

    std::vector<mask_t> masks;
    for (const scalar_t r : get_dists<algebra_v>(n_surfaces)) {
        masks.emplace_back(0u, r, -100.f, 100.f);
    }

    mask_link_t mask_link{mask_ids::e_conc_cylinder3, 0u};
    material_link_t material_link{material_ids::e_slab, 0u};
    surface_desc_t<transform3_t> cyl_desc(
        transform3_t{}, mask_link, material_link, 0u, surface_id::e_sensitive);

    // Iterate through uniformly distributed momentum directions
    const auto rays = generate_rays();
    const auto cci = ray_intersector<concentric_cylinder2D, algebra_v>{};
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::size_t hit{0u};
    std::size_t miss{0u};
#endif
    for (auto _ : state) {
#ifdef DETRAY_BENCHMARK_PRINTOUTS
        hit = 0u;
        miss = 0u;
#endif
        // Iterate through uniformly distributed momentum directions
        for (const auto& ray : rays) {

            for (const auto& cylinder : masks) {
                auto is = cci(ray, cyl_desc, cylinder, cyl_desc.transform());
#ifdef DETRAY_BENCHMARK_PRINTOUTS
                hit += is.status.count();
                miss += simd_size - is.status.count();
#endif
                benchmark::DoNotOptimize(is);
            }
        }
    }
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    std::cout << mask_t::shape::name << " SoA: hit/miss ... " << hit << " / "
              << miss << " (total: " << rays.size() * masks.size() * simd_size
              << ")" << std::endl;
#endif
}

BENCHMARK(BM_INTERSECT_CONCETRIC_CYLINDERS_SOA)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);
