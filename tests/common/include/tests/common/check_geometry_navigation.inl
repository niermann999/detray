/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <sstream>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/track_generators.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/inspectors.hpp"
#include "tests/common/tools/particle_gun.hpp"

using namespace detray;
using namespace navigation;

namespace {

using transform3_t = __plugin::transform3<scalar>;
using ray_type = detail::ray<transform3_t>;
using free_track_parameters_type = free_track_parameters<transform3_t>;

}  // anonymous namespace

/// This test runs intersection with all portals of the toy detector with a ray
/// and then compares the intersection trace with a straight line navigation.
TEST(ALGEBRA_PLUGIN, straight_line_navigation) {

    // Detector configuration
    constexpr std::size_t n_brl_layers{4};
    constexpr std::size_t n_edc_layers{7};
    vecmem::host_memory_resource host_mr;
    auto det = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    // Straight line navigation
    using detector_t = decltype(det);
    using intersection_t =
        line_plane_intersection<typename detector_t::surface_type,
                                transform3_t>;
    using object_tracer_t =
        object_tracer<intersection_t, dvector, status::e_on_module,
                      status::e_on_portal>;
    using inspector_t = aggregate_inspector<object_tracer_t, print_inspector>;
    using navigator_t = navigator<detector_t, inspector_t>;
    using stepper_t =
        line_stepper<transform3_t, unconstrained_step, always_init>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain<>>;

    // Propagator
    propagator_t prop(stepper_t{}, navigator_t{});

    constexpr std::size_t theta_steps{50};
    constexpr std::size_t phi_steps{50};

    const point3 ori{0., 0., 0.};
    // det.volume_by_pos(ori).index();

    // Iterate through uniformly distributed momentum directions
    for (const auto ray :
         uniform_track_generator<ray_type>(theta_steps, phi_steps, ori)) {

        // Shoot ray through the detector and record all surfaces it encounters
        const auto intersection_trace =
            particle_gun::shoot_particle(det, ray);  // :)

        // Now follow that ray with a track and check, if we find the same
        // volumes and distances along the way
        free_track_parameters_type track(ray.pos(), 0, ray.dir(), -1);
        propagator_t::state propagation(track, det);

        // Retrieve navigation information
        auto &inspector = propagation._navigation.inspector();
        auto &obj_tracer = inspector.template get<object_tracer_t>();
        auto &debug_printer = inspector.template get<print_inspector>();

        ASSERT_TRUE(prop.propagate(propagation)) << debug_printer.to_string();

        // Compare intersection records
        EXPECT_EQ(obj_tracer.object_trace.size(), intersection_trace.size());

        std::stringstream debug_stream;
        for (std::size_t intr_idx = 0; intr_idx < intersection_trace.size();
             ++intr_idx) {
            debug_stream << "-------Intersection trace\n"
                         << "ray gun: "
                         << "\tvol id: " << intersection_trace[intr_idx].first
                         << ", " << intersection_trace[intr_idx].second;
            debug_stream << "navig.: " << obj_tracer[intr_idx];
        }

        // Check every single recorded intersection
        for (std::size_t i = 0; i < obj_tracer.object_trace.size(); ++i) {
            if (obj_tracer[i].surface.barcode() !=
                intersection_trace[i].second.surface.barcode()) {
                // Intersection record at portal bound might be flipped
                // (the portals overlap completely)
                if (obj_tracer[i].surface.barcode() ==
                        intersection_trace[i + 1].second.surface.barcode() and
                    obj_tracer[i + 1].surface.barcode() ==
                        intersection_trace[i].second.surface.barcode()) {
                    // Have already checked the next record
                    ++i;
                    continue;
                }
            }
            EXPECT_EQ(obj_tracer[i].surface.barcode(),
                      intersection_trace[i].second.surface.barcode())
                << debug_printer.to_string() << debug_stream.str();
        }
    }
}

/// Check the Runge-Kutta based navigation against a helix trajectory as ground
/// truth
TEST(ALGEBRA_PLUGIN, helix_navigation) {
    using namespace navigation;

    // Detector configuration
    constexpr std::size_t n_brl_layers{4};
    constexpr std::size_t n_edc_layers{7};
    vecmem::host_memory_resource host_mr;

    using b_field_t = decltype(create_toy_geometry(
        std::declval<vecmem::host_memory_resource &>(), n_brl_layers,
        n_edc_layers))::bfield_type;

    const vector3 B{0. * unit<scalar>::T, 0. * unit<scalar>::T,
                    2. * unit<scalar>::T};

    auto det = create_toy_geometry(
        host_mr,
        b_field_t(b_field_t::backend_t::configuration_t{B[0], B[1], B[2]}),
        n_brl_layers, n_edc_layers);

    // Runge-Kutta based navigation
    using detector_t = decltype(det);
    using intersection_t =
        line_plane_intersection<typename detector_t::surface_type,
                                transform3_t>;
    using object_tracer_t =
        object_tracer<intersection_t, dvector, status::e_on_module,
                      status::e_on_portal>;
    using inspector_t = aggregate_inspector<object_tracer_t, print_inspector>;
    using navigator_t = navigator<detector_t, inspector_t>;
    using stepper_t = rk_stepper<b_field_t::view_t, transform3_t,
                                 unconstrained_step, always_init>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain<>>;

    // Propagator
    propagator_t prop(stepper_t{}, navigator_t{});

    constexpr std::size_t theta_steps{10};
    constexpr std::size_t phi_steps{10};

    // det.volume_by_pos(ori).index();
    const point3 ori{0., 0., 0.};
    const scalar p_mag{10. * unit<scalar>::GeV};

    // Overstepping
    constexpr scalar overstep_tol{-7. * unit<scalar>::um};

    // Iterate through uniformly distributed momentum directions
    for (auto track : uniform_track_generator<free_track_parameters_type>(
             theta_steps, phi_steps, ori, p_mag)) {
        // Prepare for overstepping in the presence of b fields
        track.set_overstep_tolerance(overstep_tol);

        // Get ground truth helix from track
        detail::helix helix(track, &B);

        // Shoot ray through the detector and record all surfaces it encounters
        const auto intersection_trace =
            particle_gun::shoot_particle(det, helix);

        // Now follow that helix with the same track and check, if we find
        // the same volumes and distances along the way
        propagator_t::state propagation(track, det.get_bfield(), det);

        // Retrieve navigation information
        auto &inspector = propagation._navigation.inspector();
        auto &obj_tracer = inspector.template get<object_tracer_t>();
        auto &debug_printer = inspector.template get<print_inspector>();

        ASSERT_TRUE(prop.propagate(propagation)) << debug_printer.to_string();

        std::stringstream debug_stream;
        for (std::size_t intr_idx = 0; intr_idx < intersection_trace.size();
             ++intr_idx) {
            debug_stream << "-------Intersection trace\n"
                         << "helix gun: "
                         << "\tvol id: " << intersection_trace[intr_idx].first
                         << ", " << intersection_trace[intr_idx].second;
            debug_stream << "navig.: " << obj_tracer[intr_idx];
        }

        // Compare intersection records
        EXPECT_EQ(obj_tracer.object_trace.size(), intersection_trace.size())
            << debug_printer.to_string() << debug_stream.str();

        // Check every single recorded intersection
        for (std::size_t i = 0; i < obj_tracer.object_trace.size(); ++i) {
            if (obj_tracer[i].surface.barcode() !=
                intersection_trace[i].second.surface.barcode()) {
                // Intersection record at portal bound might be flipped
                // (the portals overlap completely)
                if (obj_tracer[i].surface.barcode() ==
                        intersection_trace[i + 1].second.surface.barcode() and
                    obj_tracer[i + 1].surface.barcode() ==
                        intersection_trace[i].second.surface.barcode()) {
                    // Have already checked the next record
                    ++i;
                    continue;
                }
            }
            EXPECT_EQ(obj_tracer[i].surface.barcode(),
                      intersection_trace[i].second.surface.barcode());
        }
    }
}