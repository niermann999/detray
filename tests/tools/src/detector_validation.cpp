/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/options/detector_io_options.hpp"
#include "detray/options/parse_options.hpp"
#include "detray/options/propagation_options.hpp"
#include "detray/options/track_generator_options.hpp"
#include "detray/test/detail/register_checks.hpp"
#include "detray/test/detail/whiteboard.hpp"
#include "detray/test/detector_consistency.hpp"
#include "detray/test/detector_scan.hpp"
#include "detray/test/navigation_validation.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// Boost
#include <boost/program_options.hpp>

// System include(s)
#include <sstream>
#include <stdexcept>
#include <string>

namespace po = boost::program_options;
using namespace detray;

int main(int argc, char** argv) {

    // Use the most general type to be able to read in all detector files
    using detector_t = detray::detector<>;
    using scalar_t = typename detector_t::scalar_type;

    // Filter out the google test flags
    ::testing::InitGoogleTest(&argc, argv);

    // Specific options for this test
    po::options_description desc("\ndetray detector validation options");

    desc.add_options()("write_volume_graph", "Write the volume graph to file")(
        "write_scan_data", "Write the ray/helix scan data to file");

    // Configs to be filled
    detray::io::detector_reader_config reader_cfg{};
    reader_cfg.do_check(false);  // < Don't run consistency check twice
    detray::test::consistency_check<detector_t>::config con_chk_cfg{};
    detray::test::ray_scan<detector_t>::config ray_scan_cfg{};
    detray::test::helix_scan<detector_t>::config hel_scan_cfg{};
    detray::test::straight_line_navigation<detector_t>::config str_nav_cfg{};
    detray::test::helix_navigation<detector_t>::config hel_nav_cfg{};

    po::variables_map vm = detray::options::parse_options(
        desc, argc, argv, reader_cfg, hel_scan_cfg.track_generator(),
        hel_nav_cfg.propagation());

    // General options
    if (vm.count("write_volume_graph")) {
        con_chk_cfg.write_graph(true);
        throw std::invalid_argument("Writing of volume graph not implemented");
    }
    if (vm.count("write_scan_data")) {
        ray_scan_cfg.write_intersections(true);
        hel_scan_cfg.write_intersections(true);
    }

    // For now: Copy the options to the other tests
    ray_scan_cfg.track_generator() = hel_scan_cfg.track_generator();
    str_nav_cfg.propagation() = hel_nav_cfg.propagation();

    vecmem::host_memory_resource host_mr;

    const auto [det, names] =
        detray::io::read_detector<detector_t>(host_mr, reader_cfg);
    const std::string& det_name = det.name(names);

    // Create the whiteboard for data transfer between the steps
    auto white_board = std::make_shared<test::whiteboard>();
    ray_scan_cfg.name(det_name + "_ray_scan");
    ray_scan_cfg.whiteboard(white_board);
    ray_scan_cfg.intersection_file(det_name + "_ray_scan_intersections.csv");
    ray_scan_cfg.track_param_file(det_name + "_ray_scan_track_parameters.csv");

    hel_scan_cfg.name(det_name + "_helix_scan");
    hel_scan_cfg.whiteboard(white_board);
    // Let the Newton algorithm dynamically choose tol. based on approx. error
    hel_scan_cfg.mask_tolerance({detray::detail::invalid_value<scalar_t>(),
                                 detray::detail::invalid_value<scalar_t>()});
    hel_scan_cfg.intersection_file(det_name + "_helix_scan_intersections.csv");
    hel_scan_cfg.track_param_file(det_name +
                                  "_helix_scan_track_parameters.csv");

    str_nav_cfg.whiteboard(white_board);
    hel_nav_cfg.whiteboard(white_board);

    // General data consistency of the detector
    detray::detail::register_checks<detray::test::consistency_check>(
        det, names, con_chk_cfg);

    // Navigation link consistency, discovered by ray intersection
    ray_scan_cfg.name(det_name + "_ray_scan");
    detray::detail::register_checks<detray::test::ray_scan>(det, names,
                                                            ray_scan_cfg);

    // Comparision of straight line navigation with ray scan
    str_nav_cfg.name(det_name + "_straight_line_navigation");
    // Ensure that the same mask tolerance is used
    auto mask_tolerance = ray_scan_cfg.mask_tolerance();
    str_nav_cfg.propagation().navigation.min_mask_tolerance =
        static_cast<float>(mask_tolerance[0]);
    str_nav_cfg.propagation().navigation.max_mask_tolerance =
        static_cast<float>(mask_tolerance[1]);
    detray::detail::register_checks<detray::test::straight_line_navigation>(
        det, names, str_nav_cfg);

    // Navigation link consistency, discovered by helix intersection
    hel_scan_cfg.name(det_name + "_helix_scan");
    detray::detail::register_checks<detray::test::helix_scan>(det, names,
                                                              hel_scan_cfg);

    // Comparision of navigation in a constant B-field with helix
    hel_nav_cfg.name(det_name + "_helix_navigation");
    detray::detail::register_checks<detray::test::helix_navigation>(
        det, names, hel_nav_cfg);

    // Run the checks
    return RUN_ALL_TESTS();
}
