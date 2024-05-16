/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/build_toy_detector.hpp"
#include "detray/io/frontend/detector_writer.hpp"
#include "detray/options/detector_io_options.hpp"
#include "detray/options/parse_options.hpp"
#include "detray/options/toy_detector_options.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Boost
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char **argv) {

    using namespace detray;

    // Configuration
    detray::toy_det_config toy_cfg{};
    detray::io::detector_writer_config writer_cfg{};
    writer_cfg.format(detray::io::format::json).replace_files(false);
    // Default output path
    writer_cfg.path("./toy_detector/");

    // Specific options for this test
    po::options_description desc("\nToy detector generation options");

    desc.add_options()("write_volume_graph", "Write the volume graph to file");

    po::variables_map vm =
        detray::options::parse_options(desc, argc, argv, toy_cfg, writer_cfg);

    if (vm.count("homogeneous_material") && vm.count("material_maps")) {
        std::cout << "Please specify only one material description"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (vm.count("homogeneous_material")) {
        toy_cfg.use_material_maps(false);
        writer_cfg.write_material(true);
    }
    if (vm.count("material_maps")) {
        toy_cfg.use_material_maps(true);
        writer_cfg.write_material(true);
    }

    // Build the geometry
    vecmem::host_memory_resource host_mr;
    auto [toy_det, toy_names] = build_toy_detector(host_mr, toy_cfg);

    // Write to file
    detray::io::write_detector(toy_det, toy_names, writer_cfg);

    // General options
    if (vm.count("write_volume_graph")) {
        throw std::invalid_argument("Writing of volume graph not implemented");
    }
}
