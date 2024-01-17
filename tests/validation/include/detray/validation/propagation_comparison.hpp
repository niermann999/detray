/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/detectors/bfield.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/inspectors.hpp"
#include "detray/validation/detail/svg_display.hpp"
#include "tests/common/test_base/fixture_base.hpp"
#include "tests/common/test_base/propagator_test.hpp"
#include "propagator_cuda_kernel.hpp"

// System include(s)
#include <iostream>
#include <string>

namespace detray {

/// @brief Test class that runs the helix navigation check on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class propagation_comparison : public test::fixture_base<> {

    using scalar_t = typename detector_t::scalar_type;
    using transform3_t = typename detector_t::transform3;
    using free_track_parameters_t = free_track_parameters<transform3_t>;

    public:
    using fixture_type = test::fixture_base<>;

    struct config : public fixture_type::configuration {
        using trk_gen_config_t = typename uniform_track_generator<
            free_track_parameters_t>::configuration;

        std::string m_name{"propagation_comparison"};
        trk_gen_config_t m_trk_gen_cfg{};
        // Visualization style to be applied to the svgs
        detray::svgtools::styling::style m_style =
            detray::svgtools::styling::tableau_colorblind::style;

        /// Getters
        /// @{
        const std::string &name() const { return m_name; }
        trk_gen_config_t &track_generator() { return m_trk_gen_cfg; }
        const trk_gen_config_t &track_generator() const {
            return m_trk_gen_cfg;
        }
        const auto &svg_style() const { return m_style; }
        /// @}

        /// Setters
        /// @{
        config &name(const std::string n) {
            m_name = n;
            return *this;
        }
        /// @}
    };

    template <typename config_t>
    explicit propagation_comparison(const detector_t &det,
                              const typename detector_t::name_map &names,
                              const config_t &cfg = {})
        : m_det{det}, m_names{names} {
        m_cfg.name(cfg.name());
        m_cfg.track_generator() = cfg.track_generator();
        m_cfg.propagation() = cfg.propagation();
    }

    /// Run the check
    void TestBody() override {
        using namespace detray;
        using namespace navigation;

        // VecMem memory resource(s)
        vecmem::cuda::managed_memory_resource mng_mr;

        // Get the magnetic field
        const typename fixture_type::point3 B{0.f * unit<scalar_t>::T,
                                              0.f * unit<scalar_t>::T,
                                              2.f * unit<scalar_t>::T};
        bfield_t hom_bfield = bfield::create_const_field(B);

        // Create the vector of initial track parameterizations
        auto tracks_host = generate_tracks<typename config::track_generator_t>(mr, m_cfg.track_generator());
        vecmem::vector<free_track_parameters_t> tracks_device(tracks_host, mng_mr);

        // Host propagation
        auto &&[host_path_lengths, host_positions, host_jac_transports] =
            run_propagation_host(mng_mr, m_det, field, tracks_host);


        // Device propagation
        // Helper object for performing memory copies.
        vecmem::copy copy;

        // Get tracks data
        auto tracks_data = vecmem::get_data(tracks_host);

        // Create navigator candidates buffer
        auto candidates_buffer =
            create_candidates_buffer(m_det, tracks_host.size(), mng_mr);
        copy.setup(candidates_buffer);

        // Create vector buffer for track recording
        std::vector<std::size_t> sizes(tracks_host.size(), 0);
        std::vector<std::size_t> capacities;
        for (auto &r : host_positions) {
            capacities.push_back(r.size());
        }

        vecmem::data::jagged_vector_buffer<scalar> path_lengths_buffer(
            capacities, mng_mr, nullptr, vecmem::data::buffer_type::resizable);
        vecmem::data::jagged_vector_buffer<vector3_t> positions_buffer(
            capacities, mng_mr, nullptr, vecmem::data::buffer_type::resizable);
        vecmem::data::jagged_vector_buffer<free_matrix> jac_transports_buffer(
            capacities, mng_mr, nullptr, vecmem::data::buffer_type::resizable);

        copy.setup(path_lengths_buffer);
        copy.setup(positions_buffer);
        copy.setup(jac_transports_buffer);

        // Run the propagator test for GPU device
        covfie::field<bfield::const_bknd_t> device_field(field);
        typename detector_t::view_type det_view = detray::get_data(m_det);

        run_propagation_device<bfield::const_bknd_t, detector_t>(
            &mng_mr, m_det, det_view, device_field, tracks_data, candidates_buffer,
            path_lengths_buffer, positions_buffer, jac_transports_buffer);

        vecmem::jagged_vector<scalar> device_path_lengths(&mng_mr);
        vecmem::jagged_vector<vector3_t> device_positions(&mng_mr);
        vecmem::jagged_vector<free_matrix> device_jac_transports(&mng_mr);

        copy(path_lengths_buffer, device_path_lengths);
        copy(positions_buffer, device_positions);
        copy(jac_transports_buffer, device_jac_transports);

        /// Error statistic
        /*std::size_t n_close_miss{0u}, n_fatal{0u};

        std::ios_base::openmode io_mode = std::ios::trunc | std::ios::out;
        detray::io::detail::file_handle debug_file{"./propagation_comparison", ".txt",
                                                   io_mode};

        bool success &= detail::compare_traces(intersection_trace,
                                            obj_tracer, helix, n_tracks,
                                            trk_state_generator.size());
            if (not success) {
                // Write debug info to file
                *debug_file << "HELIX " << n_tracks << ":\n\n"
                            << nav_printer.to_string()
                            << step_printer.to_string();

                // Create the svg for failed tracks.
                detray::svgtools::illustrator il{m_det, m_names,
                                                 m_cfg.svg_style()};
                il.show_info(true);
                il.hide_eta_lines(true);
                il.hide_portals(false);
                il.hide_passives(false);

                detail::svg_display(gctx, il, intersection_trace, helix,
                                    "helix_" + std::to_string(n_tracks),
                                    m_cfg.name(), obj_tracer.object_trace);

                // Keep a statistic on the errors that occured
                if (!propagation._navigation.is_complete()) {
                    ++n_fatal;
                } else {
                    // @TODO: Check mask boundaries
                    ++n_close_miss;
                }
            }

            EXPECT_TRUE(success)
                << "\nFailed on helix " << tracks_host << "/"
                << trk_state_generator.size() << ": " << helix << "\n\n";

        if (n_close_miss > 0u || n_fatal > 0u) {
            std::cout << "-----------------------------------"
                      << "Error Statistic:\n\n"
                      << "\n (close misses: " << n_close_miss
                      << ", fatal failures: " << n_fatal << ")\n"
                      << "-----------------------------------\n"
                      << std::endl;
        }*/
    }

    private:
    /// The configuration of this test
    config m_cfg;
    /// The detector to be checked
    const detector_t &m_det;
    /// Volume names
    const typename detector_t::name_map &m_names;
};

}  // namespace detray
