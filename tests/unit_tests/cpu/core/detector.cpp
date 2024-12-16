/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"

#include "detray/definitions/detail/indexing.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/utils/type_list.hpp"

// Detray test include(s)
#include "detray/test/utils/prefill_detector.hpp"
#include "detray/test/utils/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

namespace detray {

template <typename detector_t, std::size_t I = 0, typename... Fs>
consteval auto make_frame_type_set(
    const types::list<Fs...>& list = {},
    std::array<dindex, detector_t::masks::n_types> id_array = {0}) {
    using frame_list_t = types::list<Fs...>;

    if constexpr (I == detector_t::masks::n_types) {
        return std::make_tuple(list, id_array);
    } else {
        using algebra_t = typename detector_t::algebra_type;
        using next_mask_t =
            typename detector_t::mask_container::template get_type<
                static_cast<typename detector_t::masks::id>(I)>;
        using frame_t = ray_intersector<typename next_mask_t::shape, algebra_t>;

        // Map mask position in mask store to frame type id

        // Coordinate frame type already registered?
        if constexpr (types::contains<frame_t, frame_list_t>) {
            id_array[I] = types::position<frame_t, frame_list_t>;
            return make_frame_type_set<detector_t, I + 1u>(list, id_array);
        } else {
            id_array[I] = sizeof...(Fs);
            return make_frame_type_set<detector_t, I + 1u>(
                types::push_back<frame_list_t, frame_t>{}, id_array);
        }
    }
}

}  // namespace detray

/// This tests the functionality of a detector as a data store manager
GTEST_TEST(detray_core, detector) {

    using namespace detray;

    using metadata_t = test::default_metadata;
    using detector_t = detector<metadata_t>;
    using mask_id = typename detector_t::masks::id;
    using material_id = typename detector_t::materials::id;
    using finder_id = typename detector_t::accel::id;

    vecmem::host_memory_resource host_mr;
    detector_t d1(host_mr);
    auto geo_ctx = typename detector_t::geometry_context{};

    auto [type_list, id_array] =
        make_frame_type_set<detector_t>(types::list<>{});

    for (std::size_t i = 0; i < id_array.size(); ++i) {
        std::cout << "i: " << i << ", id: " << id_array[i] << std::endl;
    }

    types::print<typename detector_t::masks::types>();
    types::print<decltype(type_list)>();

    // Helper lambda for checking the contents of an "empty" detector object.
    auto check_empty_detector = [](auto& d) {
        EXPECT_TRUE(d.volumes().empty());
        EXPECT_TRUE(d.portals().empty());
        EXPECT_TRUE(d.transform_store().empty());
        EXPECT_TRUE(d.mask_store().template empty<mask_id::e_rectangle2>());
        EXPECT_TRUE(
            d.mask_store().template empty<mask_id::e_portal_rectangle2>());
        EXPECT_TRUE(d.mask_store().template empty<mask_id::e_trapezoid2>());
        EXPECT_TRUE(d.mask_store().template empty<mask_id::e_annulus2>());
        EXPECT_TRUE(d.mask_store().template empty<mask_id::e_cylinder2>());
        EXPECT_TRUE(
            d.mask_store().template empty<mask_id::e_portal_cylinder2>());
        EXPECT_TRUE(d.mask_store().template empty<mask_id::e_ring2>());
        EXPECT_TRUE(d.mask_store().template empty<mask_id::e_portal_ring2>());
        EXPECT_TRUE(d.mask_store().template empty<mask_id::e_straw_tube>());
        EXPECT_TRUE(d.mask_store().template empty<mask_id::e_drift_cell>());
        EXPECT_TRUE(d.material_store().template empty<material_id::e_slab>());
        EXPECT_TRUE(d.material_store().template empty<material_id::e_rod>());
        EXPECT_TRUE(
            d.accelerator_store().template empty<finder_id::e_brute_force>());
        EXPECT_TRUE(
            d.accelerator_store().template empty<finder_id::e_disc_grid>());
        EXPECT_TRUE(d.accelerator_store()
                        .template empty<finder_id::e_cylinder2_grid>());
        EXPECT_TRUE(
            d.accelerator_store().template empty<finder_id::e_irr_disc_grid>());
        EXPECT_TRUE(d.accelerator_store()
                        .template empty<finder_id::e_irr_cylinder2_grid>());
        EXPECT_TRUE(
            d.accelerator_store().template empty<finder_id::e_default>());
    };

    // Check the empty detector object.
    check_empty_detector(d1);

    // Add some geometrical data
    prefill_detector(d1, geo_ctx);

    // Helper lambda for checking the contents of a "filled" detector object.
    auto check_filled_detector = [](auto& d) {
        // TODO: add B-field check
        EXPECT_EQ(d.volumes().size(), 1u);
        EXPECT_EQ(d.portals().size(), 3u);
        EXPECT_EQ(d.transform_store().size(), 4u);
        EXPECT_EQ(d.mask_store().template size<mask_id::e_rectangle2>(), 1u);
        EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_rectangle2>(),
                  1u);
        EXPECT_EQ(d.mask_store().template size<mask_id::e_trapezoid2>(), 1u);
        EXPECT_EQ(d.mask_store().template size<mask_id::e_annulus2>(), 1u);
        EXPECT_EQ(d.mask_store().template size<mask_id::e_cylinder2>(), 0u);
        EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_cylinder2>(),
                  0u);
        EXPECT_EQ(d.mask_store().template size<mask_id::e_ring2>(), 0u);
        EXPECT_EQ(d.mask_store().template size<mask_id::e_portal_ring2>(), 0u);
        EXPECT_EQ(d.mask_store().template size<mask_id::e_straw_tube>(), 0u);
        EXPECT_EQ(d.mask_store().template size<mask_id::e_drift_cell>(), 0u);
        EXPECT_EQ(d.material_store().template size<material_id::e_slab>(), 2u);
        EXPECT_EQ(d.material_store().template size<material_id::e_rod>(), 1u);
        EXPECT_EQ(
            d.accelerator_store().template size<finder_id::e_brute_force>(),
            1u);
        EXPECT_EQ(d.accelerator_store().template size<finder_id::e_disc_grid>(),
                  0u);
        EXPECT_EQ(
            d.accelerator_store().template size<finder_id::e_cylinder2_grid>(),
            0u);
        EXPECT_EQ(
            d.accelerator_store().template size<finder_id::e_irr_disc_grid>(),
            0u);
        EXPECT_EQ(d.accelerator_store()
                      .template size<finder_id::e_irr_cylinder2_grid>(),
                  0u);
        EXPECT_EQ(d.accelerator_store().template size<finder_id::e_default>(),
                  1u);
    };

    // Check the filled detector object.
    check_filled_detector(d1);

    // Move construct a detector object.
    detector_t d2{std::move(d1)};
    check_filled_detector(d2);

    // Create a new, empty detector.
    detector_t d3{host_mr};
    check_empty_detector(d3);

    // Move assign the filled detector to the empty one.
    d3 = std::move(d2);
    check_filled_detector(d3);
}
