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

template <bool do_debug = false>
struct select_ray_intersector {
    template <typename mask_t>
    using type = ray_intersector<typename mask_t::shape,
                                 typename mask_t::algebra_type, do_debug>;
};

template <typename registry_t, class type_selector, std::size_t I = 0,
          typename... Fs>
consteval auto map_types(
    const types::list<Fs...>& list = detray::types::list<>{}) {
    using list_t = types::list<Fs...>;

    if constexpr (I == registry_t::n_types) {
        return list;
    } else {
        using next_type = typename registry_t::template get_type<
            static_cast<typename registry_t::id>(I)>::type;
        // Map the mask type to another type
        using mapped_t = typename type_selector::template type<next_type>;

        // Map mask position in mask store to frame type id
        if constexpr (types::contains<mapped_t, list_t>) {
            return map_types<registry_t, type_selector, I + 1u>(list);
        } else {
            return map_types<registry_t, type_selector, I + 1u>(
                types::push_back<list_t, mapped_t>{});
        }
    }
}

template <typename registry_t, class type_selector, std::size_t I = 0,
          typename... Fs>
consteval auto make_id_lookup(
    const types::list<Fs...>& list,
    std::array<dindex, registry_t::n_types> id_array = {0}) {
    using list_t = types::list<Fs...>;

    if constexpr (I == registry_t::n_types) {
        return id_array;
    } else {
        using next_type = typename registry_t::template get_type<
            static_cast<typename registry_t::id>(I)>::type;
        // Map the mask type to another type
        using mapped_t = typename type_selector::template type<next_type>;

        // Map mask position in mask store to frame type id
        if constexpr (types::contains<mapped_t, list_t>) {
            id_array[I] = types::position<mapped_t, list_t>;
            return make_id_lookup<registry_t, type_selector, I + 1u>(list,
                                                                     id_array);
        } else {
            id_array[I] = sizeof...(Fs);
            return make_id_lookup<registry_t, type_selector, I + 1u>(
                types::push_back<list_t, mapped_t>{}, id_array);
        }
    }
}

template <typename registry_t, class type_selector_t>
class mapped_type_registry {
    public:
    /// Make the type ids accessible
    using id = typename registry_t::id;
    /// Make the registered types accessible
    using types = decltype(map_types<registry_t, type_selector_t>());

    /// Conventions for some basic info
    enum : std::size_t {
        n_types = detray::types::size<types>,
        e_any = detray::types::size<types>,
        e_unknown = detray::types::size<types> + 1u,
    };

    static constexpr std::array<dindex, registry_t::n_types> id_map =
        make_id_lookup<registry_t, type_selector_t>(types{});

    DETRAY_HOST_DEVICE
    static constexpr std::size_t map(id i) {
        return id_map[static_cast<std::size_t>(i)];
    }

    /// Get the index for a type.
    /*template <typename object_t>
    DETRAY_HOST_DEVICE static consteval id get_id() {
        return to_id(detray::types::position<std::decay_t<object_t>, types>);
    }

    /// Get the index for a type. Use template parameter deduction.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval id get_id(const object_t&) {
        return get_id<object_t>();
    }*/

    /// Checks whether a given types is known in the registry.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval bool contains() {
        return detray::types::contains<std::decay_t<object_t>, types>;
    }

    /// Checks whether a given types is known in the registry.
    /// Use template parameter deduction.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval bool contains(const object_t&) {
        return contains<object_t>();
    }

    /// Checks whether a given index can be mapped to a type.
    DETRAY_HOST_DEVICE static constexpr bool is_valid(
        const std::size_t type_id) {
        return map(type_id) < n_types;
    }

    /// Convert index to ID and do some (limited) checking.
    ///
    /// @tparam ref_idx matches to index arg to perform static checks
    /// @param index argument to be converted to valid id type
    ///
    /// @return the matching ID type.
    /*template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr id to_id(const std::size_t index) {
        if (ref_idx == index) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                is_valid(ref_idx),
                "Index out of range: Please make sure that indices and type "
                "enums match the number of types in container.");
            return static_cast<id>(ref_idx);
        }
        if constexpr (ref_idx < detray::types::size<types> - 1) {
            return to_id<ref_idx + 1>(index);
        }
        // This produces a compiler error when used in type unrolling code
        return static_cast<id>(detray::types::size<types>);
    }

    /// Convert index to ID and do some (limited) checking.
    ///
    /// @tparam ref_idx matches to index arg to perform static checks
    /// @param index argument to be converted to valid id type
    ///
    /// @return the matching ID type.
    template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr std::size_t to_index(const id i) {
        if (to_id(ref_idx) == i) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                is_valid(ref_idx),
                "Index out of range: Please make sure that indices and type "
                "enums match the number of types in container.");
            return ref_idx;
        }
        if constexpr (ref_idx < detray::types::size<types> - 1) {
            return to_index<ref_idx + 1>(i);
        }
        // This produces a compiler error when used in type unrolling code
        return detray::types::size<types>;
    }

    /// Extract an index and check it.
    template <typename object_t>
    struct get_index {
        static constexpr id value = get_id<object_t>();
        DETRAY_HOST_DEVICE
        consteval bool operator()() const noexcept { return is_valid(value); }
    };*/

    /// Return a type for an index. If the index cannot be mapped, there will be
    /// a compiler error.
    template <id type_id>
    struct get_type {
        using type = detray::types::at<types, map(type_id)>;
    };
};

}  // namespace detray

/// This tests the functionality of a detector as a data store manager
GTEST_TEST(detray_core, detector) {

    using namespace detray;

    using metadata_t = test::toy_metadata;  // test::default_metadata;
    using detector_t = detector<metadata_t>;
    using mask_id = typename detector_t::masks::id;
    using material_id = typename detector_t::materials::id;
    using finder_id = typename detector_t::accel::id;

    vecmem::host_memory_resource host_mr;
    detector_t d1(host_mr);
    auto geo_ctx = typename detector_t::geometry_context{};

    using mapped_registry_t =
        mapped_type_registry<typename detector_t::masks,
                             select_ray_intersector<true>>;

    constexpr auto id_array = mapped_registry_t::id_map;

    for (std::size_t i = 0; i < id_array.size(); ++i) {
        std::cout << "i: " << i << ", id: " << id_array[i] << std::endl;
    }

    using intersector_t = typename mapped_registry_t::template get_type<
        mask_id::e_cylinder2>::type;

    types::print<types::list<intersector_t>>();

    // types::print<typename detector_t::masks::types>();

    // Helper lambda for checking the contents of an "empty" detector object.
    /*auto check_empty_detector = [](auto& d) {
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
    check_filled_detector(d3);*/
}
