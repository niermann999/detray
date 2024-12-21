/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/type_list.hpp"

// System include(s)
#include <type_traits>
#include <utility>

namespace detray {

/// @brief match types with indices and vice versa.
///
/// @tparam IDs enum that references the types (not used in base class)
/// @tparam registered_types the types that can be mapped to indices
template <class ID, typename... registered_types>
class type_registry {
    public:
    /// Make the type IDs accessible
    using id = ID;
    /// Make the registered types accessible
    using types = detray::types::list<registered_types...>;

    /// Conventions for some basic info
    enum : std::size_t {
        n_types = sizeof...(registered_types),
        e_any = sizeof...(registered_types),
        e_unknown = sizeof...(registered_types) + 1,
    };

    /// Get the index for a type.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval ID get_id() {
        return to_id(detray::types::position<std::decay_t<object_t>, types>);
    }

    /// Get the index for a type. Use template parameter deduction.
    template <typename object_t>
    DETRAY_HOST_DEVICE static consteval ID get_id(const object_t&) {
        return get_id<object_t>();
    }

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
        return type_id < n_types;
    }

    /// Convert index to ID and do some (limited) checking.
    ///
    /// @tparam ref_idx matches to index arg to perform static checks
    /// @param index argument to be converted to valid id type
    ///
    /// @return the matching ID type.
    template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr ID to_id(const std::size_t index) {
        if (ref_idx == index) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                is_valid(ref_idx),
                "Index out of range: Please make sure that indices and type "
                "enums match the number of types in container.");
            return static_cast<ID>(ref_idx);
        }
        if constexpr (ref_idx < sizeof...(registered_types) - 1) {
            return to_id<ref_idx + 1>(index);
        }
        // This produces a compiler error when used in type unrolling code
        return static_cast<ID>(sizeof...(registered_types));
    }

    /// Convert index to ID and do some (limited) checking.
    ///
    /// @tparam ref_idx matches to index arg to perform static checks
    /// @param index argument to be converted to valid id type
    ///
    /// @return the matching ID type.
    template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr std::size_t to_index(const ID id) {
        if (to_id(ref_idx) == id) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                is_valid(ref_idx),
                "Index out of range: Please make sure that indices and type "
                "enums match the number of types in container.");
            return ref_idx;
        }
        if constexpr (ref_idx < sizeof...(registered_types) - 1) {
            return to_index<ref_idx + 1>(id);
        }
        // This produces a compiler error when used in type unrolling code
        return sizeof...(registered_types);
    }

    /// Extract an index and check it.
    template <typename object_t>
    struct get_index {
        static constexpr ID value = get_id<object_t>();
        DETRAY_HOST_DEVICE
        consteval bool operator()() const noexcept { return is_valid(value); }
    };

    /// Return a type for an index. If the index cannot be mapped, there will be
    /// a compiler error.
    template <ID type_id>
    struct get_type {
        using type = detray::types::at<types, to_index(type_id)>;
    };
};

}  // namespace detray
