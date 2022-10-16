/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/core/detail/container_views.hpp"
#include "detray/core/detail/new_tuple_container.hpp"
#include "detray/core/detail/type_registry.hpp"
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <type_traits>

namespace detray {

/// Placeholder context type
struct empty_context {};

/// @brief Wraps a vecmem enabled tuple and adds functionality to handle data
/// collections.
///
/// @tparam An enum of type IDs that needs to match the value types of the
/// @c Ts pack.
/// @tparam context_t How to retrieve data according to e.g. conditions data
/// @tparam tuple_t The type of the underlying tuple container.
/// @tparam container_t The type of container to use for the respective
///                     data collections.
/// @tparam Ts the data types (value types of the collections)
template <typename ID = std::size_t, typename context_t = empty_context,
          template <typename...> class tuple_t = dtuple,
          template <typename...> class container_t = dvector, typename... Ts>
class multi_type_store {

    public:
    using size_type = typename container_t<detail::first_t<Ts...>>::size_type;
    using context_type = context_t;

    /// How to find a data collection in the store
    /// @{
    using ids = ID;
    template <typename index_t>
    using link_type = dtyped_index<ID, index_t>;
    using single_link = link_type<dindex>;
    using range_link = link_type<dindex_range>;
    /// @}

    /// Allow matching between IDs and collection value types
    /// @{
    using value_types = detail::registry_base<ID, true, Ts...>;
    template <ID id>
    using get_type = typename value_types::template get_type<id>::type;
    /// @}

    /// Underlying tuple container that can handle vecmem views
    using tuple_type = detail::tuple_container<tuple_t, container_t<Ts>...>;
    /// Vecmem view types
    using view_type = typename tuple_type::view_type;
    using const_view_type = typename tuple_type::const_view_type;

    /// Empty container
    constexpr multi_type_store() = default;

    // Delegate constructors to tuple container, which handles the memory

    /// Copy construct from element types
    constexpr explicit multi_type_store(const Ts &...args)
        : m_tuple_container(args...) {}

    /// Construct with a specific vecmem memory resource @param resource
    /// (host-side only)
    template <typename allocator_t = vecmem::memory_resource,
              std::enable_if_t<not detail::is_device_view_v<allocator_t>,
                               bool> = true>
    DETRAY_HOST explicit multi_type_store(allocator_t &resource)
        : m_tuple_container(resource) {}

    /// Copy Construct with a specific (vecmem) memory resource @param resource
    /// (host-side only)
    template <
        typename allocator_t = vecmem::memory_resource,
        typename T = tuple_t<Ts...>,
        std::enable_if_t<std::is_same_v<T, std::tuple<Ts...>>, bool> = true>
    DETRAY_HOST explicit multi_type_store(allocator_t &resource,
                                          const Ts &...args)
        : m_tuple_container(resource, args...) {}

    /// Construct from the container @param view . Mainly used device-side.
    template <
        typename tuple_view_t,
        std::enable_if_t<detail::is_device_view_v<tuple_view_t>, bool> = true>
    DETRAY_HOST_DEVICE multi_type_store(tuple_view_t &view)
        : m_tuple_container(view) {}

    /// @returns a pointer to the underlying tuple container - const
    DETRAY_HOST_DEVICE
    constexpr auto data() const -> const tuple_type * {
        return &m_tuple_container;
    }

    /// @returns a pointer to the underlying tuple container - non-const
    DETRAY_HOST_DEVICE
    constexpr auto data() -> tuple_type * { return &m_tuple_container; }

    /// Convenience functionality: Only available for a store with a single
    /// collection.
    /// @{

    /// Elementwise access. Needs @c operator[] for storage type - non-const
    DETRAY_HOST_DEVICE
    constexpr auto operator[](const dindex i) {
        return detail::get<0>(m_tuple_container)[i];
    }

    /// Elementwise access. Needs @c operator[] for storage type - const
    DETRAY_HOST_DEVICE
    constexpr auto operator[](const dindex i) const {
        return detail::get<0>(m_tuple_container)[i];
    }

    /// @}

    /// @returns the size of the underlying tuple
    DETRAY_HOST_DEVICE
    constexpr auto n_collections() const -> std::size_t {
        return sizeof...(Ts);
    }

    /// @returns the collections iterator at the start position.
    template <ID id = ID{0}>
    DETRAY_HOST_DEVICE constexpr auto begin(const context_type & /*ctx*/ = {}) {
        return detail::get<value_types::to_index(id)>(m_tuple_container)
            .begin();
    }

    /// @returns the collections iterator sentinel.
    template <ID id = ID{0}>
    DETRAY_HOST_DEVICE constexpr auto end(const context_type & /*ctx*/ = {}) {
        return detail::get<value_types::to_index(id)>(m_tuple_container).end();
    }

    /// @returns a data collection by @tparam ID - const
    template <ID id>
    DETRAY_HOST_DEVICE constexpr decltype(auto) get(
        const context_type & /*ctx*/ = {}) const noexcept {
        return detail::get<value_types::to_index(id)>(m_tuple_container);
    }

    /// @returns a data collection by @tparam ID - non-const
    template <ID id>
    DETRAY_HOST_DEVICE constexpr decltype(auto) get(
        const context_type & /*ctx*/ = {}) noexcept {
        return detail::get<value_types::to_index(id)>(m_tuple_container);
    }

    /// @returns the size of a data collection by @tparam ID
    template <ID id>
    DETRAY_HOST_DEVICE constexpr auto size(
        const context_type & /*ctx*/ = {}) const -> std::size_t {
        return detail::get<value_types::to_index(id)>(m_tuple_container).size();
    }

    /// @returns true if the collection given by @tparam ID is empty
    template <ID id>
    DETRAY_HOST_DEVICE constexpr auto empty(
        const context_type & /*ctx*/ = {}) const -> bool {
        return detail::get<value_types::to_index(id)>(m_tuple_container)
            .empty();
    }

    /// Reserve memory of size @param n for a collection given by @tparam id
    template <ID id>
    DETRAY_HOST void reserve(std::size_t n, const context_type & /*ctx*/) {
        detail::get<value_types::to_index(id)>(m_tuple_container).reserve(n);
    }

    /// Add a new element to a collection
    ///
    /// @tparam ID is the id of the collection
    /// @tparam Args are the types of the constructor arguments
    ///
    /// @param args is the list of constructor arguments
    ///
    /// @note in general can throw an exception
    template <ID id, typename T>
    DETRAY_HOST constexpr auto push_back(
        const T &arg, const context_type & /*ctx*/ = {}) noexcept(false)
        -> void {
        auto &coll = detail::get<value_types::to_index(id)>(m_tuple_container);
        return coll.push_back(arg);
    }

    /// Add a new element to a collection in place
    ///
    /// @tparam ID is the id of the collection
    /// @tparam Args are the types of the constructor arguments
    ///
    /// @param args is the list of constructor arguments
    ///
    /// @note in general can throw an exception
    template <ID id, typename... Args>
    DETRAY_HOST constexpr decltype(auto) emplace_back(
        const context_type & /*ctx*/ = {}, Args &&...args) noexcept(false) {
        auto &coll = detail::get<value_types::to_index(id)>(m_tuple_container);
        return coll.emplace_back(std::forward<Args>(args)...);
    }

    /// Add a collection - copy
    ///
    /// @tparam T is the value type of the collection
    ///
    /// @param new_data is the new collection to be added
    ///
    /// @note in general can throw an exception
    template <typename T>
    DETRAY_HOST auto insert(container_t<T> &new_data,
                            const context_type & /*ctx*/ = {}) noexcept(false)
        -> void {

        static_assert((std::is_same_v<T, Ts> || ...) == true,
                      "The type is not included in the parameter pack.");

        auto &coll = const_cast<container_t<T> &>(
            detail::get<container_t<T>>(m_tuple_container));

        coll.reserve(coll.size() + new_data.size());
        coll.insert(coll.end(), new_data.begin(), new_data.end());
    }

    /// Add a new collection - move
    ///
    /// @tparam T is the value type of the collection
    ///
    /// @param new_data is the new collection to be added
    ///
    /// @note in general can throw an exception
    template <typename T>
    DETRAY_HOST auto insert(container_t<T> &&new_data,
                            const context_type & /*ctx*/ = {}) noexcept(false)
        -> void {

        static_assert((std::is_same_v<T, Ts> || ...) == true,
                      "The type is not included in the parameter pack.");

        auto &coll = detail::get<container_t<T>>(m_tuple_container);

        coll.reserve(coll.size() + new_data.size());
        coll.insert(coll.end(), std::make_move_iterator(new_data.begin()),
                    std::make_move_iterator(new_data.end()));
    }

    /// Append another store to the current one
    ///
    /// @tparam current_idx is the index to start unrolling
    ///
    /// @param other The other container
    ///
    /// @note in general can throw an exception
    template <std::size_t current_idx = 0>
    DETRAY_HOST void append(multi_type_store &other) noexcept(false) {
        auto &coll = other.template get<value_types::to_id(current_idx)>();
        insert(coll);

        if constexpr (current_idx < sizeof...(Ts) - 1) {
            append<current_idx + 1>(other);
        }
    }

    /// Calls a functor with a specific data collection (given by ID).
    ///
    /// @tparam functor_t functor that will be called on the group.
    /// @tparam Args argument types for the functor
    ///
    /// @param id the element id
    /// @param args additional functor arguments
    ///
    /// @return the functor output
    template <typename functor_t, typename... Args>
    DETRAY_HOST_DEVICE typename functor_t::output_type call(const ID id,
                                                            Args &&...args) {
        return m_tuple_container.template call<functor_t>(
            static_cast<size_type>(id), std::forward<Args>(args)...);
    }

    /// Calls a functor with a specific element of a data collection
    /// (given by a link).
    ///
    /// @tparam functor_t functor that will be called on the group.
    /// @tparam Args argument types for the functor
    ///
    /// @param id the element id
    /// @param args additional functor arguments
    ///
    /// @return the functor output
    template <typename functor_t, typename link_t, typename... Args>
    DETRAY_HOST_DEVICE typename functor_t::output_type call(
        const link_t link, Args &&...args) const {
        return m_tuple_container.template call<functor_t>(
            value_types::to_index(detail::get<0>(link)), detail::get<1>(link),
            std::forward<Args>(args)...);
    }

    private:
    /// The underlying tuple container implementation
    tuple_type m_tuple_container;
};

/// A stand-alone function to get the vecmem view of the tuple container
///
/// @note the @c view_type typedef will not be available, if one of the element
/// types does not define a vecmem view.
///
/// @return the view on this tuple container
template <typename ID, typename context_t, template <typename...> class tuple_t,
          template <typename...> class container_t, typename... Ts>
inline typename multi_type_store<ID, context_t, tuple_t, container_t,
                                 Ts...>::view_type
get_data(
    multi_type_store<ID, context_t, tuple_t, container_t, Ts...> &container) {
    return get_data(*container.data());
}

}  // namespace detray