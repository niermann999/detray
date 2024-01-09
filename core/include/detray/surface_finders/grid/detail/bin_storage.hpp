/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/ranges.hpp"

namespace detray::detail {

template <bool is_owning, typename bin_t, typename containers>
class bin_storage {};

/// @brief bin data state of a grid - owning
///
/// The data that is owned and managed by the grid. Does not contain the
/// data of the axes, as that is managed by the multi-axis type directly.
template <typename bin_t, typename containers>
class bin_storage<true, bin_t, containers>
    : public detray::ranges::view_interface<
          bin_storage<true, bin_t, containers>> {
    public:
    using bin_type = bin_t;
    /// Backend storage type for the grid
    using bin_container_type = typename containers::template vector_type<bin_t>;

    // Vecmem based view type
    using view_type = dvector_view<bin_type>;
    using const_view_type = dvector_view<const bin_type>;

    // Vecmem based buffer type
    using buffer_type = dvector_buffer<bin_type>;

    /// Default constructor
    bin_storage() = default;

    /// Construct containers using a memory resources
    DETRAY_HOST
    bin_storage(vecmem::memory_resource& resource) : m_bin_data(&resource) {}

    /// Construct grid data from containers - move
    DETRAY_HOST_DEVICE
    bin_storage(bin_container_type&& bin_data)
        : m_bin_data(std::move(bin_data)) {}

    /// Construct containers from a vecmem view
    DETRAY_HOST_DEVICE
    bin_storage(
        const dvector_view<typename bin_container_type::value_type>& view)
        : m_bin_data(view) {}

    DETRAY_HOST_DEVICE
    auto begin() -> bin_type* { return m_bin_data.data(); }
    DETRAY_HOST_DEVICE
    auto begin() const -> const bin_type* { return m_bin_data.data(); }
    DETRAY_HOST_DEVICE
    auto end() -> bin_type* { return m_bin_data.data() + m_bin_data.size(); }
    DETRAY_HOST_DEVICE
    auto end() const -> const bin_type* {
        return m_bin_data.data() + m_bin_data.size();
    }

    /// @returns view of the bin storage
    DETRAY_HOST auto get_data() -> view_type {
        return detray::get_data(m_bin_data);
    }

    /// @returns view of the bin storage - const
    DETRAY_HOST auto get_data() const -> const_view_type {
        return detray::get_data(m_bin_data);
    }

    private:
    /// Contains all bin data
    bin_container_type m_bin_data{};
};

/// @brief bin data state of a grid - non-owning
///
/// The grid holds a view onto a global collection of data. This type is NOT
/// used to move @c bin_storage from host to device.
/// @todo Use spans
template <typename bin_t, typename containers>
class bin_storage<false, bin_t, containers>
    : public detray::ranges::subrange<
          typename containers::template vector_type<bin_t>> {
    public:
    using bin_type = bin_t;
    /// Backend storage type for the grid
    using bin_container_type = typename containers::template vector_type<bin_t>;

    // Vecmem based view type (needed only for compilation)
    using view_type = dvector_view<bin_type>;
    using const_view_type = dvector_view<const bin_type>;

    // Vecmem based buffer type (needed only for compilation)
    using buffer_type = dvector_buffer<bin_type>;

    using base_type = detray::ranges::subrange<bin_container_type>;
    using base_type::base_type;
};

}  // namespace detray::detail
