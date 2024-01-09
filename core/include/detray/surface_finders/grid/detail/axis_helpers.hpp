/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray::n_axis {

/// @brief Multi-bin: contains bin indices from multiple axes
template <std::size_t DIM, typename index_t = dindex>
struct multi_bin : public dmulti_index<index_t, DIM> {};

/// @brief Helper to tie two bin indices to a range.
/// @note Cannot use dindex_range for signed integer bin indices.
using bin_range = std::array<int, 2>;

/// @brief Multi-bin-range: contains bin index ranges from multiple axes
template <std::size_t DIM>
struct multi_bin_range : public dmulti_index<bin_range, DIM> {};

namespace detail {

// TODO: Replace by iterator based approach, along the lines of std::ranges

/// @brief data state of a multi-axis - owning
///
/// The multi-axis is owning the axis and binning data: Directly holds the
/// vector containers for the axis data and bin edges.
template <typename containers = host_container_types,
          typename scalar_t = scalar>
struct multi_axis_data {

    // Extract container types
    template <typename T>
    using vector_type = typename containers::template vector_type<T>;

    /// Data that the axes keep: index ranges in the edges container
    vector_type<dindex_range> m_edge_offsets{};
    /// Contains all bin edges for all axes
    vector_type<scalar_t> m_edges{};

    /// Default constructor (empty data)
    multi_axis_data() = default;

    /// Construct containers using a specific memory resources
    DETRAY_HOST
    multi_axis_data(vecmem::memory_resource &resource)
        : m_edge_offsets(&resource), m_edges(&resource) {}

    /// Construct from containers - move
    DETRAY_HOST_DEVICE
    multi_axis_data(vector_type<dindex_range> &&edge_offsets,
                    vector_type<scalar_t> &&edges)
        : m_edge_offsets(std::move(edge_offsets)), m_edges(std::move(edges)) {}

    /// Construct containers from vecmem based view types
    ///
    /// @param axes_view vecmem view on the edge offsets
    /// @param edges_view vecmem view on the bin edges
    DETRAY_HOST_DEVICE
    multi_axis_data(const dvector_view<dindex_range> &axes_view,
                    const dvector_view<scalar_t> &edges_view)
        : m_edge_offsets(axes_view), m_edges(edges_view) {}

    /// @returns pointer to all of the edge offsets
    DETRAY_HOST_DEVICE
    auto edge_offsets() const -> const vector_type<dindex_range> * {
        return &m_edge_offsets;
    }

    /// @returns pointer to all of the edge offsets
    // TODO: Don't do...
    DETRAY_HOST_DEVICE
    auto edge_offsets() -> vector_type<dindex_range> * {
        return &m_edge_offsets;
    }

    /// @returns pointer to the data for one particular @param i axis
    DETRAY_HOST_DEVICE
    auto edge_offsets(const unsigned int i) const -> const dindex_range * {
        return &(m_edge_offsets[i]);
    }

    /// @returns pointer to the entire bin edges collection for every axis
    DETRAY_HOST_DEVICE
    auto edges() const -> const vector_type<scalar_t> * { return &m_edges; }

    /// @returns pointer to the entire bin edges collection for every axis
    DETRAY_HOST_DEVICE
    auto edges() -> vector_type<scalar_t> * { return &m_edges; }

    /// @returns offset of this multi-axis in the global data container (not
    /// needed when data is owned)
    DETRAY_HOST_DEVICE
    constexpr auto offset() const -> dindex { return 0; }
};

/// @brief data state of a multi-axis - non-owning (i.e. container view)
///
/// The multi-axis obtains a view of a global collection of data
/// which it uses to provide the bin index lookup as a pure service to the
/// grid. This type is NOT used to move the @c multi_axis_data between host
/// and device when the data is owned by the multi-axis.
/// @todo Use spans
template <typename containers = host_container_types,
          typename scalar_t = scalar>
struct multi_axis_view {

    // Extract container types
    template <typename T>
    using vector_type = typename containers::template vector_type<T>;

    /// Data that the axes keep: index ranges in the edges container
    const vector_type<dindex_range> *m_edge_offsets{nullptr};
    /// Contains all bin edges for all the axes
    const vector_type<scalar_t> *m_edges{nullptr};
    /// Offset for this multi-axis into the global container
    unsigned int m_offset{0};

    /// Default constructor (no concrete memory access)
    multi_axis_view() = default;

    /// Construct from explicitly given pointers to global containers
    ///
    /// @param edge_offsets_ptr const pointer to edge offsets
    /// @param edges_ptr const pointer to bin edges
    /// @param offset in the edge offsets and bin edges collections pointed to
    DETRAY_HOST_DEVICE
    multi_axis_view(const vector_type<dindex_range> *edge_offsets_ptr,
                    const vector_type<scalar_t> *edges_ptr,
                    const unsigned int offset)
        : m_edge_offsets(edge_offsets_ptr),
          m_edges(edges_ptr),
          m_offset{offset} {}

    /// @returns pointer to all(!) of the edge offsets - const
    DETRAY_HOST_DEVICE
    auto edge_offsets() const -> const vector_type<dindex_range> * {
        return m_edge_offsets;
    }

    /// @returns pointer to all(!) of the edge offsets - non-const
    // TODO: Don't do...
    DETRAY_HOST_DEVICE
    auto edge_offsets() -> vector_type<dindex_range> * {
        return m_edge_offsets;
    }

    /// @returns pointer to the data for one particular @param i axis
    DETRAY_HOST_DEVICE
    auto edge_offsets(const unsigned i) const -> const dindex_range * {
        return &((*m_edge_offsets)[m_offset + i]);
    }

    /// @returns pointer to the entire(!) bin edges collection for every axis
    DETRAY_HOST_DEVICE
    auto edges() const -> const vector_type<scalar_t> * { return m_edges; }

    /// @returns pointer to the entire(!) bin edges collection for every axis
    DETRAY_HOST_DEVICE
    auto edges() -> vector_type<scalar_t> * { return m_edges; }

    /// @returns offset of this multi-axis in the global data container
    DETRAY_HOST_DEVICE
    auto offset() const -> dindex { return m_offset; }
};

}  // namespace detail

}  // namespace detray::n_axis
