/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/surface_finders/grid/axis.hpp"
#include "detray/surface_finders/grid/detail/bin_view.hpp"
#include "detray/surface_finders/grid/detail/grid_helpers.hpp"
#include "detray/surface_finders/grid/populator.hpp"
#include "detray/surface_finders/grid/serializer.hpp"
#include "detray/utils/ranges.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <array>
#include <cstddef>
#include <type_traits>

namespace detray {

/// @brief An N-dimensional grid for object storage.
///
/// @tparam multi_axis_t the types of grid axes
/// @tparam value_t type of data to be stored in the grid's bins.
/// @tparam serializer_t how to serialize axis-local bin indices into global bin
///                      indices in the grid backend storage and vice versa.
/// @tparam populator_impl_t is a prescription what to do when a bin gets
///                          populated or read. It thus also defines the backend
///                          storage type/value type.
template <typename multi_axis_t, typename value_t,
          template <std::size_t> class serializer_t, typename populator_impl_t>
class grid {

    public:
    /// Grid dimension
    static constexpr unsigned int Dim = multi_axis_t::Dim;
    static constexpr bool is_owning = multi_axis_t::is_owning;

    /// Single value in a bin entry
    using value_type = value_t;
    template <std::size_t DIM>
    using serializer_type = serializer_t<DIM>;
    /// Interface to the populator (determines the bin content and value type).
    using populator_impl = populator_impl_t;
    using bin_type =
        typename populator<populator_impl>::template bin_type<value_type>;

    /// The type of the multi-axis is tied to the type of the grid: a non-
    /// owning grid holds a non-owning multi-axis member.
    using axes_type = multi_axis_t;
    using glob_bin_index = dindex;
    using loc_bin_index = typename axes_type::loc_bin_index;
    using local_frame_type = typename axes_type::local_frame_type;
    using point_type = typename axes_type::point_type;
    using scalar_type = typename axes_type::scalar_type;

    /// How to define a neighborhood for this grid
    template <typename neighbor_t>
    using neighborhood_type = std::array<neighbor_t, Dim>;

    /// Backend storage type for the grid
    using bin_storage_type =
        typename axes_type::container_types::template vector_type<bin_type>;

    /// Vecmem based grid view type
    using view_type =
        dmulti_view<dvector_view<bin_type>, typename axes_type::view_type>;
    /// Vecmem based grid view type - const
    using const_view_type = dmulti_view<dvector_view<const bin_type>,
                                        typename axes_type::const_view_type>;
    /// Vecmem based buffer type
    using buffer_type = dmulti_buffer<dvector_buffer<bin_type>,
                                      typename axes_type::buffer_type>;

    /// Grid backend can be owning (single grid) or non-owning (grid collection)
    using storage_type =
        std::conditional_t<is_owning, detail::grid_data<bin_storage_type>,
                           detail::grid_view<bin_storage_type>>;
    /// Find the corresponding (non-)owning grid type
    template <bool owning>
    using type = grid<typename multi_axis_t::template type<owning>, value_t,
                      serializer_t, populator_impl_t>;

    /// Make grid default constructible: Empty grid with empty axis
    grid() = default;

    /// Create empty grid with empty axes from specific vecmem memory resource
    DETRAY_HOST
    explicit grid(vecmem::memory_resource &resource)
        : m_data(resource), m_axes(resource) {}

    /// Create grid with well defined @param axes and @param bins_data - move
    DETRAY_HOST_DEVICE
    grid(bin_storage_type &&bin_data, axes_type &&axes)
        : m_data(std::move(bin_data)), m_axes(std::move(axes)) {}

    /// Create grid from container pointers - non-owning (both grid and axes)
    DETRAY_HOST_DEVICE
    grid(const bin_storage_type *bin_data_ptr, const axes_type &axes,
         const dindex offset = 0)
        : m_data(bin_data_ptr, offset), m_axes(axes) {}

    /// Create grid from container pointers - non-owning (both grid and axes)
    DETRAY_HOST_DEVICE
    grid(bin_storage_type *bin_data_ptr, axes_type &axes,
         const dindex offset = 0u)
        : m_data(bin_data_ptr, offset), m_axes(axes) {}

    /// Create grid from container pointers - non-owning (both grid and axes)
    // TODO: Remove
    DETRAY_HOST_DEVICE
    grid(const bin_storage_type *bin_data_ptr, axes_type &&axes,
         const dindex offset = 0)
        : m_data(const_cast<bin_storage_type *>(bin_data_ptr), offset),
          m_axes(axes) {}

    /// Create grid from container pointers - non-owning (both grid and axes)
    DETRAY_HOST_DEVICE
    grid(bin_storage_type *bin_data_ptr, axes_type &&axes,
         const dindex offset = 0u)
        : m_data(bin_data_ptr, offset), m_axes(axes) {}

    /// Device-side construction from a vecmem based view type
    template <typename grid_view_t,
              typename std::enable_if_t<detail::is_device_view_v<grid_view_t>,
                                        bool> = true>
    DETRAY_HOST_DEVICE grid(grid_view_t &view)
        : m_data(detray::detail::get<0>(view.m_view)),
          m_axes(detray::detail::get<1>(view.m_view)) {}

    /// @returns the underlying bin storage - const
    DETRAY_HOST_DEVICE
    auto bin_data() const -> const bin_storage_type & {
        return *(m_data.bin_data());
    }

    /// @returns the underlying bin storage
    DETRAY_HOST_DEVICE
    auto bin_data() -> bin_storage_type & { return *(m_data.bin_data()); }

    /// @returns the offset into a global collection - const
    DETRAY_HOST_DEVICE
    auto offset() const -> dindex { return m_data.offset(); }

    /// @returns the multi-axis used by the grid - const
    DETRAY_HOST_DEVICE
    auto axes() const -> const axes_type & { return m_axes; }

    /// @returns the full range of bins - const
    DETRAY_HOST_DEVICE
    auto bins() const {
        dindex first_bin{offset()};
        return detray::ranges::subrange(
            bin_data(), dindex_range{first_bin, first_bin + nbins()});
    }

    /// @returns the full range of bins
    DETRAY_HOST_DEVICE
    auto bins() {
        dindex first_bin{offset()};
        return detray::ranges::subrange(
            bin_data(), dindex_range{first_bin, first_bin + nbins()});
    }

    /// @returns the grid local coordinate system
    DETRAY_HOST_DEVICE
    static constexpr auto local_frame() -> local_frame_type { return {}; }

    /// @returns an axis object, corresponding to the index.
    template <std::size_t index>
    DETRAY_HOST_DEVICE inline constexpr auto get_axis() const {
        return m_axes.template get_axis<index>();
    }

    /// @returns an axis object, corresponding to the label.
    template <n_axis::label L>
    DETRAY_HOST_DEVICE inline constexpr auto get_axis() const {
        return m_axes.template get_axis<L>();
    }

    /// @returns an axis object of the given type.
    template <typename axis_t>
    DETRAY_HOST_DEVICE inline constexpr axis_t get_axis() const {
        return m_axes.template get_axis<axis_t>();
    }

    /// @returns the total number of bins in the grid
    DETRAY_HOST_DEVICE inline constexpr auto nbins() const -> dindex {
        const auto n_bins_per_axis = m_axes.nbins();
        dindex n_bins{1u};
        for (dindex i = 0u; i < Dim; ++i) {
            n_bins *= n_bins_per_axis[i];
        }
        return n_bins;
    }

    /// @returns the total number of values in the grid
    /// @note this has to query every bin for the number of elements
    DETRAY_HOST_DEVICE inline constexpr auto size() const -> dindex {
        return static_cast<dindex>(all().size());
    }

    /// Transform a point in global cartesian coordinates to local coordinates
    ///
    /// @param trf the placement transform of the grid (e.g. from a volume or
    ///            a surface).
    /// @param p   the point in global coordinates
    ///
    /// @returns a point in the coordinate system that is spanned by the grid's
    /// axes.
    template <typename transform_t, typename point3_t, typename vector3_t>
    DETRAY_HOST_DEVICE point_type global_to_bound(const transform_t &trf,
                                                  const point3_t &p,
                                                  const vector3_t &d) const {
        return local_frame().global_to_bound(trf, p, d);
    }

    /// @returns the iterable view of the bin content
    /// @{
    /// @param indices the single indices corresponding to a multi_bin
    template <typename... I, std::enable_if_t<sizeof...(I) == Dim, bool> = true>
    DETRAY_HOST_DEVICE auto bin(I... indices) const {
        return bin(loc_bin_index{indices...});
    }

    /// @param mbin the multi-index of bins over all axes - const
    DETRAY_HOST_DEVICE
    auto bin(const loc_bin_index &mbin) const {
        return bin(m_serializer(m_axes, mbin));
    }

    /// @param mbin the multi-index of bins over all axes
    DETRAY_HOST_DEVICE
    auto bin(const loc_bin_index &mbin) {
        return bin(m_serializer(m_axes, mbin));
    }

    /// @param gbin the global bin index - const
    DETRAY_HOST_DEVICE
    auto bin(const glob_bin_index gbin) const {
        return m_populator.view(bin_data(), gbin + offset());
    }

    /// @param gbin the global bin index
    DETRAY_HOST_DEVICE
    auto bin(const glob_bin_index gbin) {
        return m_populator.view(bin_data(), gbin + offset());
    }
    /// @}

    /// Access a single entry in a bin from the global bin index, as well as
    /// the index of the entry in the bin.
    ///
    /// @param idx the index of a specific grid entry
    ///
    /// @returns a single bin entry.
    /// @{
    DETRAY_HOST_DEVICE
    auto at(const glob_bin_index gbin, const dindex pos) {
        return bin(gbin)[pos];
    }
    DETRAY_HOST_DEVICE
    auto at(const glob_bin_index gbin, const dindex pos) const {
        return bin(gbin)[pos];
    }
    DETRAY_HOST_DEVICE
    auto at(const loc_bin_index &mbin, const dindex pos) {
        return bin(mbin)[pos];
    }
    DETRAY_HOST_DEVICE
    auto at(const loc_bin_index &mbin, const dindex pos) const {
        return bin(mbin)[pos];
    }
    /// @}

    /// @returns a view over the flatened bin content by joining the bin ranges
    DETRAY_HOST_DEVICE auto all() { return detray::views::join(bins()); }

    /// @returns a view over the flatened bin content by joining the bin ranges
    DETRAY_HOST_DEVICE auto all() const { return detray::views::join(bins()); }

    /// Interface for the navigator
    template <typename detector_t, typename track_t>
    DETRAY_HOST_DEVICE auto search(
        const detector_t & /*det*/,
        const typename detector_t::volume_type & /*volume*/,
        const track_t & /*track*/) const {
        // Track position in grid coordinates
        /*const auto &trf = det.transform_store()[volume.transform()];
        const auto loc_pos = global_to_local(trf, track.pos(), track.dir());
        // Grid lookup
        return search(loc_pos);*/
        // @todo: Implement local neighborhood lookup
        return all();
    }

    /// Find the value of a single bin - const
    ///
    /// @param p is point in the local (bound) frame
    ///
    /// @return the iterable view of the bin content
    DETRAY_HOST_DEVICE auto search(const point_type &p) const {
        return bin(m_axes.bins(p));
    }

    /// Find the value of a single bin
    ///
    /// @param p is point in the local (bound) frame
    ///
    /// @return the iterable view of the bin content
    DETRAY_HOST_DEVICE auto search(const point_type &p) {
        return bin(m_axes.bins(p));
    }

    /// @brief Return a neighborhood of values from the grid
    ///
    /// The lookup is done with a search window around the bin
    ///
    /// @param p is point in the local frame
    /// @param win_size size of the binned/scalar search window
    ///
    /// @return the sequence of values
    template <typename neighbor_t>
    DETRAY_HOST_DEVICE auto search(
        const point_type &p, const std::array<neighbor_t, 2> &win_size) const {

        // Return iterable over bins in the search window
        auto search_window = axes().bin_ranges(p, win_size);
        auto search_area = n_axis::detail::bin_view(*this, search_window);

        // Join the respective bins to a single iteration
        return detray::views::join(std::move(search_area));
    }

    /// Poupulate a bin with a single one of its corresponding values @param v
    /// @{
    /// @param mbin the multi bin index to be populated
    DETRAY_HOST_DEVICE auto populate(const loc_bin_index mbin,
                                     const value_type &v)
        -> n_axis::multi_bin<Dim> {
        populate(m_serializer(m_axes, mbin), v);
        return mbin;
    }

    /// @param mbin the multi bin index to be populated
    DETRAY_HOST_DEVICE
    auto populate(const loc_bin_index mbin, value_type &&v) -> loc_bin_index {
        populate(m_serializer(m_axes, mbin), std::move(v));
        return mbin;
    }

    /// @param gbin the global bin index to be populated
    DETRAY_HOST_DEVICE
    auto populate(const glob_bin_index gbin, const value_type &v)
        -> loc_bin_index {
        m_populator(*(m_data.bin_data()), gbin + m_data.offset(), v);
        return m_serializer(m_axes, gbin);
    }

    /// @param gbin the global bin index to be populated
    DETRAY_HOST_DEVICE
    auto populate(const glob_bin_index gbin, value_type &&v) -> loc_bin_index {
        m_populator(*(m_data.bin_data()), gbin + m_data.offset(), std::move(v));
        return m_serializer(m_axes, gbin);
    }

    /// @param p the point in local coordinates that defines the bin to be
    ///          populated
    DETRAY_HOST_DEVICE
    auto populate(const point_type &p, value_type &&v) -> loc_bin_index {
        auto mbin = m_axes.bins(p);
        populate(m_serializer(m_axes, mbin), std::move(v));
        return mbin;
    }

    /// @param p the point in local coordinates that defines the bin to be
    ///          populated
    DETRAY_HOST_DEVICE
    auto populate(const point_type &p, const value_type &v) -> loc_bin_index {
        auto mbin = m_axes.bins(p);
        populate(m_serializer(m_axes, mbin), v);
        return mbin;
    }
    /// @}

    /// @return the maximum number of surface candidates during a
    /// neighborhood lookup
    DETRAY_HOST_DEVICE constexpr auto n_max_candidates() const -> unsigned int {
        // @todo: Hotfix for the toy geometry
        return 20u;
    }

    /// @returns an instance of the grid serializer
    static constexpr auto serializer() -> serializer_t<Dim> { return {}; }

    /// @returns a local multi-bin index from a global bin index @param gid
    constexpr auto deserialize(const glob_bin_index gid) const
        -> loc_bin_index {
        return m_serializer(axes(), gid);
    }

    /// @returns a global bin index from a local bin index @param mbin
    constexpr auto serialize(const loc_bin_index &mbin) const
        -> glob_bin_index {
        return m_serializer(axes(), mbin);
    }

    /// @returns view of a grid, including the grids multi_axis. Also valid if
    /// the value type of the grid is cv qualified (then value_t propagates
    /// quialifiers) - non-const
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST auto get_data() -> view_type {
        return view_type{detray::get_data(bin_data()),
                         detray::get_data(m_axes)};
    }

    /// @returns view of a grid, including the grids multi_axis. Also valid if
    /// the value type of the grid is cv qualified (then value_t propagates
    /// quialifiers) - const
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST auto get_data() const -> const_view_type {
        return const_view_type{detray::get_data(bin_data()),
                               detray::get_data(axes())};
    }

    private:
    /// Struct that contains the grid's data state
    storage_type m_data{};
    /// The axes of the grid
    axes_type m_axes{};
    /// How to write and fetch bin values from the backend storage
    populator<populator_impl> m_populator{};
    /// Serialization/Deserialization of multi-bin indices to the global bin
    /// data storage
    serializer_t<Dim> m_serializer{};
};

}  // namespace detray
