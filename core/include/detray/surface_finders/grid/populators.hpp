/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/invalid_values.hpp"

namespace detray {

/// A replace populator that overrides whatever current content is in the bin
/// with a new entry.
///
/// @note entry type and bin content type may not be identicial for all bin
/// types
struct replace {

    /// Sort the entries contained in a bin content when fetched: not needed,
    /// as there is only a single entry
    static constexpr bool do_sort = false;

    /// Replace the bin content with a new entry - forwarding
    ///
    /// @param bin the bin for which to replace the content
    /// @param content new content to be added
    template <typename bin_t, typename content_t>
    DETRAY_HOST_DEVICE void operator()(bin_t &bin, content_t &&entry) const {
        bin.init(std::forward<content_t>(entry));
    }
};

/// A regular attach populator that adds a new entry to a given bin.
///
/// @tparam kSORT sort the entries in the bin
template <bool kSORT = false>
struct attach {

    /// Sort the entries contained in a bin entry when viewed
    static constexpr bool do_sort = kSORT;

    /// Append a new entry to the bin - forwarding
    ///
    /// @param bin the bin for which to replace the content
    /// @param content new content to be added
    template <typename bin_t, typename entry_t>
    DETRAY_HOST_DEVICE void operator()(bin_t &bin, entry_t &&entry) const {
        bin.push_back(std::forward<entry_t>(entry));
    }
};

/// A complete populator that adds elements to a bin content which contains an
/// array of entries, until it is completed - ignored afterwards.
///
/// @tparam kSORT sort the entries in the bin content
template <bool kSORT = false>
struct complete {

    /// Sort the entries contained in a bin content when viewed
    static constexpr bool do_sort = kSORT;

    /// Complete the bin content with a new entry - copy
    ///
    /// @param bin the bin for which to replace the content
    /// @param content new content to be added
    template <typename bin_t, typename entry_t>
    DETRAY_HOST_DEVICE void operator()(bin_t &bin, entry_t &&entry) const {
        for (dindex i{bin.size()}; i < bin.capacity(); ++i) {
            bin.push_back(std::forward<entry_t>(entry));
        }
    }
};

}  // namespace detray
