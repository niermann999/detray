/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024W CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <limits>
#include <memory>

namespace detray {

/// Simple geometry that consists of a few surfaces
template <concepts::algebra algebra_t, typename mask_t>
struct simple_geometry {
    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;

    using surface_type = detray::surface_descriptor<>;
    using material_t = detray::material<dvalue<algebra_t>>;

    /// Construct from surface data:
    DETRAY_HOST_DEVICE
    simple_geometry(const std::vector<transform3_type> &trf,
                    std::vector<mask_t> &&mask,
                    const std::vector<material_t> &mat)
        : m_trf{std::move(trf)},
          m_mask{std::move(mask)},
          m_material{std::move(mat)} {}

    /// Threadsafe interface
    /// @{
    const std::vector<transform3_type> &transform() const { return m_trf; }
    const std::vector<mask_t> &mask() const { return m_mask; }
    const std::vector<material_t> &material() const { return m_material; }
    /// @}

    /// The surfaces data
    std::vector<transform3_type> m_trf;
    std::vector<mask_t> m_mask;
    std::vector<material_t> m_material;
};

}  // namespace detray
