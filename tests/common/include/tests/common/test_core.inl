/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>
#include <cmath>

#include "detray/core/intersection.hpp"
#include "detray/core/type_registry.hpp"
#include "detray/geometry/surface_base.hpp"
#include "detray/masks/unmasked.hpp"
#include "detray/tools/intersection_kernel.hpp"

/// @note __plugin has to be defined with a preprocessor command

using namespace detray;

using point2 = __plugin::point2<detray::scalar>;

// Three-dimensional definitions
using transform3 = __plugin::transform3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;
using point3 = __plugin::point3<detray::scalar>;
using mask_defs = default_mask_registry<unmasked<>>;
using mask_link_t = typename mask_defs::link_type;

constexpr scalar epsilon = std::numeric_limits<scalar>::epsilon();

// This tests the construction of a surface_base object
TEST(ALGEBRA_PLUGIN, surface_base) {
    // Preparatioon work, create a transform
    vector3 z = vector::normalize(vector3{3., 2., 1.});
    vector3 x = vector::normalize(vector3{2., -3., 0.});
    point3 t{2., 3., 4.};
    transform3 trf(t, z, x);

    mask_link_t mask_id{mask_defs::e_unknown, 0};
    surface_base<mask_defs, intersection_kernel, transform3, dindex> s(
        std::move(trf), std::move(mask_id), -1, false);
}

// This tests the construction of a intresection
TEST(ALGEBRA_PLUGIN, intersection) {

    intersection i0 = {2., point3{0.3, 0.5, 0.7}, point2{0.2, 0.4},
                       intersection_status::e_hit};

    intersection i1 = {1.7, point3{0.2, 0.3, 0.}, point2{0.2, 0.4},
                       intersection_status::e_inside};

    intersection invalid;
    ASSERT_TRUE(invalid.status == intersection_status::e_missed);

    dvector<intersection> intersections = {invalid, i0, i1};
    std::sort(intersections.begin(), intersections.end());

    ASSERT_NEAR(intersections[0].path, 1.7, epsilon);
    ASSERT_NEAR(intersections[1].path, 2, epsilon);
    ASSERT_TRUE(std::isinf(intersections[2].path));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
