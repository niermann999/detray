/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/coordinates/cartesian3.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/tracks.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;
using point2 = test::point2;
using point3 = test::point3;
using vector3 = test::vector3;
using transform3 = test::transform3;
using matrix_operator = typename transform3::matrix_actor;

const scalar isclose{1e-5f};

// This test cartesian3 coordinate
GTEST_TEST(detray_coordinates, cartesian3) {

    // Preparation work
    const vector3 z = {0.f, 0.f, 1.f};
    const vector3 x = {0.f, 1.f, 0.f};
    const point3 t = {2.f, 3.f, 4.f};
    const transform3 trf(t, z, x);
    const cartesian3<transform3> c3;
    const point3 global1 = {4.f, 7.f, 5.f};
    const vector3 glob_vec = {4.f, 7.f, 4.f};
    const vector3 mom = {1.f, 2.f, 3.f};

    // Global to local transformation
    const point3 local = c3.global_to_local(trf, global1);

    // Check if the local position is correct
    ASSERT_NEAR(local[0], 4.f, isclose);
    ASSERT_NEAR(local[1], -2.f, isclose);
    ASSERT_NEAR(local[2], 1.f, isclose);

    // Local to global transformation
    const point3 global2 = c3.local_to_global(trf, local);

    // Check if the same global position is obtained
    ASSERT_NEAR(global1[0], global2[0], isclose);
    ASSERT_NEAR(global1[1], global2[1], isclose);
    ASSERT_NEAR(global1[2], global2[2], isclose);

    // Vector to local transformation
    const vector3 loc_vec = c3.vector_to_local(trf, glob_vec);

    // Check if the local vector is correct
    ASSERT_NEAR(loc_vec[0], 7.f, isclose);
    ASSERT_NEAR(loc_vec[1], -4.f, isclose);

    // Vector to global transformation
    const vector3 glob_vec2 = c3.vector_to_global(trf, loc_vec);

    // Check if the same global vector is obtained
    ASSERT_NEAR(glob_vec[0], glob_vec2[0], isclose);
    ASSERT_NEAR(glob_vec[1], glob_vec2[1], isclose);
    ASSERT_NEAR(glob_vec[2], glob_vec2[2], isclose);
}
