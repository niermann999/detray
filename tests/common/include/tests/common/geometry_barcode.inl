/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

// Detray include(s)
#include "detray/geometry/barcode.hpp"

using namespace detray;

/// Test retrieval of surface from collection using brute force searching
TEST(geometry, barcode) {

    auto bcd = geometry::barcode{};

    // Check a empty barcode
    EXPECT_EQ(bcd.volume(), (1UL << 12) - 1UL);
    EXPECT_EQ(bcd.id(), static_cast<surface_id>((1UL << 4) - 1UL));
    EXPECT_EQ(bcd.index(), (1UL << 40) - 1UL);
    EXPECT_EQ(bcd.extra(), (1UL << 8) - 1UL);

    bcd.set_volume(2UL)
        .set_id(surface_id::e_passive)
        .set_index(42UL)
        .set_extra(24UL);

    // Check the values after setting them
    EXPECT_EQ(bcd.volume(), 2UL);
    EXPECT_EQ(bcd.id(), surface_id::e_passive);
    EXPECT_EQ(bcd.index(), 42UL);
    EXPECT_EQ(bcd.extra(), 24UL);
}