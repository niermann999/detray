# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Look for openMP, which is used for the CPU benchmark
find_package(OpenMP)

detray_add_unit_test(image
   "test_image.cpp"
   LINK_LIBRARIES detray::tracer detray::io
)

detray_add_unit_test(render_masks
   "render_masks.cpp"
   LINK_LIBRARIES detray::algebra_vc_aos detray::algebra_vc_soa detray::tracer detray::io
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(
        detray_unit_test_render_masks
        PRIVATE OpenMP::OpenMP_CXX
    )
endif()
