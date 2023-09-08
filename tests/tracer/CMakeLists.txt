# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Look for openMP, which is used for the CPU benchmark
find_package(OpenMP)

detray_add_executable(test_image
   "test_image.cpp"
   LINK_LIBRARIES detray::test detray::io detray::utils)

detray_add_executable(render_masks
   "single_shapes.cpp"
   LINK_LIBRARIES detray_algebra_array detray_algebra_vc_soa detray::test detray::io detray::utils)

if( OpenMP_CXX_FOUND )
  target_link_libraries( detray_render_masks PRIVATE OpenMP::OpenMP_CXX )
endif()
