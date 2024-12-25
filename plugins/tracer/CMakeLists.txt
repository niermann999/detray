# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building 'detray::tracer' component")

# Set up the core library.
file(
    GLOB _detray_tracer_public_headers
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "include/detray/plugins/tracer/definitions/*.hpp"
    "include/detray/plugins/tracer/renderer/*.hpp"
    "include/detray/plugins/tracer/shaders/detail/*.hpp"
    "include/detray/plugins/tracer/shaders/*.hpp"
    "include/detray/plugins/tracer/texture/*.hpp"
)
detray_add_library( detray_tracer tracer
    ${_detray_tracer_public_headers}
)
target_link_libraries(
    detray_tracer
    INTERFACE vecmem::core detray::core detray::io
)

# Set up the libraries that use specific algebra plugins.
detray_add_library( detray_tracer_array tracer_array )
target_link_libraries(
    detray_tracer_array
    INTERFACE detray::tracer detray::algebra_array
)

if(DETRAY_EIGEN_PLUGIN)
    detray_add_library( detray_tracer_eigen tracer_eigen )
    target_link_libraries(
        detray_tracer_eigen
        INTERFACE detray::tracer detray::algebra_eigen
    )
endif()

if(DETRAY_SMATRIX_PLUGIN)
    detray_add_library( detray_tracer_smatrix tracer_smatrix )
    target_link_libraries(
        detray_tracer_smatrix
        INTERFACE detray::tracer detray::algebra_smatrix
    )
endif()

if(DETRAY_VC_AOS_PLUGIN)
    detray_add_library( detray_tracer_vc_aos tracer_vc_aos )
    target_link_libraries(
        detray_tracer_vc_aos
        INTERFACE detray::tracer detray::algebra_vc_aos
    )
endif()

if(DETRAY_VC_SOA_PLUGIN)
    detray_add_library( detray_tracer_vc_soa tracer_vc_soa )
    target_link_libraries(
        detray_tracer_vc_soa
        INTERFACE detray::tracer detray::algebra_vc_soa
    )
endif()

# Test the public headers of the detray tracer libraries.
if(BUILD_TESTING AND DETRAY_BUILD_TESTING)
    string(
        REPLACE
        "include/"
        ""
        _detray_tracer_public_headers
        "${_detray_tracer_public_headers}"
    )
    detray_test_public_headers(detray_tracer_array
      ${_detray_tracer_public_headers}
    )
    if(DETRAY_EIGEN_PLUGIN)
        detray_test_public_headers(detray_tracer_eigen
         ${_detray_tracer_public_headers}
        )
    endif()
    if(DETRAY_SMATRIX_PLUGIN)
        detray_test_public_headers(detray_tracer_smatrix
         ${_detray_tracer_public_headers}
        )
    endif()
    if(DETRAY_VC_PLUGIN)
        detray_test_public_headers(detray_tracer_vc
         ${_detray_tracer_public_headers}
        )
    endif()
endif()

# Clean up.
unset(_detray_tracer_public_headers)
