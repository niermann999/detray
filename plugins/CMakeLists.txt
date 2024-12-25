# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Add all subdirectories.
add_subdirectory(algebra)

if(DETRAY_SVG_DISPLAY)
    add_subdirectory(svgtools)
endif()

if(DETRAY_BUILD_TRACER)
    add_subdirectory(tracer)
endif()
