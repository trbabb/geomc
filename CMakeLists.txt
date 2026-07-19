cmake_minimum_required(VERSION 3.21)

project(geomc
    VERSION 1.0.0
    DESCRIPTION "Header-only C++ library for geometry and linear algebra"
    HOMEPAGE_URL "https://github.com/trbabb/geomc"
    LANGUAGES CXX
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------

# geomc's iterator adaptors are built on Boost.Iterator (header-only).
find_package(Boost REQUIRED CONFIG)

# ---------------------------------------------------------------------------
# The library (header-only INTERFACE target)
# ---------------------------------------------------------------------------

add_library(geomc INTERFACE)
add_library(geomc::geomc ALIAS geomc)

target_compile_features(geomc INTERFACE cxx_std_23)

# Includes are of the form `#include <geomc/linalg/Vec.h>`, so the include
# root is the directory that *contains* the `geomc/` folder.
target_include_directories(geomc INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(geomc INTERFACE Boost::headers)

# ---------------------------------------------------------------------------
# Install rules
# ---------------------------------------------------------------------------

install(TARGETS geomc EXPORT geomcTargets)

# Install the headers, preserving the `geomc/...` directory layout.
install(DIRECTORY geomc
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

set(GEOMC_CMAKE_CONFIG_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/geomc)

install(EXPORT geomcTargets
    FILE geomcTargets.cmake
    NAMESPACE geomc::
    DESTINATION ${GEOMC_CMAKE_CONFIG_DESTINATION}
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/geomcConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/geomcConfig.cmake
    INSTALL_DESTINATION ${GEOMC_CMAKE_CONFIG_DESTINATION}
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/geomcConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
    ARCH_INDEPENDENT
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/geomcConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/geomcConfigVersion.cmake
    DESTINATION ${GEOMC_CMAKE_CONFIG_DESTINATION}
)

# ---------------------------------------------------------------------------
# Tests (only built when geomc is the top-level project, or on request)
# ---------------------------------------------------------------------------

option(GEOMC_BUILD_TESTS "Build the geomc regression tests" ${PROJECT_IS_TOP_LEVEL})

if(GEOMC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(regression)
endif()

# ---------------------------------------------------------------------------
# Documentation
# ---------------------------------------------------------------------------

option(GEOMC_BUILD_DOCS "Build the geomc API documentation (requires Doxygen)" OFF)

if(GEOMC_BUILD_DOCS)
    find_package(Doxygen REQUIRED)
    add_custom_target(docs
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/doc/gen
        COMMAND ${DOXYGEN_EXECUTABLE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
endif()
