# Regression tests, one executable per .cpp, registered with CTest.

find_package(GTest REQUIRED)

# pcg-cpp is a header-only PRNG used by several tests. It has no CMake package,
# so locate the header directly.
find_path(PCG_INCLUDE_DIR pcg_random.hpp)
if(NOT PCG_INCLUDE_DIR)
    message(WARNING
        "pcg_random.hpp not found; regression tests that need it will fail to build. "
        "Install pcg-cpp or set PCG_INCLUDE_DIR.")
endif()

file(GLOB test_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach(src IN LISTS test_sources)
    get_filename_component(name "${src}" NAME_WE)
    add_executable(${name} "${src}")
    # GTest::gtest_main pulls in GTest::gtest transitively.
    target_link_libraries(${name} PRIVATE
        geomc::geomc
        GTest::gtest_main
    )
    # Tests include sibling helper headers (e.g. shape_generation.h) directly.
    target_include_directories(${name} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    if(PCG_INCLUDE_DIR)
        target_include_directories(${name} PRIVATE "${PCG_INCLUDE_DIR}")
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(${name} PRIVATE -Wno-character-conversion)
    endif()
    # Some geometry tests rely on FMA/native rounding for their tolerances to
    # hold; the tests are built for the local machine, so tune for it.
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${name} PRIVATE -march=native)
    endif()
    add_test(NAME ${name} COMMAND ${name})
endforeach()
