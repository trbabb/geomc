Geomc linear algebra and geometry library

Tim Babb  
tr.babb@gmail.com

Usage
=====

Geomc is a header-only C++23 library. Includes are of the form:

    #include <geomc/linalg/Matrix.h>

It depends on [Boost](https://www.boost.org/) (header-only, for its iterator
adaptors).

Documentation
=============

[Geomc topics](http://trbabb.github.io/geomc/html/topics.html)

Installing
==========

Geomc builds and installs with CMake. To install the headers and CMake package
files (to `/usr/local` by default):

    cmake -B build
    cmake --install build

To install to a custom location:

    cmake -B build
    cmake --install build --prefix /path/to/prefix

Because the library is header-only, "installing" simply copies the headers and
the CMake package config; there is nothing to compile.

Using it from CMake
===================

Once installed, consumers can find geomc with `find_package` and link the
imported target, which carries the include paths and the Boost dependency:

    find_package(geomc REQUIRED)

    add_executable(my_app main.cpp)
    target_link_libraries(my_app PRIVATE geomc::geomc)

If geomc was installed to a non-standard prefix, point CMake at it with
`-DCMAKE_PREFIX_PATH=/path/to/prefix`.

You can also consume geomc directly from a source checkout without installing,
via `add_subdirectory(geomc)` or `FetchContent`; the same `geomc::geomc` target
is available.

Tests
=====

The regression tests are built when geomc is the top-level project (or with
`-DGEOMC_BUILD_TESTS=ON`). They require [GoogleTest](https://github.com/google/googletest)
and [pcg-cpp](https://github.com/imneme/pcg-cpp):

    cmake -B build -DGEOMC_BUILD_TESTS=ON
    cmake --build build
    ctest --test-dir build

Building the docs
=================

To build the API documentation (requires Doxygen):

    cmake -B build -DGEOMC_BUILD_DOCS=ON
    cmake --build build --target docs
