Geomc linear algebra and geometry library

Tim Babb  
tr.babb@gmail.com

Usage
=====

Includes are of the form:

    #include <geomc/linalg/Matrix.h>

Documentation
=============

[Geomc topics](http://trbabb.github.io/geomc/html/topics.html)

Building
========

To build the library, run:

    scons install

To make the documentation, run:

    scons docs

These build options are also available:
    
- Target webassembly:
    `scons --wasm`
- Enable address sanitization for debugging:
    `scons --sanitize`
- Disable optimization:
    `scons noopt=1`
- Enable debug symbols:
    `scons debug=1`
- Rebuild `compile_commands.json` for clang tools:
    `scons compile_commands`
