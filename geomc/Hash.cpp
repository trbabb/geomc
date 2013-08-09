#include <geomc/Hash.h>

//TODO: this is slower. :(
size_t geom::general_hash(const void *data, size_t bytes) {
    const static size_t k = sizeof(size_t);
    const static size_t m = detail::HashParams<k>::multiplier;
    const static size_t a = detail::HashParams<k>::increment;
    size_t h = a;
    size_t chunks = bytes / sizeof(size_t);
    for (size_t chunk = 0; chunk < chunks; chunk++) {
        h = m * h + ((size_t*)data)[chunk];
    }
    size_t remaining = bytes % sizeof(size_t);
    size_t h_0 = 0;
    const unsigned char *last = ((const unsigned char*) data) + chunks * k;
    for (size_t b = 0; b < remaining; b++, last++) {
        h_0 = h_0 << 8 | (*last);
    }
    h = m * h + h_0;
    return h;
}
