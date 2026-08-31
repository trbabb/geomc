#define TEST_MODULE_NAME Rect

#include <cstdint>
#include <gtest/gtest.h>
#include <geomc/shape/Rect.h>

using namespace geom;

// These tests focus on unsigned-coordinate Rects, where quantities produced by
// differencing two positions (extents, displacements, distances) must be
// represented in a signed type (`diff_t`) rather than wrapping around. For
// signed and floating point coordinate types `diff_t == T`, so behavior there
// is unchanged and is covered by the general shape tests.

typedef Rect<uint32_t,2> URect2;
typedef Vec<uint32_t,2>  UVec2;
typedef Vec<int32_t,2>   SVec2;

// The signed companion type must be the signed integer of the same width for
// unsigned coordinates, and identical to `T` otherwise.
TEST(TEST_MODULE_NAME, diff_type_selection) {
    EXPECT_TRUE((std::is_same_v<Rect<uint32_t,3>::diff_t, int32_t>));
    EXPECT_TRUE((std::is_same_v<Rect<uint64_t,3>::diff_t, int64_t>));
    EXPECT_TRUE((std::is_same_v<Rect<int32_t, 3>::diff_t, int32_t>));
    EXPECT_TRUE((std::is_same_v<Rect<float,  3>::diff_t, float>));
    EXPECT_TRUE((std::is_same_v<Rect<double, 3>::diff_t, double>));
}

// Inclusive integer extents are `hi - lo + 1`, and are signed so that an empty
// Rect can report a negative extent instead of a huge wrapped one.
TEST(TEST_MODULE_NAME, dimensions_unsigned) {
    URect2 r {UVec2(10, 20), UVec2(40, 25)};
    SVec2 dim = r.dimensions();
    EXPECT_EQ(dim.x, 31);
    EXPECT_EQ(dim.y,  6);

    // an empty Rect reports negative extents rather than wrapping
    URect2 empty {UVec2(40, 25), UVec2(10, 20)};
    SVec2 edim = empty.dimensions();
    EXPECT_LT(edim.x, 0);
    EXPECT_LT(edim.y, 0);
}

// center() must not overflow, even for a Rect spanning nearly the whole range.
TEST(TEST_MODULE_NAME, center_no_overflow) {
    URect2 r {UVec2(10, 20), UVec2(40, 25)};
    UVec2 c = r.center();
    EXPECT_EQ(c.x, 25u); // 10 + (40 - 10 + 1)/2
    EXPECT_EQ(c.y, 23u); // 20 + (25 - 20 + 1)/2

    // (lo + hi) would overflow uint32 here; lo + (hi - lo)/2 does not.
    Rect<uint32_t,1> big {100u, 0xFFFFFFF0u};
    uint32_t expected = 100u + (0xFFFFFFF0u - 100u + 1u) / 2u;
    EXPECT_EQ(big.center(), expected);
}

// Translation by a signed displacement, including the negative direction that
// an unsigned coordinate cannot express on its own.
TEST(TEST_MODULE_NAME, translate_negative) {
    URect2 r {UVec2(10, 20), UVec2(40, 25)};
    URect2 expected {UVec2(5, 17), UVec2(35, 22)};

    EXPECT_EQ(r - SVec2(5, 3),      expected);
    EXPECT_EQ(r + SVec2(-5, -3),    expected);
    EXPECT_EQ(SVec2(-5, -3) + r,    expected); // free operator form
    URect2 s = r;
    s -= SVec2(5, 3);
    EXPECT_EQ(s, expected);
}

// Negative dilation erodes an unsigned Rect.
TEST(TEST_MODULE_NAME, erode) {
    URect2 r {UVec2(10, 20), UVec2(40, 25)};
    URect2 eroded = r.dilated(SVec2(-2, -1));
    EXPECT_EQ(eroded.lo, UVec2(12, 21));
    EXPECT_EQ(eroded.hi, UVec2(38, 24));
}

// Squared distance must be computed in signed space; a naive unsigned
// difference would wrap and produce an astronomically large result.
TEST(TEST_MODULE_NAME, dist2_unsigned) {
    URect2 r {UVec2(10, 20), UVec2(40, 25)};
    // exterior point below-left; nearest contained point is lo = (10, 20)
    // displacement (-6, -3) -> 36 + 9 = 45
    EXPECT_EQ(r.dist2(UVec2(4, 17)), 45);
    // interior point
    EXPECT_EQ(r.dist2(UVec2(25, 23)), 0);
}

// The contact vector is a displacement, so its components are signed and may be
// negative even for unsigned coordinates.
TEST(TEST_MODULE_NAME, contact_vector_signed) {
    URect2 r {UVec2(10, 20), UVec2(40, 25)};
    // b lies to the left of r: disjoint on x, overlapping on y
    URect2 b {UVec2(0, 20), UVec2(5, 25)};
    std::pair<Vec<int32_t,2>,bool> result = r.contact_vector(b);
    EXPECT_FALSE(result.second); // not overlapping
    // to bring r into contact with b it must move left: negative x displacement
    EXPECT_LT(result.first.x, 0);
}

// Volume clamps empty (negative) extents to zero without a wrapped extent
// inflating the result.
TEST(TEST_MODULE_NAME, measure_interior_unsigned) {
    URect2 r {UVec2(10, 20), UVec2(40, 25)};
    EXPECT_EQ(r.measure_interior(), 31u * 6u);

    URect2 empty {UVec2(40, 25), UVec2(10, 20)};
    EXPECT_EQ(empty.measure_interior(), 0u);
}

// set_dimensions resizes about the center and honors the requested (signed)
// lengths exactly.
TEST(TEST_MODULE_NAME, set_dimensions_unsigned) {
    URect2 r {UVec2(10, 20), UVec2(40, 25)};
    UVec2 pre_center = r.center();
    r.set_dimensions(SVec2(10, 4));
    SVec2 dim = r.dimensions();
    EXPECT_EQ(dim.x, 10);
    EXPECT_EQ(dim.y,  4);
    EXPECT_EQ(r.center(), pre_center); // center preserved
}
