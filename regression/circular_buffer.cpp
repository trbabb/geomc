#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Hello

#include <boost/test/unit_test.hpp>
#include <geomc/CircularBuffer.h>

using namespace geom;
using namespace std;


// todo: validate objects with copy/move/resource destruction semantics


template <typename T, index_t N>
void validate_stepfill(CircularBuffer<T, N>& buf, index_t ct) {
    BOOST_CHECK(buf.size() == ct);
    BOOST_CHECK(buf.front() == 10);
    BOOST_CHECK(buf.back() == 10 * ct);
    
    for (index_t i = 0; i < ct; ++i) {
        // did the insert work?
        BOOST_CHECK(buf[i] == (i + 1) * 10);
        // does "wraparound" work?
        BOOST_CHECK(buf[i + buf.size()] == buf[i]);
        // does "negative indexing" work?
        BOOST_CHECK(buf[i - buf.size()] == buf[i]);
    }
}


template <typename T, index_t N>
void exercise(CircularBuffer<T, N>& buf, index_t ct) {
    // forwards insertion:
    for (index_t i = 0; i < ct; ++i) {
        buf.push_back((i + 1) * 10);
    }
    
    validate_stepfill(buf, ct);
    
    buf.clear();
    BOOST_CHECK_EQUAL(buf.size(), 0);
    
    // try backwards insertion:
    for (index_t i = 0; i < ct; ++i) {
        buf.push_front((ct - i) * 10);
    }
    
    validate_stepfill(buf, ct);
    
    // clear the buffer by popping
    index_t removed;
    for (removed = 0; buf.size() > 0; ++removed) {
        int w = buf.pop_front();
        BOOST_CHECK(w == (removed + 1) * 10);
    }
    BOOST_CHECK_EQUAL(buf.size(), 0);
    BOOST_CHECK_EQUAL(removed, ct);
}


BOOST_AUTO_TEST_SUITE(circular_buffer)


BOOST_AUTO_TEST_CASE(create_circular_buffer) {
    CircularBuffer<int, 8> buf;
    
    BOOST_CHECK_EQUAL(buf.size(), 0);
    BOOST_CHECK_EQUAL(buf.capacity(), 8);
    
    buf.push_back(111);
    BOOST_CHECK_EQUAL(buf.size(), 1);
    BOOST_CHECK_EQUAL(buf.front(), 111);
    BOOST_CHECK_EQUAL(buf.back(), 111);
    
    int z = buf.pop_front();
    BOOST_CHECK_EQUAL(z, 111);
    BOOST_CHECK_EQUAL(buf.size(), 0);
    
    buf.push_front(222);
    BOOST_CHECK_EQUAL(buf.size(), 1);
    BOOST_CHECK_EQUAL(buf.front(), 222);
    BOOST_CHECK_EQUAL(buf.back(), 222);
    
    z = buf.pop_back();
    BOOST_CHECK_EQUAL(z, 222);
    BOOST_CHECK_EQUAL(buf.size(), 0);
}


BOOST_AUTO_TEST_CASE(exercise_small_buffer) {
    CircularBuffer<int, 8> buf;
    // try the buffer with a small number of items.
    exercise(buf, 6);
    // try the buffer again, but with the head at the middle of the buffer:
    exercise(buf, 6);
}


BOOST_AUTO_TEST_CASE(exercise_large_buffer) {
    CircularBuffer<int, 8> buf;
    // try the buffer with a large number of items.
    exercise(buf, 19);
    // try again, but with the head in the middle of the buffer:
    exercise(buf, 19);
}


BOOST_AUTO_TEST_SUITE_END()
