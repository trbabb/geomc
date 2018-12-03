#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Hello

// #include <iostream>

#include <boost/test/unit_test.hpp>
#include <geomc/CircularBuffer.h>

using namespace geom;
using namespace std;


// todo: validate item objects with copy/move/resource destruction semantics
// todo: validate that this class itself is a good citizen when it comes to 
//       destruction/move/copy.


template <typename T, index_t N>
void validate_stepfill(CircularBuffer<T, N>& buf, index_t ct) {
    BOOST_CHECK(buf.size() == ct);
    BOOST_CHECK(buf.front() == 10);
    BOOST_CHECK(buf.back() == 10 * ct);
    
    for (index_t i = 0; i < ct; ++i) {
        // did the insert work?
        BOOST_CHECK_EQUAL(buf[i], (i + 1) * 10);
        // does "wraparound" work?
        BOOST_CHECK_EQUAL(buf[i + buf.size()], buf[i]);
        // does "negative indexing" work?
        BOOST_CHECK_EQUAL(buf[i - buf.size()], buf[i]);
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

template <typename T, index_t N>
void exercise_copy_and_move(index_t ct) {
    CircularBuffer<T, N> buf0;
    CircularBuffer<T, N> buf1;
    
    // stepfill
    for (index_t i = 0; i < ct; ++i) {
        buf0.push_back((i + 1) * 10);
    }
    validate_stepfill(buf0, ct);
    
    // copy-assign
    buf1 = buf0;
    validate_stepfill(buf1, ct);
    
    // copy construct
    CircularBuffer<T, N> buf2(buf0);
    validate_stepfill(buf2, ct);
    buf2.clear();
    
    // move assign
    buf2 = std::move(buf1);
    validate_stepfill(buf2, ct);
    BOOST_CHECK_EQUAL(buf1.size(), 0);
    exercise(buf1, N - 1);
    exercise(buf1, N * 2);
    buf2.clear();
    exercise(buf2, N - 1);
    exercise(buf2, N * 2);
    
    // move construct
    CircularBuffer<T, N> buf3(std::move(buf0));
    validate_stepfill(buf3, ct);
    BOOST_CHECK_EQUAL(buf0.size(), 0);
    exercise(buf0, N - 1);
    exercise(buf0, N * 2);
    buf3.clear();
    exercise(buf3, N - 1);
    exercise(buf3, N * 2);
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


BOOST_AUTO_TEST_CASE(exercise_small_copy_and_move) {
    exercise_copy_and_move<int, 8>(7);
}

BOOST_AUTO_TEST_CASE(exercise_large_copy_and_move) {
    exercise_copy_and_move<int, 8>(15);
}


BOOST_AUTO_TEST_SUITE_END()
