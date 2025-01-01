#define TEST_MODULE_NAME CircularBuffer

// #include <iostream>

#include <gtest/gtest.h>
#include <geomc/Deque.h>

using namespace geom;
using namespace std;


// todo: validate item objects with copy/move/resource destruction semantics


template <typename T, index_t N>
void validate_stepfill(Deque<T, N>& buf, index_t ct) {
    EXPECT_EQ(buf.size(), ct);
    EXPECT_EQ(buf.front(), 10);
    EXPECT_EQ(buf.back(), 10 * ct);
    
    for (index_t i = 0; i < ct; ++i) {
        // did the insert work?
        EXPECT_EQ(buf[i], (i + 1) * 10);
        // does "wraparound" work?
        EXPECT_EQ(buf[i + buf.size()], buf[i]);
        // does "negative indexing" work?
        EXPECT_EQ(buf[i - buf.size()], buf[i]);
    }
}


template <typename T, index_t N>
void exercise(Deque<T, N>& buf, index_t ct) {
    // forwards insertion:
    for (index_t i = 0; i < ct; ++i) {
        buf.push_back((i + 1) * 10);
    }
    
    validate_stepfill(buf, ct);
    
    buf.clear();
    EXPECT_EQ(buf.size(), 0);
    
    // try backwards insertion:
    for (index_t i = 0; i < ct; ++i) {
        buf.push_front((ct - i) * 10);
    }
    
    validate_stepfill(buf, ct);
    
    // clear the buffer by popping
    index_t removed;
    for (removed = 0; buf.size() > 0; ++removed) {
        std::optional<int> w = buf.pop_front();
        EXPECT_TRUE(w.has_value());
        EXPECT_TRUE(w == (removed + 1) * 10);
    }
    EXPECT_EQ(buf.size(), 0);
    EXPECT_EQ(removed, ct);
}


template <typename T, index_t N>
void exercise_copy_and_move(index_t ct) {
    Deque<T, N> buf0;
    Deque<T, N> buf1;
    
    // stepfill
    for (index_t i = 0; i < ct; ++i) {
        buf0.push_back((i + 1) * 10);
    }
    validate_stepfill(buf0, ct);
    
    // manually empty it
    while (buf0.size() > 0) { buf0.pop_front(); }
    // stepfill again, but with the head somewhere in the middle
    for (index_t i = 0; i < ct; ++i) {
        buf0.push_back((i + 1) * 10);
    }
    validate_stepfill(buf0, ct);
    
    // copy-assign
    buf1 = buf0;
    validate_stepfill(buf1, ct);
    
    // copy construct
    Deque<T, N> buf2(buf0);
    validate_stepfill(buf2, ct);
    buf2.clear();
    
    // move assign
    buf2 = std::move(buf1);
    validate_stepfill(buf2, ct);
    EXPECT_EQ(buf1.size(), 0);
    exercise(buf1, N - 1);
    exercise(buf1, N * 2);
    buf2.clear();
    exercise(buf2, N - 1);
    exercise(buf2, N * 2);
    
    // move construct
    Deque<T, N> buf3(std::move(buf0));
    validate_stepfill(buf3, ct);
    EXPECT_EQ(buf0.size(), 0);
    exercise(buf0, N - 1);
    exercise(buf0, N * 2);
    buf3.clear();
    exercise(buf3, N - 1);
    exercise(buf3, N * 2);
}

TEST(TEST_MODULE_NAME, create_circular_buffer) {
    Deque<int, 8> buf;
    
    EXPECT_EQ(buf.size(), 0);
    EXPECT_EQ(buf.capacity(), 8);
    
    std::optional<int> z = buf.pop_front();
    EXPECT_TRUE(not z.has_value());
    z = buf.pop_back();
    EXPECT_TRUE(not z.has_value());
    
    buf.push_back(111);
    EXPECT_EQ(buf.size(), 1);
    EXPECT_EQ(buf.front(), 111);
    EXPECT_EQ(buf.back(), 111);
    
    z = buf.pop_front();
    EXPECT_TRUE(z.has_value());
    EXPECT_EQ(*z, 111);
    EXPECT_EQ(buf.size(), 0);
    
    buf.push_front(222);
    EXPECT_EQ(buf.size(), 1);
    EXPECT_EQ(buf.front(), 222);
    EXPECT_EQ(buf.back(), 222);
    
    z = buf.pop_back();
    EXPECT_TRUE(z.has_value());
    EXPECT_EQ(*z, 222);
    EXPECT_EQ(buf.size(), 0);
}


TEST(TEST_MODULE_NAME, exercise_small_buffer) {
    Deque<int, 8> buf;
    // try the buffer with a small number of items.
    exercise(buf, 6);
    // try the buffer again, but with the head at the middle of the buffer:
    exercise(buf, 6);
}


TEST(TEST_MODULE_NAME, exercise_large_buffer) {
    Deque<int, 8> buf;
    // try the buffer with a large number of items.
    exercise(buf, 19);
    // try again, but with the head in the middle of the buffer:
    exercise(buf, 19);
}


TEST(TEST_MODULE_NAME, exercise_small_copy_and_move) {
    exercise_copy_and_move<int, 8>(7);
}

TEST(TEST_MODULE_NAME, exercise_large_copy_and_move) {
    exercise_copy_and_move<int, 8>(15);
}
