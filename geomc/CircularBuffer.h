#include <geomc/Storage.h>

namespace geom {

/**
 * @addtogroup storage
 * @{
 */

/**
 * @brief A lightweight circular buffer.
 *
 * A circular buffer can accommodate adding elements to either the front or the back
 * of the list in constant time.
 *
 * @tparam T Element type.
 * @tparam N Static capacity of the buffer. Adding more than this number of elements
 * to the buffer will incur a heap memory allocation. To always use heap allocation, 
 * pass zero to this parameter.
 */
template <typename T, index_t N>
class CircularBuffer {
    
    // `head + size` scheme chosen over `head + tail` to avoid 
    // ambiguity of "full vs. empty buffer"
    
    SmallStorage<T,N> _data;
    index_t           _head;
    index_t           _size;
    
public:
    
    /// Construct a new empty circular buffer.
    CircularBuffer():
        _data(N),
        _head(0),
        _size(0) {}
    
    /// Construct a new empty circular buffer, with space for at least `capacity` items.
    CircularBuffer(index_t capacity):
        _data(capacity),
        _head(0),
        _size(0) {}
    
    /**
     * @brief Get the `i`th element in the buffer. 
     * 
     * Indicies beyond the end of the buffer will wrap around again to the beginning.
     * Negative indices are permitted and count from the end of the buffer. -1 denotes the 
     * last element in the buffer.
     */
    inline const T& operator[](index_t i) const {
        i = positive_mod(i, _size);
        return _data.get()[(_head + i) % _data.size()];
    }
    
    /**
     * @brief Get the `i`th element in the buffer. 
     * 
     * Indicies beyond the end of the buffer will wrap around again to the beginning.
     * Negative indices are permitted and count from the end of the buffer. -1 denotes the 
     * last element in the buffer.
     */
    inline T& operator[](index_t i) {
        i = positive_mod(i, _size);
        return _data.get()[(_head + i) % _data.size()];
    }
    
    /// Return the number of items in the buffer.
    inline index_t size() const {
        return _size;
    }
    
    /// Return the total number of items that can be accommodated without an additional memory allocation.
    inline index_t capacity() const {
        return _data.size();
    }
    
    /// Add an element to the end of the buffer.
    inline void push_back(const T& t) {
        check_size();
        _data.get()[(_head + _size) % _data.size()] = t;
        _size += 1;
    }
    
    /**
     * @brief Add an element to the front of the buffer.
     * 
     * This increases the indices of all the existing elements by one.
     * The new element will have index zero.
     */
    inline void push_front(const T& t) {
        check_size();
        const index_t n = _data.size();
        _head = _data.get() + (_head + n - 1) % n; 
        _data.get()[_head] = t;
        _size += 1;
    }
    
    /**
     * @brief Remove the first element in the buffer.
     * 
     * This decreases the indices of all the remaining elements by one.
     */
    T pop_front() {
        const index_t old = _head;
        _head  = (_head + 1) % _data.size();
        _size -= 1;
        return _data.get()[old];
    }
    
    /// Remove the element at the end of the buffer.
    T pop_back() {
        index_t i = (_head + _size) % _data.size();
        size -= 1;
        return _data.get()[i];
    }
    
    /// Empty the buffer of all items.
    inline void clear() {
        _head = 0;
        _size = 0;
    }
    
    
protected:
    
    // check if a memory reallocation is necessary, and do so if needed.
    inline void check_size() {
        const index_t n = _data.size();
        if (_size >= n) {
            _data.resize(n * 2);
            // make the "wraparound" part contiguous in the new buffer
            // todo: SmallStorage does its own copy. wbn to avoid that.
            const index_t extra = (_head + _size) % n;
            std::copy(
                std::make_move_iterator(_data.get()),
                std::make_move_iterator(_data.get() + extra),
                _data.get() + n);
        }
    }
    
};

/// @} // end group "storage"
    
} // end namespace geom
