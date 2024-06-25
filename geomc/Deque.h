#include <type_traits>
#include <utility>
#include <geomc/Storage.h>
#include <geomc/function/Utils.h>

// todo: this guy wants an iterator
// todo: head + tail might be more elegant if tail points to last and not just-beyond-last.


namespace geom {

/**
 * @addtogroup storage
 * @{
 */

/**
 * @brief A lightweight double ended queue backed by a flat array.
 *
 * A deque can accommodate adding elements to either the front or the back
 * of the list in constant time. Random access is also constant time.
 *
 * Furthermore, indexing off the end of the deque wraps around to the
 * beginning again, as in a circular buffer. This works in both directions.
 *
 * @tparam T Element type.
 * @tparam N Static capacity of the buffer. Adding more than this number of elements
 * to the buffer will incur a heap memory allocation. To always use heap allocation, 
 * pass zero to this parameter, in which case memory will be allocated when the first
 * element is added.
 */
template <typename T, index_t N>
class Deque {
    
    typedef typename std::aligned_storage<sizeof(T), alignof(T)>::type storage_t;
    
    // `head + size` scheme chosen over `head + tail` to avoid 
    // ambiguity of "full vs. empty buffer"
    
    // We usage storage_t both to admit objects with no default constructor,
    // and also to avoid unnecessary constructions/destructions of items 
    // in the empty area of the buffer.
    
    // We manage our own static/dynamic buffer instead of delegating to SmallStorage
    // because we need to manage the liveness of the elements in the buffer.
    // A SmallStorage<storage_t> won't properly call move, copy, or assignment
    // operators because storage_t hides the overlying class, and of course 
    // SmallStorage does not know which items are live.
    
    storage_t  _buf[N];
    storage_t* _data;
    index_t    _capacity;
    index_t    _head;
    index_t    _size;
    
public:
    
    /// Construct a new empty circular buffer, with capacity for `N` items.
    Deque():
        _data(_buf),
        _capacity(N),
        _head(0),
        _size(0) {}
    
    /// Construct a new empty circular buffer, with space for at least `capacity` items.
    Deque(index_t capacity):
        _data((capacity > N) ? (new storage_t[capacity]) : _buf),
        _capacity(std::max(capacity, N)),
        _head(0),
        _size(0) {}
    
    /**
     * @brief Construct a new circular buffer by copying `count` items in the sequence
     * starting at `begin`.
     */
    template <typename InputIterator>
    Deque(InputIterator begin, index_t count):
            _data((count > N) ? (new storage_t[count]) : _buf),
            _capacity(std::max(count, N)),
            _head(0),
            _size(count)
    {
        // placement construct all the items in the provided range:
        for (index_t i = 0; i < count; ++i, ++begin) {
            T* p = get() + i;
            new (p) T(*begin);
        }
    }
    
    /// Construct a new circular buffer containing copies of all the items in `other`.
    Deque(const Deque<T,N>& other):
            _data((other._size > N) ? (new storage_t[other._capacity]) : _buf),
            _capacity(std::max(other._capacity, N)),
            _head(0),
            _size(other._size)
    {
        // copy-construct the populated items:
        for (index_t i = 0; i < _size; ++i) {
            T* p = get() + i;
            new (p) T(other[i]);
        }
    }
    
    /// Move the contents of `other` to a new Deque
    Deque(Deque<T,N>&& other):
            _data((other._capacity > N) ? other._data : _buf),
            _capacity(other._capacity),
            _head((other._capacity > N) ? other._head : 0),
            _size(other._size)
    {
        if (other._capacity <= N) {
            // move-construct individual items.
            for (index_t i = 0; i < _size; ++i) {
                T* p = get() + i;
                new (p) T(std::move(other[i]));
                other[i].~T();
            }
        } else {
            // we took ownership of `other`'s allocated array.
            // turn `other` back into a static buffer.
            other._capacity = N;
            other._data = other._buf;
        }
        other._size = 0;
        other._head = 0;
    }
    
    // Destroy this buffer and all the items within.
    virtual ~Deque() {
        for (index_t i = 0; i < size(); ++i) {
            item(i)->~T();
        }
        if (_capacity > N and _capacity > 0) {
            delete [] _data;
            _capacity = 0;
            _data = nullptr;
        }
        _size = 0;
    }
    
    /// Assignment operator
    Deque& operator=(const Deque& other) {
        index_t old_size = size();
        index_t new_size = other.size();
        index_t n = std::min(old_size, new_size);
        reserve(other._size);
        // delegate assignment to T
        for (index_t i = 0; i < n; ++i) {
            this->operator[](i) = other[i];
        }
        // copy-construct new items
        for (index_t i = n; i < new_size; ++i) {
            T* p = item(i);
            new (p) T(other[i]);
        }
        // destroy extra items
        for (index_t i = n; i < old_size; ++i) {
            T* p = item(i);
            p->~T();
        }
        _size = other._size;
        
        return *this;
    }
    
    /// Move the contents of `other` to this buffer
    Deque& operator=(Deque&& other) {
        // destroy all our stuff
        for (index_t i = 0; i < _size; ++i) {
            item(i)->~T();
        }
        if (other._capacity > N) {
            // `other` has a dynamic array, just take ownership of it
            if (_capacity > N) {
                // delete our own dynamic array
                delete [] _data;
            }
            _data     = other._data;
            _capacity = other._capacity;
            _head     = other._head;
            // turn `other` back into a static buffer
            other._data     = other._buf;
            other._capacity = N;
        } else {
            // `other`'s array cannot be moved; copy items instead.
            // we do not need to reserve(), obviously, because _size < N.
            // if we have space for more than N, that's fine too.
            _head = 0;
            for (index_t i = 0; i < other._size; ++i) {
                T* p = get() + i;
                // move-construct a new item in our own buffer
                new (p) T(std::move(other.operator[](i)));
                // destroy its source
                other.operator[](i).~T();
            }
        }
        
        _size = other._size;
        other._size = 0;
        other._head = 0;
        
        return *this;
    }
    
    /**
     * @brief Get the `i`th element in the buffer. 
     * 
     * Indicies beyond the end of the buffer will wrap around again to the beginning.
     * Negative indices are permitted and count from the end of the buffer, with -1
     * denoting the last element in the buffer.
     */
    T& operator[](index_t i) {
        i = positive_mod(i, _size);
        return *item(i);
    }
    
    /**
     * @brief Get the `i`th (const) element in the buffer. 
     * 
     * Indicies beyond the end of the buffer will wrap around again to the beginning.
     * Negative indices are permitted and count from the end of the buffer, with -1
     * denoting the last element in the buffer.
     */
    const T& operator[](index_t i) const {
        i = positive_mod(i, _size);
        return *item(i);
    }
    
    /**
     * @brief Equality operator.
     *
     * Two Deques are equal iff they contain equal elements in equal order.
     */
    bool operator==(const Deque<T, N>& other) {
        if (this == &other)       return true;  // literally the same object
        if (_size != other._size) return false; // different number of items; unequal
        for (index_t i = 0; i < _size; ++i) {
            if ((*this)[i] != other[i]) return false;
        }
        return true;
    }
    
    /// Inequality operator.
    bool operator!=(const Deque<T, N>& other) {
        return not ((*this) == other);
    }
    
    /// Return the number of items in the buffer.
    index_t size() const {
        return _size;
    }
    
    /**
     * @brief Return the total number of items that can be accommodated without an additional
     * memory allocation.
     */
    index_t capacity() const {
        return _capacity;
    }
    
    /// Add an element to the end of the buffer.
    template <typename U> // templated so that we can deduce and collapse lvalues
    T& push_back(U&& t) {
        reserve(_size + 1);
        T* p = item(_size);
        // placement new, using move semantics if applicable:
        new (p) T(std::forward<U>(t));
        _size += 1;
        return *p;
    }
    
    /**
     * @brief Add an element to the beginning of the buffer.
     * 
     * This increases the indices of all the existing elements by one.
     * The new element will have index zero.
     */
    template <typename U> // templated so that we can deduce and collapse lvalues
    T& push_front(U&& t) {
        reserve(_size + 1);
        const index_t n = capacity();
        _head = (_head + n - 1) % n;
        T* place = get() + _head;
        // placement new, using move semantics if applicable:
        new (place) T(std::forward<U>(t));
        _size += 1;
        return *place;
    }
    
    /**
     * @brief Remove the first element in the buffer, if there is one.
     * 
     * This decreases the indices of all the remaining elements by one.
     */
    std::optional<T> pop_front() {
        if (_size == 0) return std::nullopt;
        const index_t old = _head;
        _head  = (_head + 1) % capacity();
        _size -= 1;
        // move-construct a new T into the return value
        T ret(std::move(get()[old]));
        // destroy the old source of the return value
        get()[old].~T();
        return ret;
    }
    
    /// Remove the element at the end of the buffer, if there is one.
    std::optional<T> pop_back() {
        if (_size == 0) return std::nullopt;
        index_t i = (_head + _size - 1) % capacity();
        _size -= 1;
        // move-construct a new T into the return value
        T ret(std::move(get()[i]));
        // destory the source of the return value
        get()[i].~T();
        return ret;
    }
    
    /// Return a reference to the item at the beginning of the buffer.
    T& front() {
        return get()[_head];
    }
    
    /// Return a const reference to the item at the beginning of the buffer.
    const T& front() const {
        return get()[_head];
    }
    
    /// Return a reference to the item at the end of the buffer.
    T& back() {
        return *item(_size - 1);
    }
    
    /// Return a const reference to the item at the end of the buffer.
    const T& back() const {
        return *item(_size - 1);
    }
    
    /// Empty the buffer of all items.
    void clear() {
        // destroy all items.
        for (index_t i = 0; i < _size; ++i) {
            item(i)->~T();
        }
        _head = 0;
        _size = 0;
    }
    
    /** 
     * @brief Increase the capacity of the circular buffer to a value that's
     * greater or equal to `new_cap`. If `new_cap` is greater than the current
     * capacity, then new storage is allocated; otherwise the method does nothing.
     */
    void reserve(index_t new_cap) {
        if (new_cap > _capacity) {
            index_t resize_to = std::max(_capacity * 2, new_cap);
            storage_t* array  = new storage_t[resize_to];
            
            // move the items into the new array
            for (index_t i = 0; i < _size; ++i) {
                T* p = reinterpret_cast<T*>(array + i);
                // move-construct
                new (p) T(std::move(*item(i)));
                // destroy the source
                item(i)->~T();
            }
            
            if (_capacity > N) {
                delete [] _data;
            }
            
            _data = array;
            _head = 0;
            _capacity = resize_to;
        }
    }
        

    
protected:
    
    // get the buffer
    inline T* get() {
        return reinterpret_cast<T*>(_data);
    }
    
    // get the const buffer
    inline const T* get() const {
        return reinterpret_cast<const T*>(_data);
    }
    
    // get the place for the `i`th item
    inline T* item(index_t i) {
        return get() + (_head + i) % capacity();
    }
    
    // get the place for the `i`th const item
    inline const T* item(index_t i) const {
        return get() + (_head + i) % capacity();
    }
    
};

/// @} // end group "storage"
    
} // end namespace geom
