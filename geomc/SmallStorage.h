#pragma once

#include <geomc/geomc_defs.h>

namespace geom {

/**
 * @brief Array storage which does not allocate from the heap unless the requested 
 * size is larger than a threshold, `N`. 
 * 
 * @tparam T Element type.
 * @tparam N Maximum stack-allocated array size.
 *
 * Useful when a required buffer is expected to be small, but there is no explicit
 * upper bound on the possible size. This allows the stack to be used in the majority
 * case, only deferring to heap allocation when unusually large buffers are needed.
 *
 * On assignment or copy, reallocations will be avoided if the existing buffer is already 
 * the correct size.
 *
 * This buffer will *always* stack-allocate a buffer of size `N` (it will simply not be used
 * when a heap allocation is necessary).
 * 
 * @ingroup storage
 */
template <typename T, size_t N>
struct SmallStorage {
private:
    
    using storage_t = std::aligned_storage_t<sizeof(T), alignof(T)>;
    
    size_t     _capacity = N;
    size_t     _size     = 0;
    storage_t  _buffer[N];
    storage_t* _data     = nullptr;
    
    void _destroy(storage_t* items, size_t n_alive) {
        if constexpr (not std::is_trivially_destructible_v<T>) {
            for (size_t i = 0; i < n_alive; ++i) {
                reinterpret_cast<T*>(items + i)->~T();
            }
        }
    }
    
    void _ensure_capacity(size_t n) {
        if (n > _capacity) {
            size_t new_capacity = std::max(n, _capacity * 2);
            storage_t* new_data = new storage_t[new_capacity];
            T* src = begin();
            for (size_t i = 0; i < n; ++i) {
                new (new_data + i) T{std::move(src[i])};
                src[i].~T(); // destroy the source
            }
            if (_capacity > N) delete [] _data;
            _data     = new_data;
            _capacity = new_capacity;
        }
    }
    
    T* _item(size_t i) {
        return reinterpret_cast<T*>(_data) + i;
    }
    
    const T* _item(size_t i) const {
        return reinterpret_cast<const T*>(_data) + i;
    }
    
public:
    
    static constexpr size_t StaicCapacity = N;

    /// Construct an array with zero items, but capcity for N.
    SmallStorage():_data(_buffer) {}
    
    /// Construct an array containing `n` default-constructed items.
    explicit SmallStorage(size_t n):
        _capacity(std::max(n,N)),
        _size(n),
        _data(n > N ? new storage_t[n] : _buffer)
    {
        for (size_t i = 0; i < n; ++i) {
            // default construct
            new (_data + i) T{};
        }
    }
    
    /// Construct an array by move-constructing the items in the given initializer list.
    SmallStorage(std::initializer_list<T> items):
        _capacity(std::max(items.size(), N)),
        _size(items.size()),
        _data(_capacity > N ? new storage_t[_size] : _buffer)
    {
        for (size_t i = 0; i < _size; ++i) {
            // move construct
            new (_data + i) T{std::move(items.begin()[i])};
        }
    }
    
    /// Construct an array of `n` items by copy-constructing the items from the given array.
    SmallStorage(const T* items, size_t n):
        _capacity(std::max(n, N)),
        _size(n),
        _data(_capacity > N ? new storage_t[_size] : _buffer)
    {
        for (size_t i = 0; i < _size; ++i) {
            new (_data + i) T {items[i]};
        }
    }
    
    /// Construct an array containing `count` copies of `value`.
    SmallStorage(const T& value, size_t count):
        _capacity(std::max(count, N)),
        _size(count),
        _data(_capacity > N ? new storage_t[_size] : _buffer)
    {
        for (size_t i = 0; i < _size; ++i) {
            new (_data + i) T {value};
        }
    }
    
    /// Construct an array containing the elements of `other`.
    SmallStorage(const SmallStorage<T,N>& other):
        _capacity(std::max(other._size, N)),
        _size(other._size),
        _data(_capacity > N ? new storage_t[_size] : _buffer)
    {
        const T* src = other.begin();
        for (size_t i = 0; i < _size; ++i) {
            // copy construct
            new (_data + i) T {src[i]};
        }
    }
    
    /// Move the contents of `other` to this array.
    SmallStorage(SmallStorage<T,N>&& other):
        _capacity(other._capacity),
        _size(other._size),
        _data(_capacity > N ? other._data : _buffer)
    {
        if (_capacity <= N) {
            // move the items from the other container to this one
            T* src = other.begin();
            for (size_t i = 0; i < _size; ++i) {
                // move construct
                new (_data + i) T {std::move(src[i])};
                src[i].~T(); // destroy the source
            }
        } else {
            // array was stolen.
            // empty the other container
            other._capacity = N;
            other._data     = other._buffer;
            other._size     = 0;
        }
    }
    
    /// Destroy this array and its contents.
    ~SmallStorage() {
        _destroy(_data, _size);
        if (_capacity > N) delete[] _data;
        _size = 0;
        _data = nullptr;
    }
    
    /// Move the contents of `other` into this array.
    SmallStorage<T,N>& operator=(SmallStorage<T,N>&& other) {
        _destroy(_data, _size);
        // if we have allocated storage, free it. we're going to clobber
        if (_capacity > N) delete [] _data;
        _capacity = other._capacity;
        _size     = other._size;
        if (other._capacity > N) {
            // the donor container has allocated storage, just grab it
            _data = other._data;
        } else {
            // the donor container has static storage. copy the items.
            _data = _buffer;
            T* src = other.begin();
            for (size_t i = 0; i < _size; ++i) {
                new (_data + i) T {std::move(src[i])};
                src[i].~T();
            }
        }
        other._data     = other._buffer;
        other._capacity = N;
        other._size     = 0;
        return *this;
    }
    
    /// Copy the contents of `other` into this array.
    SmallStorage<T,N>& operator=(const SmallStorage<T,N>& other) {
        if (&other == this) return *this; // no self-assign
        _destroy(_data, _size);
        if (_capacity < other._size) {
            // we cannot hold the incoming items without allocating.
            // create the array anew
            if (_capacity > N) delete [] _data;
            _capacity = std::min(other._capacity, 2 * other._size);
            _data     = new storage_t[_capacity];
        } else {
            // todo: we could copy-assign to the live items,
            //   but we'd probably also want to statically check that
            //   copy-assign exists to avoid splitting the loop into live/non-live
            //   when it isn't necessary
        }
        const T* src = other.begin();
        for (size_t i = 0; i < other._size; ++i) {
            new (_data + i) T {src[i]};
        }
        _size = other._size;
        return *this;
    }
    
    /// Access the `i`th item.
    T& operator[](size_t i) {
        return reinterpret_cast<T*>(_data)[i];
    }
    
    /// Access the `i`th const item.
    const T& operator[](size_t i) const {
        return reinterpret_cast<const T*>(_data)[i];
    }
    
          T* begin()       { return reinterpret_cast<      T*>(_data); }
    const T* begin() const { return reinterpret_cast<const T*>(_data); }
    
          T* end()         { return reinterpret_cast<      T*>(_data) + _size; }
    const T* end()   const { return reinterpret_cast<const T*>(_data) + _size; }
    
    /// Add a copy of the given item to the end of the array, increasing its size by 1.
    T& push_back(const T& item) {
        _ensure_capacity(_size + 1);
        T* place = reinterpret_cast<T*>(_data + _size);
        new (place) T {item};
        _size += 1;
        return *place;
    }
    
    /// Construct a new item at the end of the array from the given constructor arguments.
    template <typename... Args>
    T& emplace_back(Args&&... args) {
        _ensure_capacity(_size + 1);
        T* place = reinterpret_cast<T*>(_data + _size);
        new (place) T {std::forward<Args>(args)...};
        _size += 1;
        return *place;
    }
    
    /// Remove the last item from the array and return it, if one exists.
    std::optional<T> pop_back() {
        if (_size == 0) return std::nullopt;
        T* last = _item(_size - 1);
        T out = std::move(*last);
        last->~T();
        _size -= 1;
        return out;
    }
    
    /// Access the last item in the array.
    T& back() {
        return begin()[_size - 1];
    }
    
    /// Access the last const item in the array.
    const T& back() const {
        return begin()[_size - 1];
    }
    
    /// Access the first item in the array.
    T& front() {
        return *begin();
    }
    
    /// Access the first const item in the array.
    const T& front() const {
        return *begin();
    }
    
    /// Return a pointer to the underlying array that backs this `SmallStorage`.
          T* data()       { return reinterpret_cast<      T*>(_data); }
    /// Return a const pointer to the underlying array that backs this `SmallStorage`.
    const T* data() const { return reinterpret_cast<const T*>(_data); }
    
    /// Return the number of items in the array.
    size_t size() const { return _size; }
    
    /// Return the number of items the array can hold without allocating from the heap.
    static constexpr size_t static_capacity() {
        return N;
    }
    
    /// Change the number of items in the array; constructing new items from the arguments
    /// if any are provided, and destroying any excess items.
    template <typename... Args>
    void resize(size_t count, const Args&... args) {
        if (count > _size) {
            // construct the missing items
            _ensure_capacity(count);
            for (size_t i = _size; i < count; ++i) {
                new (_data + i) T {args...};
            }
        } else if constexpr (not std::is_trivially_destructible_v<T>) {
            // destroy the excess items
            for (size_t i = count; i < _size; ++i) {
                _item(i)->~T();
            }
        }
        _size = count;
    }
    
    /// Ensure there is capacity for at least `count` items without reallocating.
    void reserve(size_t count) {
        _ensure_capacity(count);
    }
    
    /// Return true iff the array contains no items.
    bool empty() const {
        return _size == 0;
    }
    
    /// Return the number of items the array can contain without reallocating.
    size_t capacity() const {
        return _capacity;
    }
    
    /// Remove all the items from the array.
    void clear() {
        _destroy(_data, _size);
        _size = 0;
    }
    
    /// Return true iff all the items in `this` and `other` are equal and in equal order.
    bool operator==(const SmallStorage<T,N>& other) const {
        if (_data == other._data) return true;
        if (_size != other._size) return false;
        for (size_t i = 0; i < _size; ++i) {
            if ((*this)[i] != other[i]) return false;
        }
        return true;
    }
    
}; // struct SmallStorage

} // namespace geom
