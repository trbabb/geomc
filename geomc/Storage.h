#pragma once
/* 
 * File:   Storage.h
 * Author: tbabb
 *
 * Created on December 24, 2014, 4:06 PM
 */

#include <geomc/geomc_defs.h>

namespace geom {

    
/** @defgroup storage Storage
 *  @brief Classes for arrays with templated static or dynamic size.
 */
    
/**
 * @addtogroup storage
 * @{
 */


enum StoragePolicy {
    /// Dynamically-sized storage will use reference counting to manage memory; all copies of the Storage object share the same underlying array.
    STORAGE_SHARED,
    /// Dynamically-sized storage is assumed to have a single owner, and the underlying array will be duplicated on copy or assignment.
    STORAGE_UNIQUE,
    /// Backing storage is provided and managed by the caller.
    STORAGE_WRAPPED
};


/**
 * @brief Defines a type for storing a length or element count.
 * 
 * Primarily for use in objects which may be either dynamically or statically 
 * sized according to a template parameter. 
 * 
 * For statically-sized objects, no storage is needed for holding the length, 
 * since it is known at compile time. In the static case the storage type is a 
 * zero-length array. Dynamically-sized objects have runtime-determined length, 
 * and so one `index_t` is needed to hold this value. 
 * 
 * @tparam N Count, or 0 if the count is to be determined at runtime.
 * 
 * Example
 * =======
 * 
 * @code#include <geomc/Storage.h>@endcode
 * 
 *     template <index_t N>
 *     class C {
 *         // this member will consume 0 bytes unless N is dynamic:
 *         typename Dimension<N>::storage_t size;
 *         
 *         // ...        
 *         
 *         void setSize(index_t n) { 
 *             // set `size` to `n` if N is dynamic, else do nothing.
 *             Dimension<N>::set(size,n);
 *         }
 * 
 *         index_t getSize(index_t n) {
 *             // return contents of `size` if N is dynamic, else N.
 *             return Dimension<N>::get(size);
 *         }
 *     }
 * 
 */

template <index_t N>
struct Dimension {
    typedef index_t storage_t[0];
    /// Set the count (if stored in `s`) to `v`.
    static inline void set(storage_t s, index_t v) {}
    /// Retrieve the count from `s`.
    static inline index_t get(const storage_t s) { return N; }
};

template <>
struct Dimension<DYNAMIC_DIM> {
    typedef index_t storage_t[1];
    static inline void set(storage_t s, index_t v) { s[0] = v; }
    static inline index_t get(const storage_t s) { return *s; }  
};


//////////////// Storage ////////////////


/**
 * @brief Array storage with templated static or dynamic size.
 * 
 * @tparam T Element type.
 * @tparam N Size of the array, or 0 for dynamic size.
 * 
 * The array will be dynamically allocated if the dimension template argument `N` is 
 * zero; otherwise the storage will be local. Dynamically allocated
 * storage uses reference counting for memory management. 
 * 
 * This is the simplest type of storage, in which the client is responsible for 
 * keeping track of the size of dynamic arrays. (This may save space particularly
 * if multiple/parallel Storage arrays are in use).
 * 
 * `#include <geomc/Storage.h>`
 */
template <typename T, index_t N>
struct Storage {
    /// Data array. Reference-counted if `N` is 0 (dynamic).
    T data[N];
    
    /// Construct a new storage array of length `N` (Static size only).
    Storage() {}
    /// Construct a new storage array of length `n`. `n` is ignored if the array size is not dynamic.
    explicit Storage(index_t n) {}
    /// Construct a new storage array of length `n` and fill with `n` elements from `srcdata`.
    Storage(index_t n, const T* srcdata) { std::copy(srcdata, srcdata + N, data); }
    /// Construct a new storage array and fill with `N` elements from `srcdata` (static size only).
    explicit Storage(const T* srcdata)   { std::copy(srcdata, srcdata + N, data); }
    
    /// Return a pointer to the first element in the storage array.
    inline       T* get()       { return data; }
    /// Return a pointer to the first (const) element in the storage array.
    inline const T* get() const { return data; }
    /// Return the number of elements in the array. (Static size only).
    inline index_t size() const { return N; }
    /// Return a reference to the `i`th element in the array.
    inline T& operator[](index_t i)       { return data[i]; }
    /// Return the `i`th element in the array.
    inline T  operator[](index_t i) const { return data[i]; }
};

template <typename T>
struct Storage<T,DYNAMIC_DIM> {
    std::shared_ptr<T[]> data;
    
    explicit Storage(index_t n):
                    data(new T[n]) {}
    Storage(index_t n, const T* srcdata):
                    data(new T[n]) {
        std::copy(srcdata, srcdata + n, data.get());
    }

    inline       T* get()       { return data.get(); }
    inline const T* get() const { return data.get(); }
    inline T& operator[](index_t i)       { return data[i]; }
    inline T  operator[](index_t i) const { return data[i]; }
};


//////////////// Sized storage ////////////////


/**
 * @brief Array storage with templated static or dynamic size. If the array is
 * dynamic, its length is stored internally and can be queried. 
 * 
 * @tparam T Element type.
 * @tparam N Size of the array, or 0 for dynamic size.
 * 
 * `#include <geomc/Storage.h>`
 */
template <typename T, index_t N>
struct SizedStorage : public Storage<T,N> {
    /// Construct a new SizedStorage of size `N`. (Static size only).
             SizedStorage():Storage<T,N>() {}
    /// Construct a new SizedStorage of size `n`. `n` is ignored if the array size is not dynamic.
    explicit SizedStorage(index_t n):Storage<T,N>(n) {}
    
    /// Return the number of elements in the array.
    index_t size() const { return N; }
};

template <typename T>
struct SizedStorage<T,DYNAMIC_DIM> : public Storage<T,DYNAMIC_DIM> {
    private:
    index_t sz;
    
    public:
    SizedStorage(index_t n):
        Storage<T,DYNAMIC_DIM>(n),
        sz(n) {}
    
    index_t size() const { return sz; };
};

//////////////// Unique storage ////////////////

/**
 * @brief Array storage with templated static or dynamic size, and without 
 * reference counting. 
 *
 * Copy constructions and assignments of will result in a
 * duplication of the underlying array. (If c++11 support is enabled, then
 * move assignments and constructions of dynamically-sized arrays will be 
 * lightweight).
 *
 * Statically-sized instances use local member arrays (as opposed to 
 * heap-allocated ones), so stack-allocated `UniqueStorage`s will also 
 * allocate their underlying arrays on the stack.
 * 
 * @tparam T Element type.
 * @tparam N Size of the array, or 0 for dynamic size.
 * 
 * `#include <geomc/Storage.h>`
 */
template <typename T, index_t N>
struct UniqueStorage {
    /// Data array.
    T data[N];
    
    /// Construct a new UniqueStorage of size `N`. (Static size only).
    UniqueStorage() {}
    /// Construct a new UniqueStorage of size `n`. `n` is ignored if the array size is not dynamic.
    explicit UniqueStorage(index_t n) {}
    /// Construct a new UniqueStorage of size `n`, and copy `n` elements from `srcdata` into the new array.
    UniqueStorage(index_t n, const T* srcdata) {
        std::copy(srcdata, srcdata + n, data);
    }
    /// Construct a new UniqueStorage of size `N`, and copy `N` elements from `srcdata` into
    /// the new array.
    UniqueStorage(std::initializer_list<T> srcdata) {
        std::copy(srcdata.begin(), srcdata.end(), data);
    }
    
    /// Return a pointer to the first element in the storage array.
    inline       T* get()       { return data; }
    /// Return a pointer to the first (const) element in the storage array.
    inline const T* get() const { return data; }
    /// Return the number of elements in the array.
    inline index_t size() const { return N; }
    
    /// Return a reference to the `i`th element in the array.
    inline T& operator[](index_t i)       { return data[i]; }
    /// Return the `i`th element in the array.
    inline T  operator[](index_t i) const { return data[i]; }
};


template <typename T>
struct UniqueStorage<T, DYNAMIC_DIM> {
    T* data;
    index_t sz;
    
    explicit UniqueStorage(index_t n):
            data(new T[n]),
            sz(n) {}

    UniqueStorage(index_t n, const T* srcdata):
            data(new T[n]),
            sz(n) {
        std::copy(srcdata, srcdata + n, data);
    }
    
    UniqueStorage(std::initializer_list<T> srcdata):
            data(new T[srcdata.size()]),
            sz(srcdata.size()) {
        std::copy(srcdata.begin(), srcdata.end(), data);
    }

#if __cplusplus >= 201103L
    UniqueStorage(UniqueStorage<T,0>&& other):
            data(other.data),
            sz(other.sz) {
        other.data = NULL;
        other.sz = 0;
    }
#endif

    UniqueStorage(const UniqueStorage<T,0>& other):
            data(new T[other.sz]),
            sz(other.sz) {
        std::copy(other.data, other.data + sz, data);
    }
    
    ~UniqueStorage() {
        if (data) delete [] data;
    }

#if __cplusplus >= 201103L
    UniqueStorage& operator=(UniqueStorage<T,0>&& other) {
        if (data) delete[] data;
        data = other.data;
        sz = other.sz;
        other.data = NULL;
        other.sz = 0;
        return *this;
    }
#endif

    UniqueStorage& operator=(const UniqueStorage<T,0>& other) {
        if (&other == this) return *this;
        if (sz != other.sz) {
            if (data) delete [] data;
            data = new T[other.sz];
        }
        sz = other.sz;
        std::copy(other.data, other.data + sz, data);
        return *this;
    }
    
    inline       T* get()       { return data; }
    inline const T* get() const { return data; }
    inline index_t size() const { return sz; }
    
    inline T& operator[](index_t idx)       { return data[idx]; }
    inline T  operator[](index_t idx) const { return data[idx]; }
};


//////////////// Wrapped storage ////////////////


// fwd decl
template <typename T, index_t N, StoragePolicy P>
struct GenericStorage;


/**
 * @brief Array storage with templated static or dynamic size, acting
 * as a thin, templated wrapper around a bare array which is memory managed
 * by the caller.
 * 
 * @tparam T Element type.
 * @tparam N Size of the array, or 0 for dynamic size.
 * 
 * This class is mainly used for when ownership semantics are templated.
 * 
 * `#include <geomc/Storage.h>`
 */
template <typename T, index_t N>
struct WrappedStorage {
    T *data;
    typename Dimension<N>::storage_t dim;

    /// Construct a new `WrappedStorage` and use `srcdata` as the backing memory.
    WrappedStorage(index_t sz, T *srcdata):
            data(srcdata) {
        Dimension<N>::set(dim, sz);
    }

    /// Return a pointer to the first element in the storage array.
    inline       T* get()       { return data; }
    /// Return a const pointer to the first element in the storage array.
    inline const T* get() const { return data; }
    /// Return the number of elements in the array.
    inline index_t size() const { return Dimension<N>::get(dim); }
    
    /// Return a reference to the `i`th element in the array.
    inline T& operator[](index_t idx)       { return data[idx]; }
    /// Return the `i`th element in the array.
    inline T  operator[](index_t idx) const { return data[idx]; }
};


//////////////// Generic storage ////////////////


/**
 * @brief Array storage with templated static or dynamic size, and template-selectable
 * ownership policy.
 * 
 * @tparam T Element type.
 * @tparam N Size of the array, or 0 for dynamic size.
 *
 * The ownership policy behavior is as follows:
 * <ul>
 * <li>If the ownership policy is `STORAGE_SHARED`, then dynamically-sized
 * arrays will allocate their own reference-counted memory, and the underlying 
 * storage will be deleted when the last owner is destroyed.</li>
 * <li>If the ownership policy is `STORAGE_UNIQUE`, then dynamically-sized 
 * arrays will allocate their own memory, and the underlying arrays will be
 * duplicated on copy or assignment.</li>
 * <li>If the ownership is `STORAGE_WRAPPED`, then the array will simply wrap a pointer
 * to caller-managed memory, which must be valid for the lifetime of this
 * object, and must contain the minimum number of elements. Duplicates will
 * refer to the same backing memory. No attempt will be made to free the pointer 
 * upon destruction.</li>
 * </ul>
 * 
 * `#include <geomc/Storage.h>`
 */
template <typename T, index_t N, StoragePolicy P>
#ifdef PARSING_DOXYGEN
struct GenericStorage {
#else
struct GenericStorage : public Storage<T,N> {
#endif

    /// Storage class inherited from.
    typedef Storage<T,N> type;

#ifdef PARSING_DOXYGEN
    /// Construct a new array of size `n`, initialized with `srcdata`.
    GenericStorage(index_t n, T* srcdata) {}
#else
    inline GenericStorage(index_t n, const T* srcdata):type(n, srcdata) {}
#endif

    /// Construct a new array of size `n`. Not available for wrapped specializations.
    GenericStorage(index_t n):type(n) {}
    
    GenericStorage(std::initializer_list<T> list):type(list.size()) {
        std::copy(list.begin(), list.end(), type::data);
    }

#ifdef PARSING_DOXYGEN

    /// Return a pointer to the first element in the storage array.
    inline       T* get()       { return data; }
    /// Return a const pointer to the first element in the storage array.
    inline const T* get() const { return data; }
    /// Return the number of elements in the array.
    inline index_t size() const { return Dimension<N>::get(dim); }
    
    /// Return a reference to the `i`th element in the array.
    inline T& operator[](index_t idx)       { return data[idx]; }
    /// Return the `i`th element in the array.
    inline T  operator[](index_t idx) const { return data[idx]; }

#endif

};


// caller-owned specialization
template <typename T, index_t N>
struct GenericStorage<T,N,STORAGE_WRAPPED> : public WrappedStorage<T,N> {

    typedef WrappedStorage<T,N> type;

    GenericStorage(index_t n, T* srcdata):type(n, srcdata) {}
    
    GenericStorage(index_t n):type(0, nullptr) {}

};


// single-owner specialization
template <typename T, index_t N>
struct GenericStorage<T,N,STORAGE_UNIQUE> : public UniqueStorage<T,N> {

    typedef UniqueStorage<T,N> type;

    GenericStorage(index_t n, const T* srcdata):type(n, srcdata) {}

    GenericStorage(index_t n):type(n) {}
    
    GenericStorage(std::initializer_list<T> list):type(list.size()) {
        std::copy(list.begin(), list.end(), type::data);
    }

};


/// @} // group storage


} // end namespace geom
