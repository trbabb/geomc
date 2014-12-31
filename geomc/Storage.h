/* 
 * File:   Storage.h
 * Author: tbabb
 *
 * Created on December 24, 2014, 4:06 PM
 */

#ifndef STORAGE_H
#define	STORAGE_H

#include <boost/smart_ptr/shared_array.hpp>
#include <geomc/geomc_defs.h>

namespace geom {
    
/** @defgroup storage Storage
 *  @brief Classes for arrays with templated static or dynamic size.
 */
    
/**
 * @addtogroup storage
 * @{
 */

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

/**
 * @brief Array storage with templated static or dynamic size.
 * 
 * @tparam T Element type.
 * @tparam N Size of the array, or 0 for dynamic size.
 * 
 * The array will be dynamically allocated if the dimension template argument `N` is 
 * zero; otherwise the storage will be allocated on the stack. Dynamically allocated
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
    boost::shared_array<T> data;
    
    explicit Storage(index_t n):
                    data(new T[n]) {}
    
    inline       T* get()       { return data.get(); }
    inline const T* get() const { return data.get(); }
    inline T& operator[](index_t i)       { return data[i]; }
    inline T  operator[](index_t i) const { return data[i]; }
};


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

/**
 * @brief ArrayStorage with templated static or dynamic size, and without 
 * reference counting (allocated data is deleted upon destruction). 
 * 
 * @tparam T Element type.
 * @tparam N Size of the array, or 0 for dynamic size.
 * 
 * Preferred for when ownership never changes hands and/or memory management is 
 * not needed.
 * 
 * `#include <geomc/Storage.h>`
 */
template <typename T, index_t N>
struct UnmanagedStorage {
    /// Data array.
    T data[N];
    
    /// Construct a new UnmanagedStorage of size `N`. (Static size only).
    UnmanagedStorage() {}
    /// Construct a new UnmanagedStorage of size `n`. `n` is ignored if the array size is not dynamic.
    explicit UnmanagedStorage(index_t n) {}
    
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
struct UnmanagedStorage<T,DYNAMIC_DIM> {
    T *data;
    index_t sz;
    
    explicit UnmanagedStorage(index_t n):
        data(new T[n]),
        sz(n) {}
    
    ~UnmanagedStorage() {
        delete [] data;
    }
    
    inline       T* get()       { return data; }
    inline const T* get() const { return data; }
    inline index_t size() const { return sz; }
    
    inline T& operator[](index_t idx)       { return data[idx]; }
    inline T  operator[](index_t idx) const { return data[idx]; }
}; 

/// @} // group linalg

} // end namespace geom

#endif	/* STORAGE_H */

