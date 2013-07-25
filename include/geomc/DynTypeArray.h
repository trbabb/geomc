/*
 * DynTypeArray.h
 * 
 * Safely hold a pointer to a shared array of unknown type.
 * 
 * Pointer and indexing operator overloads are not supplied. This is unfortunately
 * not possible (or at least not worth it) because they would require template 
 * parameters at every use. Instead use get() for the raw array, and 
 * item() for elements. 
 * 
 * DynTypeArray can be copied by value; the underlying array is shared
 * by the copied instances. When the last instance is destroyed, the array is 
 * freed. Do not explicitly free the array, and do not hold the array
 * past the lifetime of its owner DynTypeArray.
 * 
 * Note that iterating on the result of get() will be faster than calling item()
 * at each iteration, since with get() the type check will happen only
 * once.
 * 
 *  Created on: Aug 8, 2011
 *      Author: tbabb
 */

#ifndef ANYARRAY_H_
#define ANYARRAY_H_

#include <typeinfo>
#include <boost/shared_array.hpp>
#include <boost/intrusive_ptr.hpp>

/*
 * Possible optimizations:
 *   - Right now, we'll see a double allocate: Typed array container + array itself.
 *   - Consider allocating from a memory pool and using placement new.
 *   - Currently storing redundant ref to component typeid.
 *     - Could save space by querying the Holder
 *     - Might save time by removing an indirection, however.
 */

// note: it is not clear to me (2/12/12) why we do not use a single-layer wrapper and a void*.
//       I am guessing it might be some complexity due to type erasure? I bet it's that 
//       delete[] can't do the right thing on void*.

using boost::shared_array;
using boost::intrusive_ptr;
using std::type_info;

class DynTypeArray;

class bad_anyarray_cast : public std::bad_cast {
public:
    virtual const char * what() const throw() {
        return "geom::bad_dyntype_array_cast: "
               "failed conversion from geom::DynTypeArray";
    }
};

/***************************
 * AnyArray class          *
 ***************************/

class DynTypeArray {
    // this is the class that holds the "outer array wrapper", and acts as the public interface.
    // it cannot actually *be* the outer wrapper, unfortunately, because the outer wrapper is
    // (and must be) pure virtual (because the get/set methods are owned by the inner wrapper). 
    // Thus we hold one instead.
protected:

    /***************************
     * Outer wrapper class     *
     ***************************/

    class ArrayHolderWrapper {
        // this outer wrapper exists so that we can store arrays
        // without knowing the type of the array.
        // we allocate first for the wrapper, then for the array.
        // an optimization could allocate both at the same time.
    public:
        
        ArrayHolderWrapper(size_t elem_size):
            _refct(0),
            _elem_size(elem_size){}
        virtual ~ArrayHolderWrapper(){}
        
        virtual const type_info& getElementType() const = 0;
        virtual       void *elem_ptr(size_t idx) = 0;
        virtual const void *elem_ptr(size_t idx) const = 0;
        
        inline size_t elem_size() const { return _elem_size; }

    private:
        friend void intrusive_ptr_add_ref(ArrayHolderWrapper *p);
        friend void intrusive_ptr_release(ArrayHolderWrapper *p);
        
    protected:
        size_t _refct;
        size_t _elem_size;
    };

    /***************************
     * Inner wrapper class     *
     ***************************/

    template <typename T> class ArrayHolder : public ArrayHolderWrapper {
        //this inner wrapper is what instantiates and carries the shared_array
    protected:
        T *_array;
    public:
        ArrayHolder(size_t n):
            ArrayHolderWrapper(sizeof(T)),
            _array(new T[n]){}
        
        virtual ~ArrayHolder(){delete [] _array;}
        
        const type_info& getElementType() const{
            return typeid(T);
        }
        
        void *elem_ptr(size_t idx){
            return &(_array[idx]);
        }
        
        const void *elem_ptr(size_t idx) const {
            return &(_array[idx]);
        }
        
    };
    
    /***************************
     * DynTypeArray guts       *
     ***************************/
    
    intrusive_ptr<ArrayHolderWrapper> arrayHolder;
    size_t len;
    const type_info *tinfo; 
    
public:
    
    DynTypeArray():
            arrayHolder(0),
            len(0){}
    
    virtual ~DynTypeArray(){}
    
    /* This is necessary (instead of templating the constructor) because constructors of non-template 
     * classes cannot be templated in c++. It is legal to declare them, but impossible to call them. Way to go, Bjarne. */
    template <typename T> void initialize(size_t n){
        len = n;
        arrayHolder = intrusive_ptr<ArrayHolderWrapper>(new ArrayHolder<T>(n));
    }
    
    template <typename T> T& item(size_t idx){
        if (matchesType<T>()){
            return *((T*)arrayHolder->elem_ptr(idx));
        } else {
            throw bad_anyarray_cast();
        }
    }
    
    template <typename T> const T& item(size_t idx) const {
        if (matchesType<T>()){
            return *((T*)arrayHolder->elem_ptr(idx));
        } else {
            throw bad_anyarray_cast();
        }
    }
    
    template <typename T> T* get(){
        if (!arrayHolder){
            return NULL;
        } else if (matchesType<T>()){
            return (T*)arrayHolder->elem_ptr(0);
        } else {
            throw bad_anyarray_cast();
        }
    }
    
    template <typename T> inline bool matchesType(){
        if (arrayHolder){
            return typeid(T) == arrayHolder->getElementType();
        } else {
            return false;
        }
    }
    
    const type_info& getElementType(){
        if (arrayHolder){
            return arrayHolder->getElementType();
        } else {
            return typeid(void);
        }
    }
    
    size_t size(){
        return len;
    }
    
    size_t bytes_held(){
        return len * arrayHolder->elem_size();
    }

private:
    friend void intrusive_ptr_add_ref(ArrayHolderWrapper *p);
    friend void intrusive_ptr_release(ArrayHolderWrapper *p);
};

/***************************
 * Reference counting      *
 ***************************/

inline void intrusive_ptr_add_ref(DynTypeArray::ArrayHolderWrapper *p){
    p->_refct++;
}

inline void intrusive_ptr_release(DynTypeArray::ArrayHolderWrapper *p){
    if ((p->_refct -= 1) <= 0){
        delete p;
    }
}

/***************************
 * Helper functions        *
 ***************************/

template <typename T> DynTypeArray makeDynTypeArray(size_t sz){
    DynTypeArray ary;
    ary.initialize<T>(sz);
    return ary;
}

#endif /* ANYARRAY_H_ */
