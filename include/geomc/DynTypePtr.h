/*
 * AnyPtr.h
 * 
 * This class stores a pointer to an object of unknown type.
 * It's safer than void*, because it stores the type_info of the
 * object that was assigned to it. This is like a boost::any, except
 * that the user provides the storage-- thus there are no malloc()s/news,
 * so this class has the potential to be much faster.
 * 
 * 
 * CAREFUL: polymorphism.
 *   - A is the base of B
 *   - We store a B*
 *   - We want to use our stored B* as an A*
 *   - Conversion fails because type(B) != type(A)
 *   - Will dynamic_cast work on void*? sounds sketchy.
 *     - A: No.
 *   - Need access to the vtable, but c++ is bad at being reflective.
 *   - This is pretty much impossible without manual input of some kind.
 *
 *  Created on: Aug 7, 2011
 *      Author: tbabb
 */

#ifndef ANYPTR_H_
#define ANYPTR_H_

#include <typeinfo>

using std::type_info;

class bad_dyntype_ptr_cast : public std::bad_cast {
public:
    virtual const char * what() const throw() {
        return "geom::bad_dyntype_ptr_cast: "
               "failed conversion from geom::DynTypePtr";
    }
};

/****************************
 * AnyPtr class             *
 ****************************/

class DynTypePtr {
public:
    
    ///// Structors /////
    
    DynTypePtr():
            ptr(0),
            tinfo(&typeid(void)){
        //do nothing else
    }
    
    template <typename T> DynTypePtr(T *var):
            ptr(var),
            tinfo(&typeid(T)){
        //do nothing else
    }
    
    virtual ~DynTypePtr(){}
    
    ///// Access /////
    
    bool isNull(){
        return ptr == 0;
    }
    
    const type_info& getType(){
        return *tinfo;
    }

    //return a reference to the pointed object
    template <typename T> T& item() {
        if (typeid(T) == *tinfo){
            return *(T*)ptr;
        } else {
            throw bad_dyntype_ptr_cast();
        }
    }
    
    //return a const reference to the pointed object
    template <typename T> const T& item() const {
        if (typeid(T) == *tinfo){
            return *(T*)ptr;
        } else {
            throw bad_dyntype_ptr_cast();
        }
    }
    
    //return the raw pointer, as the correct type.
    template <typename T> T* get() const {
        if (typeid(T) == *tinfo){
            return (T*)ptr;
        } else {
            throw bad_dyntype_ptr_cast();
        }
    }
    
    /*
     * You would really like to implement operator*(), wouldn't you?
     * But you can't, because it would be templated, and template inference
     * doesn't happen when the return type is the only way to infer the
     * template parameter. So you would have to write obj.operator*<Type>()
     * literally everwhere, which is just ugly.
     */
    
    template <typename T> DynTypePtr& operator=(T *var){
        ptr = var;
        tinfo = &typeid(T);
        return *this;
    }
    
protected:
    void *ptr;
    const type_info *tinfo;
    
};

#endif /* ANYPTR_H_ */
