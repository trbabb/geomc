Code Style
==========

K&R indent style. 

* Single space between braces and parens. 
* Single space after `if` and `else`.
* No spaces on the inside of parens.
* Single space after commas.

Example:

    if (var < 1) {
        do_code(var, var2);
    } else {
        do_something_else(var, var3);
    }

Only single-line `if`s can be without braces:

    if (test > 0) do_code(x);

If there is a long block of assignments, line up the equal signs to make it look like a more readable table:

    int      my_var = 1;
    int    some_var = 2;
    int another_var = 6;

Spaces around all basic arithmetic operators:

    x * x + y * y + sin(z);

None of this:

    x*x+y*y+sin(z); // No.

Double line breaks between functions, single line breaks to break up logical blocks within a function:

    float f(int a) {
        int z;
        if (a < 1) {
           z = do_something_with(a);
        }
        
        int b = z * a + 1;
        int c = b - g(a);
        return c * std::sqrt(b);
    }
    
    
    void g(int z) {
        int d = do_something_else_with(z);
        return h(d, z * z + 1);
    }

Classes list member variables first, followed by constructors and destructors, followed by member functions:

    class Foo {
        float x;
        float y;
        Quat<float> q;
        
        Foo() {}
        
        Vec<float,3> do() {
            return q * Vec<float,3>(x,y,1);
        }
    }


Programming
===========

Do not allocate memory if it can be avoided. Stack variables and (reasonably-size) stack-allocated arrays are always 
preferable where possible.

Precision
---------

Do not use logic with "epislons" if it can be avoided, especially since most code is templated, and the scale of a sensible epsilon would be dependent on the template parameter type. Arithmetic accuracy should fail naturally at machine precision. 

Yes:

    T x = vec.mag();
    if (x != 0) vec /= x;

No:

    // calculation becomes invalid around the scale of EPSILON, rather than 1ulp.
    T x = vec.mag();
    if (x > EPSILON) vec /= x;
    
(There is one (rare) exception, which is when choosing between two strategies that each return a more accurate result depending on scale; for example, a complete formula vs. its taylor series approximation near an unstable point. In these cases the numerical considerations should be well understood, and the two strategies should perform demonstrably better in their respective domains).

Loops
------

Loop variables should be type `index_t`, which will be 64-bit on 64-bit platforms, allowing for traversal of very large arrays:

    for (index_t i = 0; i < n; i++) {
        v[i] += k[i];
    }

Argument passing style
----------------------

Objects that are modified by a function should be passed as pointers:
    
    void foo(Object *out, int in) {
        out->x = some_function(n);
    }

Objects that are passed by reference should be *always* declared const:

    int bar(const Object &obj, int k) {
        return some_function(obj.x, k);
    }

The above makes it obvious at the call site whether a function will modify an argument:

    invert(&mtx, src);

