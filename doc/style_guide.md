Code Style
==========

Braces and indents
------------------

**Never** indent with tabs; use four spaces.

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
    
Alignment
-----------

It is good to highlight similar structure by aligning code into a visual table:

    int myvar_1   =  func(1,   2, 3);
    int myvar_123 =  func(1,  55, 0);
    int myvar_2   =  func(2, 255, 0);
    int myvar_3   = thing(1,   0, 0);

This makes it easy to examine similarities and differences, reducing mental workload. When in doubt about whether to justify left or right, prefer to keep similar elements aligned vertically.

This principle is not restricted to blocks of assignments; it works elsewhere:

    class Foo {
        int       thing_foo(int a, int b);
        int       other(int a);
        SomeClass boop(int a);
    };

Do apply common sense; tabular alignment need not be preferred where it creates excessively long lines or visual clutter.

Arithmetic
----------

Spaces around all basic arithmetic operators:

    x * x + y * y + sin(z);

None of this:

    x*x+y*y+sin(z); // No.

Line breaks
-----------

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

Classes
-------

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

Declarations
------------

Declare each variable in its own statement:

    int* a;
    int* b;

Never use comma separated declarations with decorated types. C++'s design contains a mistake where the `*` associates with the variable name instead of the type, inviting misleading constructions like:

    // NEVER:
    int* a, b;

In the above, `a` is an `int*`, but `b` is an `int`. To avoid this confusion, separate declarations must always be used instead.

Comparisons
-----------

Modern languages support `and` and `or` keywords. These are easier to read, have been part of the C/C++ standard since the early days. They should be preferred over `||` and `&&`:

    if ((cond1 and cond2) or cond3) {
        do_a_thing(x);
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
    
There is one (rare) exception, which is when choosing between two strategies that each return a more accurate result depending on scale; for example, a complete formula vs. its taylor series approximation near an unstable point. In these cases the numerical considerations should be well understood, and the two strategies should perform demonstrably better in their respective domains:

    template <typename T>
    T complicated_func(T input) {
        if (input < 1 and input > 0) {
            // unstable region better approximated with a different method:
            taylor_approx(input);
        } else {
            direct_eval(input);
        }
    }

Counters
------

Counter variables should be type `index_t`, which will be 64-bit on 64-bit platforms:

    for (index_t i = 0; i < n; i++) {
        v[i] += k[i];
    }

Using `index_t` guarantees that large counts of objects can be handled properly on 64-bit systems without overflow.

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

Optimization
------------

Structure and factor code to be readable and maintainable. These two factors *always* take precedence over optimizations that modern compilers are known to perform. In addition, `geomc` is intended to be compatible with multiple compilers, so a low-level factorization intended to help one might hurt another. 

Examples of transformations that should **not** be written in source unless they add clarity or reduce code duplication:

* Constant folding
* Loop unrolling
* Implicit casting of constants
* Anything that involves copy-pasting
* Common sub-expression elimination

Examples of **valid** optimizations include:

* Selection of an asymptotically-better algorithm
* Checks to skip execution of large blocks of code
* Refactoring array layouts for better cache performance
* Factorizations to avoid copying large objects
* Factorizations to avoid dynamic memory allocations

Limit the mental workload on programmers / maintainers to that work which compilers *cannot* do.
