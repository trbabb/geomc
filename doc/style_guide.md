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

Blank lines
-----------

It is more tidy if blank lines are indented to the same level as surrounding code; most editors can be configured to do this automatically (or at least not delete whitespace from blank lines). 

    void foo() {
    ····index_t a = // ...
    ····index_t b = // ...
    ····
    ····for (index_t i = 0; i < n; ++i) {
    ········// ...
    ····}
    }

This makes it a little nicer when adding a new comment or line of code into the blank space, especially on double blank lines between functions, as some editors only auto-indent based on the line above.

Alignment
-----------

It is good to highlight parallel structure by aligning similar elements vertically:

    int myvar_1   =  func(1,   2, 3);
    int myvar_123 =  func(1,  55, 0);
    int myvar_2   =  func(2, 255, 0);
    int myvar_3   = thing(1,   0, 0);

This makes it easy to examine similarities and differences, reducing mental workload, and making errors much easier to find. Here's an example:

    if (a_n->n_items != b_n->n_items) return false;
    if (a_n->n_children != a_n->n_children) return false;
    if (a_n->data != b_n->data) return false;

With vertical alignment, the code is much easier to read. The alignment also highlights that `b_n` has been mistyped as `a_n` on the second row:

    if (a_n->n_items    != b_n->n_items)    return false;
    if (a_n->n_children != a_n->n_children) return false;
    if (a_n->data       != b_n->data)       return false;

When in doubt about whether to justify left or right, prefer to keep digit places and matching text vertically aligned. Do not displace the beginning of a line, as it destroys the indentation.

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

Logical operators
-----------------

The `and`, `or`, and `not` keywords are preferred over `&&`, `||` and `!`:

    if ((foo == 1 and bar == 2) or baz == 3) {
        return;
    }
    
    if (not qux) {
        return;
    }

These keywords are completely portable and standards-compliant, and have been [since before QWERTY keyboards were standard](http://stackoverflow.com/questions/2393673/c-and-or-not-xor-keywords). Most syntax highlighters will correctly identify them as keywords and highlight them.

Ternary "if":
-------------

Enclose each branch of a ternary "if" in parentheses unless it is a single token, and surround the `?` and `:` with spaces:

    int x = (a or b) ? (y + 1) : z;


Wrapping lines
--------------

Long or complicatedly-nested function calls should break each top-level argument onto its own line, having an indent at least one block deeper than beginning of the fuction name. It is preferable to wrap the line before the first argument (rather than indent more deeply to align with it).

    int x = someFunctionCall(
                1,
                myVariable,
                anotherFunctionCall(5, z, 22),
                K_NUM_DOLPHINS);

Complicated function signatures should obey the same rule: Long argument lists are broken into one per line, though closely-related sets of arguments may share a line. Wrapped arguments should be indented one block beyond the function body:

    int someFunction(
            int arg1,
            int arg2,
            char** thingies,
            float x, float y,    // x, y are a pair
            SomeClass* myFoo) {
        // body goes here
        // ...
    }

Simpler function signatures which fit onto a single line should be kept that way:

    // simple function definition stays on one line:
    Matrix mult(const Matrix& a, const Matrix& b) {
        //...
    }

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

Initializer lists have one init per line. There is no space between the function definition and the inciting colon. Like long function definitions, inits are indented one deeper than the function body:

    Foo::Foo(int a, int b):
            member_one(a),
            member_two(b) {
        // body goes here
        // ...
    }

Comments
--------

A single space follows comment delimiters:

    // I have no idea what's happening here

If a line or block of code required any mental work or reasoning beyond trivial 
idioms, a comment is required there. If there is any "why was it done this way" 
or "what does this mean" question that is not trivially explained by the code 
itself without mental effort, a comment is required there. The goal is that a 
future reader-- including yourself-- should never have to reproduce that 
reasoning from scratch or from reverse engineering the code.

Comment as though you'll have all your memory of working on the code wiped
before the next time you'll touch it. This is closer to the truth than many
people tend to admit, and is also a good method for empathizing with future
maintainers who are not yourself!

Programming
===========

Memory allocation
-----------------

Assume wherever possible that code will be called in an inner loop. Consider heap memory allocations to be expensive, and prefer using stack variables and (small) arrays over dynamic allocations wherever it is reasonable. 

Where large temporary buffers are needed, consider allowing the user to provide them, or encapsulate the operation and its buffers in a class that can be re-used.

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
    
There is one (uncommon) exception, which is when choosing between two strategies that each return a more accurate result depending on scale; for example, a complete formula vs. its taylor series approximation near an unstable point. In these cases the numerical considerations should be well understood, and the two strategies should perform demonstrably better in their respective domains.

Loops
------

Loop variables should be type `index_t`, which will be 64-bit on 64-bit platforms, allowing for traversal of large arrays:

    for (index_t i = 0; i < n; i++) {
        v[i] += k[i];
    }

While it is true that some loops will not require the extra bit depth, it is better to maintain consistency, and the choice of `index_t` is unlikely to have any negative impact on performance.

Argument passing style
----------------------

Objects that are modified by a function should be passed as pointers:
    
    void foo(Object *out, int in) {
        out->x = some_function(n);
    }

Objects that are passed by reference should be *always* declared const:

    int bar(const Object &obj, int k) {
        return (obj.x + k) / obj.z;
    }

The above makes it obvious at the call site whether a function will modify an argument:

    invert(&dst, src);

