Code Style
==========

Braces and indents
------------------

**Never** indent with tabs; use four spaces. (Almost all editors can be configured to insert spaces when the `tab` key is pressed).

K&R indent style, with an exception for the case of wrapped long lines (explained later). 

* Single space between braces and parens.
* Single space after `if` and `else`.
* No spaces on the inside of parens.
* Single space after commas.

Example:

    int some_function(int var, int var2, int var3) {
        if (var < 1) {
            do_code(var, var2);
        } else {
            do_something_else(var, var3);
        }
    }

Only a single-line block body can be without braces:

    if (test > 0) short_expr(x);
    
    while (*i) i++;

Naming style
------------

Class names are `UpperCamelCase`; variables, methods, and function names are `snake_case`:

    class ThingDoer {
        void do_thing();
    };
    
    int do_a_different_thing() {
        int thing_the_first  =  1;
        int thing_the_second = 99;
        // ...
    }

*Member typedefs* and *typedefs of POD data* types are `snake_case` ending with `_t`:

    typedef std::make_signed<size_t>::type index_t;
    
    template <typename T>
    class SomeContainer {
        typedef T elem_t;
        
        // ...
    }

Wrapping lines
--------------

It is best for lines to be 95 columns or fewer.

Long or complicatedly-nested function calls should break each top-level argument onto its own line, having an indent at least one block deeper than the beginning of the fuction name. It is preferable to wrap the line before the first argument (rather than indent more deeply to align with it).

    int x = some_fn_call(
                1,
                myVariable,
                anotherFunctionCall(5, z, 22),
                NUM_DOLPHINS);

Complicated function signatures follow this pattern: 

* Long argument lists are broken into one argument per line. 
* Wrapped arguments should be indented one block beyond the function body.
* Whenever a block-opening expression is line-wrapped, the opening brace goes on its own line.

<b></b>

    int some_function(
            int arg1,
            int arg2,
            char** thingies,
            float x,
            float y,
            SomeClass* myFoo)
    {
        // body goes here
        // ...
    }

Simpler function signatures which fit onto a single line should be kept that way:

    // simple function definition stays on one line:
    Matrix mult(const Matrix& a, const Matrix& b) {
        //...
    }

Complicated `if` expressions can be line-wrapped, and as with wrapped function signatures, the opening brace gets its own line:

    if (complicated_predicate(some_expression(a)) and 
        complicated_predicate(some_other_expr(b)))
    {
        // ...
    }

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

Enclose each branch of a ternary "if" in parentheses unless it is a single token or function call, and surround the `?` and `:` with spaces:

    int x = (a or b) ? (y + 1) : z;

If the sub-expressions are long, then each clause may be on its own (parenthesized) line, with the delimiting character on the beginning of the line:

    int x = (complicated_predicate_of(x, y) and z > 0)
                ? expression_if_true(x, z, x * y * z)
                : expression_if_false(x, x * x);

This makes it easier to read which sub-expressions belong to the true and false branches.

If a ternary "if" doesn't fit into three lines as above, it is usually better to use an ordinary "if/else" block statement.

Standard guidelines about wrapping and indenting apply.

Defines
-------

All `#define`s which are not single tokens must be parenthesized:

    #define NUM_DOLPHINS (INT_MAX / 2)

All macros must parenthesize their arguments as well as their entire body expression:

    #define FOO(a, b) ((a) * (b))

This is because unexpected things may happen once tokens are substituted:

    // xxx: bad style:
    #define FOO(a, b) a * b
    
    // expands to `206`, not `1000`!
    // this substitutes to `100 * 2 + 3 * 2`— not what you expected!
    100 * FOO(2 + 3, 2)

Line breaks
-----------

Double line breaks between class and functions; single line breaks to break up logical blocks within a function:

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

Pointers and references
-----------------------

Pointer and reference decorations are part of the type. Therefore, they should be visually grouped with the type. There is **no whitespace** between a `*` or `&` type decoration and the pointee type:

    // "*" and "&" are grouped with the type:
    int*      one;
    Vec<T,N>& two;
    
    void f(const AffineTransform<T,N>& xf, T* v) {
        ...
    }

The opposite is true for "dereference" and "addressof" operators: Those operate on the object, and should be grouped with its name:

    int* foo_p = &bar;

There is one case where this policy becomes misleading (and it is the most commonly cited reason for grouping pointer decorations with the variable instead of the type), and that is in the context of comma-separated variable declarations:

    int* a, b;  // xxx: incorrect (confusing) style!
   
This misleadingly declares as `a` as an `int*`, and `b` as an `int`. This is a mistake in the design of the C++'s operator associativity/precedence, but there is a simple solution: *do not declare variables this way.*

Variables should never be declared with comma separation syntax; each variable should be **declared on its own line:**

    // this is correct style:
    int* a;
    int* b;

Increments
----------

"Pythonic" increments are preferred, and increment operators should generally not be used:

    // preferable-- looks like an assignment:
    popsicles += 1;
    
    // discouraged:
    popsicles++;
    
    // also discouraged:
    put_in_freezer(popsicles++);

It is better if every assignment or change to a variable look like an assignment.

Loops preambles are one place where increment operators are allowed. In this case pre-increment is preferred;

    // loop variable pre-increment operators are ok:
    for (index_t i = 0; i < n; ++i) {
        // ...
    }

In the case where `i` is an iterator class, pre-increments may be slightly more efficient. (Post-increments need to make a copy of their pre-increment state to return). For this reason, it is best to use pre-increment consistently in all cases.


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

Initializer lists have one init per line. There is no space between the function definition and the inciting colon. Like long function definitions, inits are indented one deeper than the function body, and line wrapping bumps the opening brace to its own line:

    Foo::Foo(int a, int b):
            member_one(a),
            member_two(b)
    {
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

Where large temporary buffers are needed, consider allowing the user to provide them, or encapsulate the operation and its buffers in a class that can be re-used across multiple calculations.

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

Loop counters should be type `index_t`, which will be 64-bit on 64-bit platforms, allowing for traversal of large arrays:

    for (index_t i = 0; i < n; i++) {
        v[i] += k[i];
    }

Unlike `size_t`, `index_t` is signed, which will prevent warnings when used in with other signed numbers in arithmetic.

While it is true that some loops will not require the extra bit depth, it is better to maintain consistency, and the choice of `index_t` is unlikely to have any negative impact on performance.

Argument passing style
----------------------

Objects that are modified by a function should be passed as pointers:
    
    void foo(Object* out, int in) {
        out->x = some_function(n);
    }

Objects that are passed by reference should be *always* declared const:

    int bar(const Object& obj, int k) {
        return (obj.x + k) / obj.z;
    }

The above makes it obvious at the call site whether a function will modify an argument:

    invert(&dst, src);

