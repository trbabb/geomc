Code Style
==========

K&R indent style. 
* Single space between braces and parens. 
* Single space after `if` and `else`.
* No spaces on the inside of parens.
* Single space after commas.

    if (var < 1) {
        do_code(var, var2);
    } else {
        do_something_else(var, var3);
    }

Only single-line `if`s can be without braces:

    if (test > 0) do_code(x);

If there is a long block of assignments, line up the equal signs to make it lookl ke a more readable table:

    int      my_var = 1;
    int    some_var = 2;
    int another_var = 6;

Programming
===========

Do not allocate memory if it can be avoided. Stack variables and stack-allocated arrays are always p
referable where possible.

Do not use logic with "epislons" if it can be avoided, especially since most code is templated.
The scale of a sensible epsilon would be dependent on the template parameter type). Arithmetic
accuracy should fail naturally at machine precision.
