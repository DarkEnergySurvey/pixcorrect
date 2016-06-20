"""Simple, limited Design by Contract (DbC) module for python

Design by Contract (DbC) is a methodology for software developement in
which each component has well defined expectations and obligations, 
enforced through testing preconditions, postconditions, and invariants.

To use the module, define python functions that enforce any assertions
required. Each enforcement function must accept the following keyword
arguments:

func_args 
  A python list listing the arguments being supplied to the function
  whose contract is being enforced.

func_kwargs
  A python dictionary with the keyword arguments being supplied to the
  function whose contract is being enforced. Be sure to handle the case 
  when the keyword arguments are not supplied.

func_result
  The return value of the function whose contract is being enforced.

For example, I can define a set of functions that enforce contracts
that specify that certain arguments, keyword values, and returns fall
in a certain range of values:

>>> def arg_in_range(min_value, max_value, arg_num, 
...                  func_args=None, func_kwargs=None, func_result=None):
...     msg_template = "Argument %d not in range (%f to %f)" 
...     msg = msg_template % (arg_num, min_value, max_value)
...     assert min_value < func_args[arg_num] < max_value, msg
... 
>>> def kw_in_range(min_value, max_value, kw,
...                 func_args=None, func_kwargs=None, func_result=None):
...     if kw in func_kwargs:
...         msg_template = "Value of keyword %s not in range (%f to %f)" 
...         msg = msg_template % (kw, min_value, max_value)
...         assert min_value < func_kwargs[kw] < max_value, msg
... 
>>> 
>>> def result_in_range(min_value, max_value,
...                     func_args=None, func_kwargs=None, func_result=None):
...     msg_template = "Result not in range (%f to %f)"
...     msg = msg_template % (min_value, max_value)
...     assert min_value < func_result < max_value, msg
... 
>>> 

I can then use the precondition and postcondition python decorators
supplied by this module to apply them to my function, and providing
the rang, argument, index, and keywords to be tested as parameters to
the decorators:

>>> @precondition(arg_in_range, 0, 10, 0)
... @precondition(kw_in_range, 1, 12, 'dx')
... @postcondition(result_in_range, 2, 9)
... def foo(x, dx=5):
...     return x+dx
... 
>>> 

If all contract conditions are met, the decorators allow the function
to be called normally:

>>> print foo(1)
6
>>> print foo(1, dx=4)
5

If the contract is violated, though, an assertion fails:

>>> print foo(-4)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "dbc.py", line 10, in func_wrapper
    check_f(*check_args, **check_kwargs)
  File "<stdin>", line 6, in arg_in_range
AssertionError: Argument 0 not in range (0.000000 to 10.000000)
>>> 
>>> print foo(8, dx=0)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "dbc.py", line 11, in func_wrapper
    result = f(*args, **kwargs) 
  File "dbc.py", line 10, in func_wrapper
    check_f(*check_args, **check_kwargs)
  File "<stdin>", line 7, in kw_in_range
AssertionError: Value of keyword dx not in range (1.000000 to 12.000000)
>>> 
>>> print foo(8, dx=4)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "dbc.py", line 11, in func_wrapper
    result = f(*args, **kwargs) 
  File "dbc.py", line 11, in func_wrapper
    result = f(*args, **kwargs) 
  File "dbc.py", line 14, in func_wrapper
    check_f(*check_args, **check_kwargs)
  File "<stdin>", line 5, in result_in_range
AssertionError: Result not in range (2.000000 to 9.000000)
>>> 
>>> 

The processess for enforcing invariants is similar:

>>> def keep_arg_len(argnum, func_args=None, func_kwargs=None):
...     n = len( func_args[argnum] )
...     return n
... 
>>> def keep_arg_max(argnum, func_args=None, func_kwargs=None):
...     m = max( func_args[argnum] )
...     return m
... 
>>> @invariant_condition(keep_arg_len, "arg length changed", 0)
... @invariant_condition(keep_arg_max, "max changed", 0)
... def foo(x, n, y):
...     if n>len(x)-1 :
...         x.append(y)
...     else:
...         x[n]=7
... 
>>> foo( [1, 9, 3, 2], 2, 5 )
>>> 
>>> foo( [1, 9, 3, 2], 1, 5 )
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "dbc.py", line 134, in func_wrapper
    result = f(*args, **kwargs) 
  File "dbc.py", line 136, in func_wrapper
    assert after==before, message
AssertionError: max changed
>>> 
>>> foo( [1, 9, 3, 2], 6, 5 )
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "dbc.py", line 136, in func_wrapper
    assert after==before, message
AssertionError: arg length changed
>>> 
>>> 


"""

from functools import wraps, partial

def contract_condition(pre, post, check_f, *check_args, **check_kwargs):
    def make_func_wrapper(f):
        @wraps(f)
        def func_wrapper(*args, **kwargs):
            check_kwargs['func_args'] = args
            check_kwargs['func_kwargs'] = kwargs
            if pre:
                check_f(*check_args, **check_kwargs)
            result = f(*args, **kwargs) 
            if post:
                check_kwargs['func_result'] = result
                check_f(*check_args, **check_kwargs)
            return result
        return func_wrapper
    return make_func_wrapper

precondition = partial(contract_condition, True, False)
postcondition = partial(contract_condition, False, True)

def invariant_condition(check_f, message, *check_args, **check_kwargs):
    # The message parameter is included here but not in pre or post condiditons,
    # because there is no good way to set the message from within check_f in the
    # case of invariants, but it is the most covenvient way for pre and post-conditions.
    
    def make_func_wrapper(f):
        @wraps(f)
        def func_wrapper(*args, **kwargs):
            check_kwargs['func_args'] = args
            check_kwargs['func_kwargs'] = kwargs
            before = check_f(*check_args, **check_kwargs)
            result = f(*args, **kwargs) 
            after = check_f(*check_args, **check_kwargs)
            assert after==before, message
            return result
        return func_wrapper
    return make_func_wrapper

def make_args_check(f):
    """Create a function for a dbc decorator from a checking function

    :Parameters:
        -`f`: A function to check a parameter

    @returns: A func that takes a func as a parameter, and returns a function suitable for application with a dbc decorator.

    """
    # Convert a checker of one arg into a function that can be
    # passed to dbc.precondition or dbc.postcondition to check
    # each argument.
    def ff(arglist, func_args=None, func_kwargs=None, func_result=None):
        if not hasattr(arglist, '__iter__'):
            arglist = [arglist]

        for i in arglist:
            f(func_args[i])
    return ff
