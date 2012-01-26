"""
Automatically parse command line args and kwargs.

This module is meant to facilitate the translation of command line arguments to
Python code to parse to methods.

Example usage:

Given an example minimalistic Python module 'example.py'

>>> def testfunc(a,b,calc='sum'):
>>> ... if calc=='sum': return a+b
>>> ... elif calc=='prod': return a*b
>>> ... return None

>>> if __name__=="__main__":
>>> ... method,args,kwargs = argkwargparser.parse()
>>> ... output = globals()[method](*args,**kwargs)

Then, in a terminal, you can do::

    $:> python example.py testfunc 3 4
    7
    $:> python example.py testfunc 3 4 calc=prod
    12

You can mix args and kwargs, they will be sorted by L{parse}. You can give lists
and more sophisticated input for args and kwargs, see L{parse}.

"""
import json
import sys
    
def parse(argv=None):
    """
    Command-line to method call arg processing.
    
        - positional args: a b -> method('a', 'b')
        - intifying args: a 123 -> method('a', 123)
        - json loading args: a '["pi", 3.14, null]' -> method('a', ['pi', 3.14, None])
        - keyword args: a foo=bar -> method('a', foo='bar')
        - using more of the above 1234 'extras=["r2"]'  -> method(1234, extras=["r2"])
    
    @param argv: Command line arg list. Defaults to `sys.argv`.
    @return: method-name, args, kwargs
    @rtype: string, list, dict
    """
    if argv is None:
        argv = sys.argv

    method_name, arg_strs = argv[1], argv[2:]
    args = []
    kwargs = {}
    for s in arg_strs:
        #-- keyword argument if '=' in string
        if s.count('=') == 1:
            key, value = s.split('=', 1)
            #-- maybe the user preceded the keyword with '--' out of habit. In
            #   that case remove it
            if key[:2]=='--':
                key = key[2:]
        #-- normal argument
        else:
            key, value = None, s
        try:
            value = json.loads(value) 
        except ValueError:
            pass
        if key:
            kwargs[key] = value
        else:
            args.append(value)
    return method_name, args, kwargs

def test(*args,**kwargs):
    print 'args',args
    print 'kwargs',kwargs

if __name__=="__main__":
    method,args,kwargs = parse()
    out = globals()[method](*args, **kwargs)
    sys.exit(out)