# -*- coding: utf-8 -*-
def lproduct(*sets,**kwargs):
    """
    Cartesian product of input iterables by Steven Taschuk
    
    Becomes obsolute for Python 2.6, since it is incorporated in itertools
    """
    # http://code.activestate.com/recipes/159975/#c2
    wheels = map(iter, sets) # wheels like in an odometer
    digits = [it.next() for it in wheels]
    digits_n = len(digits)
    while True:
        yield digits[:]
        for i in xrange(digits_n-1, -1, -1):
            try:
                digits[i] = wheels[i].next()
                break
            except StopIteration:
                wheels[i] = iter(sets[i])
                digits[i] = wheels[i].next()
        else:
            break