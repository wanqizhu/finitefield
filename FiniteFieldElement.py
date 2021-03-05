"""
Each instance of this class represents an element of the Finite Field `field`.

Most arithmetic operators are overloadded for this class, so that e.g.

```
F = FiniteField(7)
a = FiniteFieldElem(F, 3)
b = FiniteFieldElem(F, 4)
a * b
```

will output 5, which is 3*4 mod 7.
"""
class FiniteFieldElem:
    """
    Construct a field element, represented by a list of polynomial coefficients.
    If field is of prime order (m = 1), then val can also be a single integer,
    instead of a list of size 1.
    """
    def __init__(self, field, val):
        pass

    """ For printing. Formats an int (m = 1) or a list of coefficients (m > 1).
    """
    def __repr__(self):
        pass


    """
    Arithemetic operations over this field.
    """
    def __add__(self, other):
        pass

    def __radd__(self, other):
        pass

    def __sub__(self, other):
        pass

    def __rsub__(self, other):
        pass

    def __mul__(self, other):
        pass

    def __rmul__(self, other):
        pass

    def __pow__(self, other):
        pass

    def __neg__(self):
        pass

    def __eq__(self, other):
        pass
