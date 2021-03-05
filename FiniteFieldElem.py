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
        self.field = field

        if isinstance(val, int):
            val = [val]

        if len(val) != self.field.m:
            raise ValueError(f"Expected {self.field.m} coefficients, got "
                             + f"{len(val)} ({val}).")

        self.coefs = val


    """ For printing. Outputs string of the form
    [F_32] (1, 0, 0, 1, 1)
    """
    def __str__(self):
        return f"[F_{self.field.q}] ({', '.join(map(str, self.coefs))})"


    """ Allows dictionary lookup """
    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self.field.q == other.field.q and self.coefs == other.coefs

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
