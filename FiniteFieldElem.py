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
    Construct a field element, represented by a list of polynomial coefficients
    which are remainder polynomials in F_p[x].
    If field is of prime order (m = 1), then val can also be a single integer,
    instead of a list of size 1.

    Attributes:
    @self.field: 
        FiniteField instance this element lives in
    @self.coefs
        list of ints of size self.field.m, ranging over [0..p-1]^m
    """
    def __init__(self, field, val):
        self.field = field

        if isinstance(val, int):
            val = [val]

        if len(val) != self.field.m:
            raise ValueError(f"Expected {self.field.m} coefficients, got "
                             + f"{len(val)} ({val}).")

        for i, v in enumerate(val):
            if not isinstance(v, int):
                raise TypeError(f"Got non-integer coefficient: {v}.")
            val[i] = v % self.field.p

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
    Arithemetic operations; only defined over elements of the same field
    """
    def check_arith_compatibility(self, other):
        if not isinstance(other, FiniteFieldElem):
            raise TypeError("Cannot add FiniteFieldElem with object of type"
                + f"{type(other)}.")

        if self.field.q != other.field.q:
            raise ValueError("Cannot add two elements from fields of different "
                             + f"sizes ({self.field.q} and {other.field.q}).")

    def __add__(self, other):
        self.check_arith_compatibility(other)
        new_coefs = [x + y % self.field.p for x, y in zip(self.coefs, other.coefs)]
        return FiniteFieldElem(self.field, new_coefs)


    def __sub__(self, other):
        self.check_arith_compatibility(other)
        new_coefs = [x - y % self.field.p for x, y in zip(self.coefs, other.coefs)]
        return FiniteFieldElem(self.field, new_coefs)

    def __mul__(self, other):
        pass

    def __div__(self, other):
        pass

    def __pow__(self, other):
        pass

    def __neg__(self):
        return FiniteFieldElem(self.field, [-c for c in self.coefs])
