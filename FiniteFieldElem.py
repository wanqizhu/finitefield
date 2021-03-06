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

    @param field
        FiniteField instance this element lives in
    @param val
        Iterable of ints of size self.field.m, ranging over [0..p-1]^m.
        If val is a single integer, then we interpret all the higher order
        coefficients as 0.
        
    """
    def __init__(self, field, val):
        self.field = field

        if isinstance(val, int):
            val = [val] + [0] * (self.field.m - 1)

        val = list(val)

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
        if isinstance(other, int):
            # eq iff we have no higher order x term and the x^0 coef match
            return self.coefs[0] == other and sum(self.coefs[1:]) == 0

        return self.field.q == other.field.q and self.coefs == other.coefs

    """
    Binary arithemetic operations; only defined over elements of the same field
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
        self.check_arith_compatibility(other)
        if self.field.has_log_table():
            # use log table
            log_self = self.field.log(self)
            log_other = self.field.log(other)
            log_product = log_self + log_other
            return self.field.inverse_log(log_product)

        else:
            # multiply and reduce directly
            # the coefficient for x^i in the product is equal to
            # the sum of a_j * a_(i-j)
            m = self.field.m
            p = self.field.p
            new_coefs = [sum([
                              self.coefs[j] * other.coefs[i-j] % p
                              for j in range(i+1)
                              if j < m and i-j < m
                             ])
                         for i in range(2 * m - 1)]

            # compute `new_coefs % modulus` via long division
            modulus = self.field.primitive_poly
            for top_term in range(2*m-2, m-1, -1):
                multiplier = new_coefs[top_term]
                bottom_term = top_term - m
                # zero out the coefficient of x^{top_term}
                # by subtracting c_{top_term} * modulus * x^{bottom_term}
                # since we know modulus[m] = 1
                for i in range(m+1):
                    new_coefs[bottom_term + i] -= multiplier * modulus[i]
                    new_coefs[bottom_term + i] %= p

            new_coefs_modded = new_coefs[:m]
            return FiniteFieldElem(self.field, new_coefs_modded)


    def __div__(self, other):
        self.check_arith_compatibility(other)
        if self.field.has_log_table():
            # use log table
            log_self = self.field.log(self)
            log_other = self.field.log(other)
            log_div = log_self - log_other
            return self.field.inverse_log(log_div)

        else:
            # hard to divide without inverse lookup table...
            raise ValueError("FiniteField does not support division without"
                             + "a primitive element. Please call "
                             + "f.find_primitive_elem() on the underlying "
                             + "field first.")


    def __pow__(self, exponent: int):
        if self.field.has_log_table():
            # use log table
            log_self = self.field.log(self)
            return self.field.inverse_log(log_self * exponent)

        else:
            # compute power the hard way (not very efficient)
            # TODO: repeated squaring?
            power = self
            for i in range(exponent-1):
                power *= self
            return power


    def __neg__(self):
        return FiniteFieldElem(self.field, [-c for c in self.coefs])
