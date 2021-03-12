from __future__ import annotations  # allow type hints of owning class
                                    # remove for python 3.10
from typing import Union
import itertools



class FiniteField:
    """
    Represents a Finite Field of order q = p^m, where p is the prime 
    characteristic of the field. The elements of this field are polynomials 
    in F_p[x] mod g(x), where g is the degree-m prime polynomial used to 
    construct the field.

    Raises ValueError if incorrect information is specified for the field.

    @param p
        Prime characteristic of the field
    @param m
        Log_p of the size of the field [default 1]
    @param primitive_poly
        Coefficients of the degree m prime polynomial [ignored if m=1,
        necessary if m > 1]
    @param find_primite_elem
        Whether to find a primitive element of this field.
        This makes multiplication faster and enables division. [default True]

    Attributes:

    _log_table
        Mapping from field elements to int; Zero maps to -inf
        log_table[elem] = i means that primitive_elem ^ i = elem 
    _inverse_log_table
        Inverse mapping, from integer powers of primitive element to field elem
    _has_log_table
        True if the log table & its inverse are constructed

    """
    def __init__(self, p: int, m: int = 1, 
                 primitive_poly: list[int] = None,
                 find_primitive_elem: bool = True):
        # TODO: verify p is prime
        if not isinstance(m, int):
            raise ValueError(f"m must be an integer (got {m})")

        if m > 1:
            try:
                if len(primitive_poly) != m + 1:
                    raise ValueError(f"primitive_poly must have length {m+1} "
                                     + f"(got {len(primitive_poly)}).")

                if primitive_poly[m] != 1:
                    raise ValueError(f"primitive_poly must be monic.")

                # TODO: check this poly is prime?
            except TypeError:
                raise

        self.p = p
        self.q = p ** m
        self.m = m
        self.primitive_poly = primitive_poly
        self.primitive_elem = None
        self._log_table = {}
        self._inverse_log_table = {}
        self._has_log_table = False

        if find_primitive_elem:
            self.find_primitive_elem()

    
    """
    Create a representation of an element of this field.

    A field element is usually represented by polynomial 
    coefficients (a list of size m), starting from the lowest ordered term.
    For field elements which do not have higher ordered terms 
    (e.g. all elements when m == 1), it's okay to pass 
    in a single integer instead of a list of size 1.

    Raises ValueError if `value` does not match the specification above.

    @param value:
        A list of size m, or an int uniquely representing the element
    @return
        A field element
    """
    def __call__(self, value: Union[list[int], int, float]):
        try:
            return FiniteFieldElem(self, value)

        except ValueError:
            raise ValueError(f"Invalid representation of field element: {value}.")


    """
    Finds a primitive element and construct the multiplication log table.
    If a primitive element already exists, do nothing.
    
    This function searches iteratively over all field elements, so it may run 
    for O(q) time, where q is the size of the field.

    @return: The primitive element of this field
    """
    def find_primitive_elem(self):
        if self.primitive_elem:
            return self.primitive_elem

        order_of_elems = {}
        for coefs in itertools.product(range(self.p), repeat = self.m):
            elem = FiniteFieldElem(self, coefs)
            if elem == 0:
                continue

            if elem in order_of_elems:
                continue

            generated_elems = [elem]
            while generated_elems[-1] * elem != elem:
                generated_elems.append(generated_elems[-1] * elem)

            order = len(generated_elems)
            if order == self.q - 1:
                # found primitive elem
                self.primitive_elem = elem
                break

            # cache all elements cycled through to save time
            for elem in generated_elems:
                order_of_elems[elem] = order

        self._build_log_table()
        return self.primitive_elem

    """ Getters for private vars """
    def has_log_table(self):
        return self._has_log_table

    def inverse_log(self, v: int):
        if not self.has_log_table():
            raise ValueError("Field does not have a log table. Please call "
                             + "self.build_log_table() first.")

        if isinstance(v, int):
            return self._inverse_log_table[v % (self.q - 1)]

        return self._inverse_log_table[v]  # special case for -inf -> Zero

    def log(self, elem: FiniteFieldElem):
        if not self.has_log_table():
            raise ValueError("Field does not have a log table. Please call "
                             + "self.build_log_table() first.")

        if elem not in self._log_table:
            raise ValueError("Invalid field element, cannot compute log.")
        return self._log_table[elem]


    #### 
    # Private functions
    ####

    """
    Builds a multiplication log table using the primitive element.

    Raises ValueError if called without a primitive element.
    """
    def _build_log_table(self):
        if not self.primitive_elem:
            raise ValueError("Internal function 'build_log_table' called without"
                             + " first finding a primitive element.")

        self._log_table = {}
        self._inverse_log_table = {}
        self._has_log_table = False

        current_multiple = FiniteFieldElem(self, 1)

        for i in range(self.q - 1):
            self._log_table[current_multiple] = i
            self._inverse_log_table[i] = current_multiple
            current_multiple *= self.primitive_elem

        # we represent zero as alpha ^ -inf, since when we multiply
        # anything with zero, after taking logs we add -inf with some integer
        # and still get -inf, so we can reverse the log and get back zero

        # python doesn't have integer infinities, so we use float here
        zero = FiniteFieldElem(self, 0)
        self._log_table[zero] = float('-inf')
        self._inverse_log_table[float('-inf')] = zero
        self._has_log_table = True



"""
Each instance of this class represents an element of the Finite Field `field`.

Most arithmetic operators are overloadded for instances of this class, so that
e.g.

```
F = FiniteField(7)
a = F(3)  # equivalent to FiniteFieldElem(F, 3)
b = F(4)  # equivalent to FiniteFieldElem(F, 4)
a * b
```

will output 5, which is 3*4 mod 7.

This class should not be constructed directly. Instead, use the overloaded
__call__() function of the field to build an instance of this class, as per the
example above.
"""
class FiniteFieldElem:
    """
    Construct a field element, represented by a list of polynomial coefficients
    which are remainder polynomials in F_p[x], ordered from the lowest power
    of x (x^0) to the highest (x^(m-1)).

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
    def __repr__(self):
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
    Binary arithemetic operations
    Addition and subtraction are only defined over elements of the same field.
    These operations work component-wise, the same way you would add or
    subtract vectors.
    Multiplication ideally works with a primitive element,
    in which case a fast-lookup via log tables are performed. When there's no
    primitive element, multiplication can be performed directly by multiplying
    the remainder polynomials and reducing mod g(x) via long division, where
    g(x) is the prime polynomial of the underlying field. Division only works
    with a primitive element.

    Scalar operations are equivalent to operating with the field element
    corresponding to that integer Note that taking an inverse is equal to
    exponentiation by -1. This works out just fine via the log table.
    """
    def check_arith_compatibility(self, other, op_str):
        if isinstance(other, int):
            other = self.field(other)

        if not isinstance(other, FiniteFieldElem):
            raise TypeError(f"Cannot {op_str} FiniteFieldElem with object of type"
                + f"{type(other)}.")

        if self.field.q != other.field.q:
            raise ValueError(f"Cannot {op_str} two elements from fields of different "
                             + f"sizes ({self.field.q} and {other.field.q}).")

        return other

    def __add__(self, other: Union[int, FiniteFieldElem]):
        other = self.check_arith_compatibility(other, "add")
        new_coefs = [x + y % self.field.p for x, y in zip(self.coefs, other.coefs)]
        return FiniteFieldElem(self.field, new_coefs)

    def __radd__(self, other: Union[int, FiniteFieldElem]):
        return self.__add__(other)

    def __sub__(self, other: Union[int, FiniteFieldElem]):
        other = self.check_arith_compatibility(other, "subtract")
        new_coefs = [x - y % self.field.p for x, y in zip(self.coefs, other.coefs)]
        return FiniteFieldElem(self.field, new_coefs)

    def __rsub__(self, other: Union[int, FiniteFieldElem]):
        return other + (self * -1)  ## TODO: make this better?

    def __mul__(self, other: Union[int, FiniteFieldElem]):
        other = self.check_arith_compatibility(other, "multiply")
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

    def __rmul__(self, other: Union[int, FiniteFieldElem]):
        return self.__mul__(other)

    def __truediv__(self, other):
        self.check_arith_compatibility(other, "divide")
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
        if not isinstance(exponent, int):
            raise TypeError(f"Cannot raise FiniteFieldElem to exponent "
                            + f"of type {type(exponent)}.")

        if self.field.has_log_table():
            # use log table
            log_self = self.field.log(self)
            return self.field.inverse_log(log_self * exponent)

        else:
            # compute power the hard way (not very efficient)
            # TODO: repeated squaring?
            # only works for positive exponents
            if exponent < 0:
                raise ValueError("Cannot compute inverse of FiniteFieldElem "
                                 + "over a field with no primitive elements. "
                                 + "Please call find_primitive_elem() first.")
            power = FiniteFieldElem(self.field, 1)
            for i in range(exponent):
                power *= self
            return power


    def __neg__(self):
        return FiniteFieldElem(self.field, [-c for c in self.coefs])

    def __bool__(self):
        return self != 0
