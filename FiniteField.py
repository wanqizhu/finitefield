from typing import Union
import itertools
from FiniteFieldElem import FiniteFieldElem


class FiniteField:
    """
    Represents a Finite Field of order q = p^m, where p is the prime 
    characteristic of the field. The elements of this field are polynomials 
    in F_p[x] mod g(x), where g is the degree-m prime polynomial used to 
    construct the field.

    Raises ValueError if incorrect information is specified for the field.

    @param p:
        prime characteristic of the field
    @param m:
        log_p of the size of the field [default 1]
    @param primitive_poly:
        coefficients of the degree m prime polynomial [ignored if m=1,
        necessary if m > 1]
    @param primitive_elem:
        a known primitive element of the field [optional]


    Attributes:

    log_table
        mapping from field elements to int; Zero maps to -inf
        log_table[elem] = i means that primitive_elem ^ i = elem 
    _inverse_log_table
        inverse mapping, from integer powers of primitive element to field elem

    """
    def __init__(self, p: int, m: int = 1, 
                 primitive_poly: list[int] = None, 
                 primitive_elem: list[int] = None):
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
        self.primitive_elem = primitive_elem
        self._log_table = {}
        self._inverse_log_table = {}
        self._has_log_table = False

        if primitive_elem:
            self.build_log_table()

    def has_log_table(self):
        return self._has_log_table

    def inverse_log(self, v):
        if not self.has_log_table():
            raise ValueError("Field does not have a log table. Please call "
                             + "self.build_log_table() first.")

        if isinstance(v, int):
            return self._inverse_log_table[v % (self.q - 1)]

        return self._inverse_log_table[v]  # special case for -inf -> Zero

    def log(self, elem):
        if not self.has_log_table():
            raise ValueError("Field does not have a log table. Please call "
                             + "self.build_log_table() first.")

        if elem not in self._log_table:
            raise ValueError("Invalid field element, cannot compute log.")
        return self._log_table[elem]


    """
    Create a representation of an element of this field.

    A field element is usually represented by polynomial 
    coefficients (a list of size m). For m = 1, it's okay to pass in a single
    integer instead of a list of size 1.

    An alternative specification, for m > 1, 
    is by its log_alpha value (an int), where alpha is the primitive element 
    of this field. If no primitive element is specified for this
    field, then `find_primitive_elem` is called to find one. This primitive
    element will then be used for future operations in this field. Note that
    0 is represented here by -inf.

    Raises ValueError if `value` does not match the specification above.

    @param value:
        a list of size m, an int, or -inf, uniquely representing the element
    @return
        a field element
    """
    def __call__(self, value: Union[list[int], int, float]):
        try:
            if self.m > 1 and (isinstance(value, int) or isinstance(value, float)):
                if not self.primitive_elem:
                    self.find_primitive_elem()

                return self.inverse_log(value)


            # value must represent coefficients; pass directly to FiniteFieldElem
            # constructor, which chcecks for errors
            return FiniteFieldElem(self, value)

        except ValueError:
            raise ValueError(f"Invalid representation of field element: {value}.")


    """
    Finds a primitive element and construct the multiplication log table.
    If a primitive element already exists, do nothing.
    
    This function searches randomly, so it may run for a long time if the field
    order is large.

    @return: the primitive element of this field
    """
    def find_primitive_elem(self):
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

        self.build_log_table()


    """
    Builds a multiplication log table using the primitive element.

    If no primitive element exists, calls `find_primitive_elem()` and uses
    the result to bulid the table.
    """
    def build_log_table(self):
        if not self.primitive_elem:
            self.find_primitive_elem()

        self._log_table = {}
        self._inverse_log_table = {}
        self._has_log_table = False

        current_multiple = FiniteFieldElem(self, 1)

        for i in range(self.q - 2):
            self._log_table[current_multiple] = i
            self._inverse_log_table[i] = current_multiple
            current_multiple *= self.primitive_elem

        # we represent zero as alpha ^ -inf, since then when we multiply
        # anything with zero, after taking logs we add -inf with some integer
        # and still get -inf, so we can reverse the log and get back zero

        # python doesn't have integer infinities, so we use float here
        zero = FiniteFieldElem(self, 0)
        self._log_table[zero] = float('-inf')
        self._inverse_log_table[float('-inf')] = zero
        self._has_log_table = True