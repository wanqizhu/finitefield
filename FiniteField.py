from typing import Union
from FiniteFieldElem import FiniteFieldElem

class FiniteField:
    """
    Represents a Finite Field of order q = p^m, where p is the prime 
    characteristic of the field. The elements of this field are polynomials 
    in F_p[x] mod g(x), where g is the degree-m primitive polynomial used to 
    construct the field.

    Raises ValueError if incorrect information is specified for the field.

    @param p:
        prime characteristic of the field
    @param m:
        log_p of the size of the field [default 1]
    @param primitive_poly:
        coefficients of the degree m primitive polynomial [ignored if m=1,
        necessary if m > 1]
    @param primitive_elem:
        a known primitive element of the field [optional]
    """
    def __init__(self, p: int, m: int = 1, 
                 primitive_poly: list[int] = None, 
                 primitive_elem: list[int] = None):
        # TODO: verify p is prime
        if not isinstance(m, int):
            raise ValueError(f"m must be an integer (got {m})")

        self.p = p
        self.q = p ** m
        self.m = m
        self.primitive_poly = primitive_poly
        self.primitive_elem = primitive_elem
        self.log_table = {}
        self.inverse_log_table = {}

        if primitive_elem:
            self.build_log_table()


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
        if self.m > 1 and (isinstance(value, int) or isinstance(value, float)):
            if not self.primitive_elem:
                self.find_primitive_elem()

            if value not in self.log_table():
                raise ValueError("Invalid representation of field element: " +
                                 f"{value}")

            return self.log_table[value]

        # value must represent coefficients; pass directly to FiniteFieldElem
        # constructor, which chcecks for errors
        return FiniteFieldElem(self, value)


    """
    Finds a primitive element and construct the multiplication log table.
    If a primitive element already exists, do nothing.
    
    This function searches randomly, so it may run for a long time if the field
    order is large.

    @return: the primitive element of this field
    """
    def find_primitive_elem(self):
        pass
        self.build_log_table()
        pass


    """
    Builds a multiplication log table using the primitive element.

    If no primitive element exists, calls `find_primitive_elem()` and uses
    the result to bulid the table.
    """
    def build_log_table(self):
        if not self.primitive_elem:
            self.find_primitive_elem()

        self.log_table = {}
        self.inverse_log_table = {}

        current_multiple = FiniteFieldElem(self, [0] * (self.m - 1) + [1])
        multiplier = FiniteFieldElem(self, self.primitive_elem)

        for i in range(self.q - 2):
            self.log_table[i] = current_multiple
            self.inverse_log_table[current_multiple] = i
            current_multiple *= multiplier

        # we represent zero as alpha ^ -inf, since then when we multiply
        # anything with zero, after taking logs we add -inf with some integer
        # and still get -inf, so we can reverse the log and get back zero

        # python doesn't have integer infinities, so we use float here
        zero = FiniteFieldElem(self, [0] * self.m)
        self.log_table[float('-inf')] = zero
        self.inverse_log_table[zero] = float('-inf')