
class FiniteField:
    """
    Represents a Finite Field of order q = p^m, where p is the prime 
    characteristic of the field. The elements of this field are polynomials 
    in F_p[x] mod g(x), where g is the degree-m primitive polynomial used to 
    construct the field.

    Throws ValueError if incorrect information is specified for the field.

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
    def __init__(self, p, m=1, primitive_poly=None, primitive_elem=None):
        pass


    """
    Create a representation of an element of this field.

    A field element is usually represented by polynomial 
    coefficients (a list of size m). For m = 1, it's okay to pass in a single
    integer instead of a list of size 1.

    An alternative specification, for m > 1, 
    is by its log_alpha value (an int), where alpha is the primitive element 
    of this field. If no primitive element is specified for this
    field, then `find_primitive_elem` is called to find one. This primitive
    element will then be used for future operations in this field.

    @param value:
        a list of size m or an int, uniquely representing the element
    @return
        a field element
    """
    def __call__(self, value):
        pass


    """
    Finds a primitive element and construct the multiplication log table.
    If a primitive element already exists, do nothing.
    
    This function searches randomly, so it may run for a long time if the field
    order is large.

    @return: the primitive element of this field
    """
    def find_primitive_elem(self):
        pass