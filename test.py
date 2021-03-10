from FiniteField import FiniteField
from FiniteFieldElem import FiniteFieldElem
from utils import gaussian_elimination
import numpy as np

F = FiniteField(5)
print(F(3))
print(F([2]))
# print(F([3, 4]))  # error 


"""
Testing m > 1, no primitive element
1 + x + x^2 is a primitive polynomial over GF(2)
"""

F4 = FiniteField(2, 2, [1, 1, 1])
one_plus_x = F4([1, 3])
x = F4([6, -1])
one = F4([1, 0])
print(one_plus_x)
print(x)
print(one_plus_x + one)
assert one != x
assert one_plus_x + one == x
assert one == F4([115, -28])

print(f"{one_plus_x} * {x} = {one_plus_x * x} == {one}")
print(f"{one_plus_x} * {one_plus_x} = {one_plus_x * one_plus_x} == {x}")



"""
Testing primitive element finder for m > 1
1 + x^2 is a prime polynomial for any p == 3 mod 4
"""

F49 = FiniteField(7, 2, [1, 0, 1])
print(F49([3, 4]) * F49([1, 3]))  # multiply the old fashioned way

F49.find_primitive_elem()
print("primitive elem:", F49.primitive_elem)
print(F49([3, 4]) * F49([1, 3]))  # multiply fast: should be same result




"""
Testing gaussian elim
"""
A = [[2, -1, 0, 1, 0, 0],
     [-1, 2, -1, 0, 1, 0],
     [0, -1, 2, 0, 0, 1]]

A_, _ = gaussian_elimination(A)
assert np.allclose(A_, 
                   [[1, 0, 0, 3/4, 1/2, 1/4],
                   [0, 1, 0, 1/2, 1, 1/2],
                   [0, 0, 1, 1/4, 1/2, 3/4]]), (A, A_)


A = [[2, 1, -1, 8],
    [-3, -1, 2, -11],
    [-2, 1, 2, -3]]

A_, _ = gaussian_elimination(A)
assert np.allclose(A_, 
                   [[1, 0, 0, 2],
                   [0, 1, 0, 3],
                   [0, 0, 1, -1]]), (A, A_)
