from FiniteField import FiniteField
from FiniteFieldElem import FiniteFieldElem


F = FiniteField(5)
print(F(3))
print(F([2]))
# print(F([3, 4]))  # error 


"""
Testing m > 1
"""

F2 = FiniteField(2, 2, [1, 1])
print(F2([0, 1]))
print(F2([1, 3]))
print(F2([6, -1]))
print(F2([1, 3]) + F2([0, 1]))
assert F2([0, 1]) != F2([1, 0])
assert F2([1, 3]) + F2([0, 1]) == F2([1, 0])


"""
Testing primitive element initializer for m > 1
"""

# F2 = FiniteField(2, 2, [1, 1], [1, 0])
# print(F2[0])
# print(F2[1])
# print(F2[2])
# print(F2[6])
# print(F2[1] + F2[1])
