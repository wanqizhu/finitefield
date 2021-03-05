from FiniteField import FiniteField
from FiniteFieldElem import FiniteFieldElem


F = FiniteField(5)
print(F(3))
print(F([2]))
# print(F([3, 4]))  # error 

F2 = FiniteField(2, 2, [1, 1], [1, 0])