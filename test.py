from FiniteField import FiniteField
from FiniteFieldElem import FiniteFieldElem
from ReedSolomon import ReedSolomonCode
from utils import gaussian_elimination, solve_lin_sys, poly_div
import numpy as np
import unittest

class TestFiniteFields(unittest.TestCase):
    def test_field_construction(self):
        F = FiniteField(5)
        self.assertEqual(F(3).coefs, [3])
        self.assertEqual(F(-3), F([7]))  # check mod works, check both int and
                                         # list constructors work for field elem
        with self.assertRaises(ValueError):
            F([3, 4])  # too many coefficients


    """
    Testing m > 1, no primitive element
    1 + x + x^2 is a primitive polynomial over GF(2)
    """
    def test_nonprime_field_f4_arith(self):
        F4 = FiniteField(2, 2, [1, 1, 1])
        one_plus_x = F4([1, 3])
        x = F4([6, -1])
        one = F4([1, 0])
        self.assertTrue(one != x)
        self.assertEqual(one_plus_x + one, x)
        self.assertEqual(one, F4([115, -28]))

        self.assertEqual(one_plus_x * x, one)
        self.assertEqual(one_plus_x * one_plus_x, x)


    """
    Testing primitive element finder for m > 1
    1 + x^2 is a prime polynomial for any p == 3 mod 4
    """
    def test_prim_elem_generation_f49(self):
        F49 = FiniteField(7, 2, [1, 0, 1])
        self.assertEqual(F49([3, 4]) * F49([1, 3]),  # multiply the old fashioned way
                         F49([5, 6]))

        F49.find_primitive_elem()
        self.assertTrue(F49.primitive_elem)
        self.assertEqual(F49([3, 4]) * F49([1, 3]),  # multiply fast: should be same result
                         F49([5, 6]))
        self.assertEqual(F49([1, 3]) * F49([4, 2]),
                         F49(5))  # (1 + 3x) * (4 + 2x) == 5


class TestMathUtils(unittest.TestCase):
    
    """
    Testing gaussian elim with regular integers
    """
    def test_gaussian_elim(self):
        A = [[2, -1, 0, 1, 0, 0],
             [-1, 2, -1, 0, 1, 0],
             [0, -1, 2, 0, 0, 1]]

        A_ = gaussian_elimination(A)
        self.assertTrue(np.allclose(A_, 
                       [[1, 0, 0, 3/4, 1/2, 1/4],
                       [0, 1, 0, 1/2, 1, 1/2],
                       [0, 0, 1, 1/4, 1/2, 3/4]]))


        A = [[2, 1, -1, 8],
            [-3, -1, 2, -11],
            [-2, 1, 2, -3]]

        A_ = gaussian_elimination(A)
        self.assertTrue(np.allclose(A_, 
                           [[1, 0, 0, 2],
                           [0, 1, 0, 3],
                           [0, 0, 1, -1]]))


    """
    Same operations, over modular field
    """
    def test_gaussian_elim_f5(self):
        F = FiniteField(5, find_primitive_elem = True)
        A = [[2, -1, 0, 3],
             [-1, 2, -1, 0],
             [0, -1, 2, 2]]

        for i in range(len(A)):
            A[i] = [F(x) for x in A[i]]

        A_ = gaussian_elimination(A)
        A_sol = [[1, 0, 0, 4],
                 [0, 1, 0, 0],
                 [0, 0, 1, 1]]

        for i in range(len(A)):
            for j in range(len(A[i])):
                # directyl compare with int, since field is of prime order
                self.assertEqual(A_[i][j], A_sol[i][j])


    # TODO: test above on fields of higher order, using polynomial representations

    # TODO: test solve_lin_sys

    '''
    5 + 4x + 2x^2 + 3x^3 = (2 + 1x)(12 - 4x + 3x^2) + (-19)
    '''
    def test_poly_div(self):
        f = [5, 4, 2, 3]
        g = [2, 1]
        q, r = poly_div(f, g)
        self.assertEqual(q, [12, -4, 3])
        self.assertEqual(r[0], -19)



class TestRSCode(unittest.TestCase):
    def test_F113_max_error(self):
        F = FiniteField(113, find_primitive_elem = True)
        eval_points = list(map(F, range(1, 17)))
        n = 16
        k = 8

        RS = ReedSolomonCode(F, n, k, eval_points)

        m3 = list(map(F, [3, 4, 5, 6, 7, 8, 9, 10]))
        c3 = RS.encode(m3)
        c3_decoded = RS.decode(c3)
        self.assertEqual(m3, c3_decoded)
        # add random errors, up to the max decoding capacity (4)
        c3[3] += F(2)
        c3[6] += F(13)
        c3[7] += F(21)
        c3[11] += F(88)
        c3_err_decoded = RS.decode(c3)
        self.assertEqual(m3, c3_err_decoded)



    def test_F7_zero_msg(self):
        F = FiniteField(7, find_primitive_elem = True)
        eval_points = list(map(F, range(1, 5)))
        n = 4
        k = 2

        RS = ReedSolomonCode(F, n, k, eval_points)
        m1 = list(map(F, [0, 0]))
        c1 = RS.encode(m1)
        c1_decoded = RS.decode(c1)
        self.assertEqual(m1, c1_decoded)
        # add random error
        c1[3] = F(2)
        c1_err_decoded = RS.decode(c1)
        self.assertEqual(m1, c1_err_decoded)
        c1[3] = F(0)
        c1[1] = F(4)
        c1_err_decoded = RS.decode(c1)
        self.assertEqual(m1, c1_err_decoded)
        


if __name__ == '__main__':
    unittest.main(exit=False)