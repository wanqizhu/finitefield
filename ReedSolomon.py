from FiniteField import FiniteField
from FiniteFieldElem import FiniteFieldElem
from utils import gaussian_elimination, solve_lin_sys, poly_div, poly_eval
import sys


class ReedSolomonCode():
    """ Construct a specific code from the RS code family.
    """
    def __init__(self, field: FiniteField,
                 n: int,
                 k: int, 
                 eval_points: list[FiniteFieldElem]):

        if len(eval_points) != n:
            raise ValueError("Expected eval_points to be length "
                             + f"{n}, got {len(eval_points)}.")

        self.field = field
        self.n = n
        self.k = k
        self.eval_points = eval_points


    """ Maps plaintext in F_q^k to ciphertext in F_q^n. """
    def encode(self, plaintext: list[FiniteFieldElem]):
        if len(plaintext) != self.k:
            raise ValueError("Expected plaintext to be length "
                             + f"{self.k}, got {len(plaintext)}.")

        # Map our plaintext to the polynomial with matching coefficients
        def f(x: FiniteFieldElem):
            return sum([x**i * plaintext[i] for i in range(k)])

        return [f(a) for a in self.eval_points]


    """ Uses BerkekampWelch algorithm to decode a corrupted
    codeword with up to (n-k)/2 errors. """
    def decode(self, ciphertext):
        max_num_errors = (self.n - self.k)//2

        for e in range(max_num_errors, -1, -1):
            # Try to find a monic polynomial E() of degree e
            # and a polynomial Q() of degree e+k-1 such that
            # w_i * E(a_i) = Q(a_i) for 0 <= i < m, where m = 2*e+k.
            # This is equivalent to solving a system of equations Ax = b,
            # where x = (e_0, e_1, ..., e_{e-1}, q_0, q_1, ..., q_{e+k-1})
            # are the coefficients of the polynomials we're trying to find.
            # Notice that with how we constructed A and b below, the matrix
            # multiplication operation is equivalent to computing
            #   w_i * (E(a_i) - a_i^e) - Q(a_i)
            # where we used the assumption that E has leading coefficient 1,
            # and we check whether that is equal to
            #   - w_i * a_i^e

            # If the system is over-determined, then we try again with a smaller e.
            # If the actual msg has <= max_num_errors from some valid codeword,
            # then this algorithm will always succeed for some e.

            def gen_row(i):
                coef_E = [ciphertext[i] * self.eval_points[i]**j for j in range(e)]
                coef_Q = [-self.eval_points[i]**j for j in range(e+k)]
                return coef_E + coef_Q

            m = 2 * e + k
            A = [gen_row(i) for i in range(m)]
            b = [-ciphertext[i] * self.eval_points[i]**e for i in range(m)]

            try:
                x = solve_lin_sys(A, b)
            except ValueError:
                print(f"No solution found for e={e}, trying again...")
                print(sys.exc_info())
                continue

            # print("Lin System Solution:", x)

            E = x[:e] + [self.field(1)]
            Q = x[e:]

            f, r = poly_div(Q, E)
            if len(r) > 0:
                print("Error: found Q, E but E did not divide Q")
                print(f"Q: {Q}, E: {E}, f: {f}, r: {r}")
                continue

            print("Found polynomial:", f)

            codeword = [poly_eval(f, eval_point) 
                        for eval_point in self.eval_points]

            print("Decoded codeword:", codeword)
            print("Number of errors:", sum([c != m for c,m in zip(codeword, ciphertext)]))
            return codeword


        raise ValueError("No decoding found!")




F = FiniteField(113, find_primitive_elem = True)
m1 = list(map(F, [52, 84, 35, 108, 70, 78, 43, 109, 66, 20, 100, 103, 11, 41, 14, 70]))
m2 = list(map(F, [7, 64, 58, 10, 90, 89, 99, 54, 42, 55, 82, 24, 35, 95, 38, 25]))
m3 = list(map(F, [81, 81, 93, 60, 27, 12, 37, 72, 68, 58, 67, 4, 76, 105, 49, 35]))
eval_points = list(map(F, range(1, 17)))
n = 16
k = 8

RS = ReedSolomonCode(F, n, k, eval_points)

print(RS.decode(m1))
"""
Found polynomial: [12, 58, 5, 37, 58, 20, 79, 9]
Decoded codeword: [52, 84, 35, 108, 70, 78, 39, 37, 66, 20, 100, 103, 13, 20, 14, 70]
Number of errors: 4
"""
print(RS.decode(m2))
"""
No solution found for e=4, trying again...
Found polynomial: [93, 67, 48, 4, 105, 82, 13, 16]
Decoded codeword: [89, 64, 58, 33, 90, 89, 99, 54, 42, 55, 62, 24, 35, 95, 38, 25]
Number of errors: 3
"""
print(RS.decode(m3))
"""
Found polynomial: [110, 50, 12, 44, 32, 77, 112, 96]
Decoded codeword: [81, 103, 93, 60, 38, 12, 37, 72, 68, 58, 110, 4, 76, 39, 49, 35]
Number of errors: 4
"""
