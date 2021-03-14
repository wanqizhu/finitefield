from FiniteField import FiniteField, FiniteFieldElem
from utils import gaussian_elimination, solve_lin_sys, poly_div, poly_eval
import sys
import warnings

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

        if n > self.field.q:
            raise ValueError("Cannot have n > q in Reed Solomon code.")

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
            return sum([x**i * plaintext[i] for i in range(self.k)])

        return [f(a) for a in self.eval_points]


    """ Uses BerkekampWelch algorithm to decode a corrupted
    codeword with up to (n-k)/2 errors. """
    def decode(self, ciphertext: list[FiniteFieldElem]):
        max_num_errors = (self.n - self.k)//2

        for e in range(max_num_errors, -1, -1):
            # Try to find a monic polynomial E() of degree e
            # and a polynomial Q() of degree e+k-1 such that
            # w_i * E(a_i) = Q(a_i) for 0 <= i < m, where m = 2*e+k, w_i
            # is the ciphertext message, and a_i are the evaluation points.
            # This is equivalent to solving a system of equations Ax = b,
            # where x = (e_0, e_1, ..., e_{e-1}, q_0, q_1, ..., q_{e+k-1})
            # are the coefficients of the polynomials we're trying to find.
            # Each row of Ax computes
            #   w_i * (E(a_i) - a_i^e) - Q(a_i)
            # where we used the assumption that E has leading coefficient 1,
            # and we check whether that is equal to the value in b, 
            #   - w_i * a_i^e

            # If the system is over-determined, then we try again with a smaller e.
            # If the actual msg has <= max_num_errors from some valid codeword,
            # then this algorithm will always succeed for some e.

            def gen_row(i):
                coef_E = [ciphertext[i] * self.eval_points[i]**j for j in range(e)]
                coef_Q = [-self.eval_points[i]**j for j in range(e + self.k)]
                return coef_E + coef_Q

            m = 2 * e + self.k
            A = [gen_row(i) for i in range(m)]
            b = [-ciphertext[i] * self.eval_points[i]**e for i in range(m)]

            try:
                x = solve_lin_sys(A, b)
            except ValueError:
                warnings.warn(f"No solution found for e={e}, trying again...")
                continue

            E = x[:e] + [self.field(1)]
            Q = x[e:]

            f, r = poly_div(Q, E)
            if len(r) > 0:
                raise ValueError("Unknown Error: found Q, E but E did not "
                                 + " divide Q. "
                                 + "Q: {Q}, E: {E}, f: {f}, r: {r}")

            return f


        raise ValueError("No decoding found!")