from FiniteField import FiniteField
from FiniteFieldElem import FiniteFieldElem
from utils import gaussian_elimination


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
        pass
