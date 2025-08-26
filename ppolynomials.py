from sympy import poly, symbols
from sympy.polys.domains import GF
import itertools
from random import seed, sample
seed(1)
t = symbols('t')

p = 3
Fp = GF(p)
N = 1
LEN = 3

BASIS_Fpt = [t**i for i in range(p)]

def f7(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def get_random_A_subset(basis, n, len):
    Bn = f7((b**e for b, e in itertools.product(basis, range(p**n))))
    return list(sample(list(Bn),k=len)
)

def generate_ppolynomial(n, len):
    A = get_random_A_subset(BASIS_Fpt, n, len)
    print(A)
    return sum([a*symbols(f"x_{i}")**(p**n) for i, a in enumerate(A)])

if __name__ == "__main__":
    print(generate_ppolynomial(N, LEN))
