from sympy import poly, symbols, Poly
from sympy.polys.domains import GF
import itertools
from random import seed, sample
from typing import List
seed(1)
t = symbols('t')

p = 3
Fp = GF(p)
POWER = 1
NUM_OF_VAR = 3

BASIS_Fpt = [t**i for i in range(p)]

def f7(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def get_random_A_subset(basis: List, n, len):
    Bn = f7((b**e for b, e in itertools.product(basis, range(p**n))))
    return list(sample(list(Bn),k=len)
)

def generate_ppolynomial(n, len):
    A = get_random_A_subset(BASIS_Fpt, n, len)
    A = sorted(A, key=lambda p: Poly(p, t).degree())

    return sum([a*symbols(f"x_{i}")**(p**n) for i, a in enumerate(A)])

def generate_iteration(polynomial: Poly):
    """
    Take polynomial and return a system of equations by the algorithm for the Ext group
    """
    polynomial_system = []
    P = Poly(polynomial, *sorted(polynomial.free_symbols, key=lambda s: s.name))
    max_each_index = tuple(map(max, zip(*P.monoms())))[1:]
    N = max(max_each_index)
    print(polynomial.free_symbols)
    for i in range(len(polynomial.free_symbols)):
        for l in range(p**(N-max_each_index[i])):
            polynomial_system.append(create_polynomial(i,l))


def create_polynomial(i, l):
    pass


if __name__ == "__main__":
    P = generate_ppolynomial(POWER, NUM_OF_VAR) + symbols("x_0")
    print(f"starting with {P}")
    generate_iteration(P)
