from sympy import poly, symbols, Poly, collect, Expr
from sympy.polys.domains import GF
import itertools
from random import seed, sample
from typing import List
from chatgpt_util import eliminate_single_variable_polys, t_monomial_root_modp, decompose_over_basis, rename_vars_list_strict, num_vars, nth_roots_mod_prime
import math

# seed(4)
t = symbols('t')

p = 5
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
    """
    generate basic reduced ppolynomial, notice that there isn't linear part for smoothness 
    """
    A = get_random_A_subset(BASIS_Fpt, n, len)
    A = sorted(A, key=lambda p: Poly(p, t).degree())

    return sum([a*symbols(f"x_{i}")**(p**n) for i, a in enumerate(A)])

def generate_iteration(polynomial):
    """
    Take polynomial and return a system of equations by the algorithm for the Ext group
    """
    polynomial_system = []
    P = Poly(polynomial, *sorted(polynomial.free_symbols, key=lambda s: s.name))
    max_each_index = tuple(map(max, zip(*P.monoms())))[1:]
    N = int(math.log(int(max(max_each_index)), p))
    for i in range(len(max_each_index)): # number of variables
        Ni = int(math.log(int(max_each_index[i]), p))
        for l in range(p**(N-Ni)):
            polynomial_system.append(create_polynomial(i, N, Ni, l, polynomial))
    print(polynomial_system)
    eliminated = eliminate_single_variable_polys(polynomial_system, t)[0]
    removed_0 =[x for x in eliminated if x!=0]
    renamed, _ = rename_vars_list_strict(removed_0, exclude={t})
    return renamed

def get_c_ijk(polynomial, i, j, k, N, Ni):
    parts = collect(polynomial, symbols(f"x_{i}"), evaluate=False)
    cik = t*0 + parts.get(symbols(f"x_{i}") ** (p ** k), 0)
    # if i==2:
    #     print(k, cik, decompose_over_basis(cik, t, p**(k+N-Ni))[j], t_monomial_root_modp(decompose_over_basis(cik, t, p**(k+N-Ni))[j], t, p**(k+N-Ni), p)[0])
    return t_monomial_root_modp(decompose_over_basis(cik, t, p**(k+N-Ni))[j], t, p**(k+N-Ni), p)[0]

def create_polynomial(i, N, Ni, l, old_polynomial):
    polynomial = 0*t
    for k in range(Ni + 1):
        for j in range(p**(k+N-Ni)):
            c_ijk = get_c_ijk(old_polynomial, i, j, k,N,Ni)
            if c_ijk != 0:
                for r in range(p**N):

                    if (r+j+1-(l+1)*p**k) % (p**(k+N-Ni)) == 0:
                        # if i ==1:
                        #     print(i, k, j, r)
                        polynomial += c_ijk*(symbols(f"x_{r}")**(p**(Ni-k)))*t**((r+j+1-(l+1)*p**k)//(p**(k+N-Ni)))
    return polynomial

def dimension_formula(polynomial):
    P = Poly(polynomial, *sorted(polynomial.free_symbols, key=lambda s: s.name))
    max_each_index = tuple(map(max, zip(*P.monoms())))[1:]
    N = int(math.log(int(max(max_each_index)), p))
    # print(f"N of {polynomial} in {N}")
    final = p**N
    for i in range(len(max_each_index)): # number of variables
        Ni = int(math.log(int(max_each_index[i]), p))
        final -= p**(N-Ni)
    return final



def iterate_twice_check_for_non_stabilizing(P):
    first_dim = dimension_formula(P)
    first_iter = generate_iteration(P)
    
    if len(first_iter) != 1:
        print("problemo!!")
        print(first_iter)
    else:
        first_iter = first_iter[0]
        second_iter = generate_iteration(first_iter) # TODO: not always one equation
        if len(second_iter) == 1:
            second_iter = second_iter[0]
            print(f"{first_iter=}")
            print(f"{second_iter=}")
            third_dim = dimension_formula(second_iter)
            if first_dim != third_dim:
                print("completo")
                print(generate_iteration(second_iter))
                print(first_dim, third_dim)
        else:
            print(first_iter)
            print(second_iter)
            print("problemo!!!")

if __name__ == "__main__":
    P = generate_ppolynomial(POWER, NUM_OF_VAR) + symbols("x_0") #t**2* symbols("x_2")**p+t* symbols("x_1")#+ symbols("x_0")#**p + t* symbols("x_1")
    x_0, x_1, x_2 = symbols("x_0 x_1 x_2")
    # P = t**16*x_2**5 + t**6*x_1**5 + x_0**5 + x_0
    # P = t*x_1**5 + x_0**5 + x_0
    P = t**9*x_2**5 + t**4*x_1**5 + x_0**5 + x_0
    print(f"starting with {P}")
    # iterate_twice_check_for_non_stabilizing(P)
    print(generate_iteration(P))
