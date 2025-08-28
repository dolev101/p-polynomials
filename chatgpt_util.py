from sympy import sympify, Poly

def eliminate_single_variable_polys(polys):
    """
    Given a list of polynomials (SymPy expr or Poly), repeatedly:
      1) find any polynomial that is a single-variable monomial c*X_j^m
         (no other variables multiply it),
      2) remove it from the list,
      3) set that variable to 0 in all remaining polynomials (thus removing
         all monomials containing that variable).

    Returns
    -------
    new_polys : list[sympy.Expr]
        The updated list after eliminations.
    eliminated_vars : list[sympy.Symbol]
        Variables that were eliminated (sorted by name).
    """

    # normalize to expressions
    exprs = []
    for p in polys:
        if hasattr(p, "as_expr"):
            exprs.append(p.as_expr())
        else:
            exprs.append(sympify(p))

    eliminated = set()

    while True:
        to_eliminate_vars = set()
        to_drop_indices = set()

        for idx, e in enumerate(exprs):
            e = e.expand()
            if e == 0:
                continue  # nothing to do

            vars_list = sorted(e.free_symbols, key=lambda s: s.name)
            if not vars_list:
                continue  # constant, ignore

            # multivariate structure: check if exactly one monomial overall,
            # and that monomial involves exactly one variable with positive exponent
            Pf = Poly(e, *vars_list)
            monoms = Pf.monoms()
            if len(monoms) != 1:
                continue

            mono = monoms[0]  # tuple of exponents aligned with vars_list
            # count how many variables have positive exponent
            support = [i for i, exp in enumerate(mono) if exp > 0]
            if len(support) == 1:
                # it's c * X_j^m, with only one variable present
                j = support[0]
                v = vars_list[j]
                to_eliminate_vars.add(v)
                to_drop_indices.add(idx)

        if not to_eliminate_vars:
            break

        # remove those single-variable polynomials
        exprs = [e for i, e in enumerate(exprs) if i not in to_drop_indices]

        # remove all monomials containing any eliminated variable (subs var->0)
        subs_map = {v: 0 for v in to_eliminate_vars}
        exprs = [e.subs(subs_map).expand() for e in exprs]

        eliminated |= to_eliminate_vars

    return exprs, sorted(eliminated, key=lambda s: s.name)

from sympy import Symbol, sympify
from sympy import preorder_traversal

def _symbols_by_appearance(exprs):
    if not isinstance(exprs, (list, tuple)):
        exprs = [exprs]
    seen = []
    for e in exprs:
        e = sympify(e)
        for node in preorder_traversal(e):
            if getattr(node, "is_Symbol", False) and node not in seen:
                seen.append(node)
    return seen

def rename_vars_strict(expr, prefix="x", start=0, gens=None, exclude=None, exclude_names=None):
    """
    Rename all symbols in `expr` to prefix_i (i=start, ...), skipping any in `exclude`
    (set of Symbols) or `exclude_names` (set of strings).
    """
    expr = sympify(expr)
    if gens is None:
        gens = _symbols_by_appearance(expr)
    exclude = set() if exclude is None else set(exclude)
    exclude_names = set() if exclude_names is None else set(exclude_names)

    rename_syms = [s for s in gens if s not in exclude and s.name not in exclude_names]
    mapping = {s: Symbol(f"{prefix}_{start+i}") for i, s in enumerate(rename_syms)}
    return expr.xreplace(mapping), mapping

def rename_vars_list_strict(exprs, prefix="x", start=0, gens=None, exclude=None, exclude_names=None):
    """
    Same as above, but apply one consistent mapping across a list of exprs.
    """
    exprs = [sympify(e) for e in exprs]
    if gens is None:
        gens = _symbols_by_appearance(exprs)
    exclude = set() if exclude is None else set(exclude)
    exclude_names = set() if exclude_names is None else set(exclude_names)

    rename_syms = [s for s in gens if s not in exclude and s.name not in exclude_names]
    mapping = {s: Symbol(f"{prefix}_{start+i}") for i, s in enumerate(rename_syms)}
    return [e.xreplace(mapping) for e in exprs], mapping
from sympy import sympify

def num_vars(poly, exclude=None):
    """
    Count distinct SymPy symbols in `poly`.
    Optionally exclude some symbols by name or Symbol.
    """
    expr = sympify(poly)
    syms = set(expr.free_symbols)
    if exclude:
        excl = {getattr(e, "name", e) for e in exclude}  # names or Symbols
        syms = {s for s in syms if (s not in exclude and s.name not in excl)}
    return len(syms)

from math import gcd
from sympy import mod_inverse
from sympy.ntheory.residue_ntheory import primitive_root, discrete_log

def nth_roots_mod_prime(a: int, n: int, p: int):
    """
    Solve x^n ≡ a (mod p) for prime p.
    Returns a sorted list of all solutions in [0, p-1].
    """
    if p < 2:
        raise ValueError("p must be a prime ≥ 2")
    if n <= 0:
        raise ValueError("n must be ≥ 1")
    a %= p
    if a == 0:
        return [0]

    order = p - 1
    g = gcd(n, order)
    # necessary & sufficient condition for existence (a in subgroup of size (p-1)/g)
    if pow(a, order // g, p) != 1:
        return []

    # easy case: n invertible mod (p-1)
    if g == 1:
        inv = mod_inverse(n, order)
        return [pow(a, inv, p)]

    # general case: use discrete log in the cyclic group F_p^*
    gen = primitive_root(p)                  # generator of F_p^*
    k = discrete_log(p, a, gen)              # a ≡ gen^k (mod p)
    if k % g != 0:
        return []                            # (redundant check; should already hold)

    m = order // g
    n_prime = n // g
    inv_nprime = mod_inverse(n_prime, m)     # since gcd(n', m)=1
    r0 = ((k // g) * inv_nprime) % m         # one solution exponent
    x0 = pow(gen, r0, p)

    # all g solutions differ by g-th roots of unity
    zeta = pow(gen, m, p)                    # order g
    sols = [(x0 * pow(zeta, j, p)) % p for j in range(g)]
    return sorted(set(sols))

from sympy import symbols, Symbol, Poly, Matrix, fraction, together, S, sympify

def _group_by_mod_p(poly_in_t, t, p, s):
    # Poly over EX so coefficients may depend on other symbols
    P = Poly(poly_in_t, t, domain='EX')
    coeffs = [S.Zero]*p
    for k, c in P.as_dict().items():
        k = k[0]
        i = k % p
        m = k // p
        coeffs[i] += c * (s**m)
    return coeffs

def _mul_mod_tp(a, b, p, s):
    out = [S.Zero]*p
    for i in range(p):
        ai = a[i]
        if ai == 0:
            continue
        for j in range(p):
            bj = b[j]
            if bj == 0:
                continue
            r = (i + j) % p
            m = (i + j) // p
            out[r] += ai * bj * (s**m)
    return out

def decompose_over_basis(expr, t, p, return_coeffs_in='t'):
    """
    Return [g_0, ..., g_{p-1}] with expr = sum g_i(t^p) * t**i.
    """
    s = Symbol('_s')  # stands for t**p inside coefficients

    expr = together(expr)
    num, den = fraction(expr)

    a = _group_by_mod_p(num, t, p, s)   # numerator grouped by i mod p
    q = _group_by_mod_p(den, t, p, s)   # denominator grouped by i mod p

    # Build linear system for u with (∑ q_i t^i)(∑ u_i t^i) = 1 modulo t^p = s
    M = [[S.Zero]*p for _ in range(p)]
    for r in range(p):
        for j in range(p):
            acc = S.Zero
            for i in range(p):
                if (i + j) % p == r:
                    acc += q[i]*(s**((i + j)//p))
            M[r][j] = acc
    M = Matrix(M)
    e0 = Matrix([1] + [0]*(p-1))

    u = list(M.LUsolve(e0))            # u_i(s)
    g = _mul_mod_tp(a, u, p, s)        # g_i(s)

    if return_coeffs_in == 't':
        # ensure SymPy objects, then substitute s -> t**p
        return [sympify(gi).xreplace({s: t**p}) for gi in g]
    elif return_coeffs_in == 's':
        return g
    else:
        raise ValueError("return_coeffs_in must be 't' or 's'")
