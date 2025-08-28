from sympy import sympify, Poly

from sympy import sympify, Poly, S

def eliminate_single_variable_polys(polys, t):
    """
    Given a list of polynomials (SymPy expr/Poly), repeatedly:
      1) find any polynomial that is a *single monomial* using only {x_i, t}
         with a positive power of some x_i (not t-only),
      2) remove that polynomial from the list,
      3) remove every monomial containing that x_i from all remaining polynomials
         (equivalently, substitute x_i -> 0).

    Parameters
    ----------
    polys : list[Expr|Poly]
    t     : Symbol   # the distinguished variable allowed to appear together with x_i

    Returns
    -------
    new_polys : list[Expr]
    eliminated_vars : list[Symbol]   # the x_i’s that were eliminated (sorted by name)
    """
    # normalize to expressions
    exprs = []
    for p in polys:
        exprs.append(p.as_expr() if hasattr(p, "as_expr") else sympify(p))

    eliminated = set()

    while True:
        to_eliminate_vars = set()
        to_drop_indices = set()

        for idx, e in enumerate(exprs):
            e = sympify(e).expand()
            if e == 0:
                continue

            vars_list = sorted(e.free_symbols, key=lambda s: s.name)
            if not vars_list:
                continue  # constant

            # Represent as Poly over *all* symbols so any symbol in coefficients is detected
            P = Poly(e, *vars_list, domain='EX')
            if P.is_zero:
                continue

            monoms = P.monoms()
            if len(monoms) != 1:
                continue  # not a single monomial

            exps = monoms[0]
            support_idx = [i for i, exp in enumerate(exps) if exp > 0]
            support_syms = {vars_list[i] for i in support_idx}

            # must be subset of {x_i, t} and must include exactly one non-t variable
            non_t_syms = {s for s in support_syms if s != t}
            if len(non_t_syms) == 1 and support_syms.issubset(non_t_syms | {t}):
                v = next(iter(non_t_syms))  # the x_i to eliminate
                to_eliminate_vars.add(v)
                to_drop_indices.add(idx)

        if not to_eliminate_vars:
            break

        # drop the eliminable monomial-polynomials
        exprs = [e for i, e in enumerate(exprs) if i not in to_drop_indices]

        # remove all monomials containing any eliminated x_i: substitute x_i -> 0
        subs_map = {v: S.Zero for v in to_eliminate_vars}
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
from sympy import Poly, sympify, S, Rational
from sympy.polys.polyerrors import PolynomialError

def t_monomial_root_modp(expr, t, n, p, allow_fractional=False):
    """
    If expr is a single monomial in t, c*t**k (with c independent of t),
    return all n-th roots over the form r * t**(k/n), where:
      - r runs over the n-th roots of c modulo prime p, computed by
        `nth_roots_mod_prime(c_mod, n, p)` (assumed available in scope),
      - the t-exponent is divided by n. If k % n != 0 and allow_fractional=False,
        raise ValueError; otherwise use Rational(k, n).

    Returns
    -------
    list[sympy.Expr]
        A list of expressions (possibly empty if no coefficient root exists).

    Notes
    -----
    - Coefficient c must be an integer modulo p. If not, ValueError is raised.
    - If expr == 0, returns [0].
    - If c == 0, the only n-th root mod p is 0; returns [0] (since 0 * t**anything == 0).
    """
    expr = sympify(expr)

    if n <= 0:
        raise ValueError("n must be a positive integer")
    if expr.is_zero:
        return [S.Zero]

    # Treat as polynomial in t; coefficients may involve other symbols (domain='EX')
    try:
        P = Poly(expr, t, domain='EX')
    except PolynomialError as e:
        raise ValueError("Expression is not a polynomial in t") from e

    monoms = P.monoms()
    coeffs = P.coeffs()

    if len(monoms) != 1:
        raise ValueError("Expression is not a single monomial in t")

    k = monoms[0][0]   # exponent of t
    c = coeffs[0]      # coefficient (must be independent of t)

    # Determine t-exponent root
    if k % n == 0:
        new_exp = k // n
    else:
        if not allow_fractional:
            raise ValueError(f"t-exponent {k} not divisible by n={n}")
        new_exp = Rational(k, n)

    # Coefficient must be an integer modulo p
    if not c.is_Integer:
        # Allow exact rationals that reduce mod p (denominator invertible mod p)
        if c.is_Rational:
            num = int(c.p)  # numerator
            den = int(c.q)  # denominator
            if den % p == 0:
                raise ValueError("Coefficient denominator not invertible modulo p")
            # Reduce num/den mod p
            c_mod = (num % p) * pow(den % p, -1, p) % p
        else:
            raise ValueError("Coefficient must be an integer (or rational) for mod-p n-th roots")
    else:
        c_mod = c % p

    # Handle c == 0 separately: the only n-th root mod p is 0
    if c_mod == 0:
        return [S.Zero]

    # Compute all n-th roots of the coefficient modulo p (assumes function exists)
    roots_mod_p = nth_roots_mod_prime(c_mod, n, p)  # -> iterable of ints in [0, p-1]

    # Build expressions
    results = [S(int(r)) * (t ** new_exp) for r in roots_mod_p]
    return results
