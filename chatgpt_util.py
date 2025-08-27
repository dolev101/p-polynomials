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
