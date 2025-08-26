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
