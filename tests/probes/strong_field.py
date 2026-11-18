"""Static strong fields of the traceless equations: the solution, its horizon, and where it holds.

THE QUESTION.  Outside a static spherical source the induced action
gives gamma = beta = 1, "the expansion of the Schwarzschild metric" to
second order (post_newtonian.py N2).  What are the exact static
spherical vacuum solutions of the lane's traceless equations, what
does their horizon carry, and down to what mass does the
two derivative action that gives them hold?

PREMISES.  P1 the vacuum traceless equations R_mn - (1/4) g_mn R = 0.
P2 the metric -f(r) dt^2 + dr^2/h(r) + r^2 dOmega^2.  P3 the entropy of
a horizon of area A on the one list, "A/(4G) = A Sigma/(8 pi t0)", with
the G of the anchor (horizon_entropy.py).  P4 the kernel's departures
from the Einstein value, 10 and 50 percent at q = 0.129 and 0.430 M_b
(horizon_amount.py A4), and G = M_b^-2 (root.py).

STATEMENTS (the numbers below check each one).

  F1 (h = f).  The combination h E_rr + E_tt/f of P1 is proportional
     to f h' - h f': the equations force h = k f, and a rescaling of t
     sets k = 1.
  F2 (one equation, two constants).  With h = f every component of P1
     is proportional to r^2 f'' - 2 f + 2, whose general solution is
     f = 1 + A/r + B r^2.  The curvature scalar on it is -12 B, so
     B = -Lambda/3 with R = 4 Lambda, the constant of
     expansion_form.py E2, and the Newtonian limit sets A = -2 G M:
         f = 1 - 2 G M/r - Lambda r^2/3.
     The vacuum equations fix neither constant.
  F3 (the horizon's temperature).  Near a zero r_h of f the Euclidean
     metric is rho^2 (f'(r_h)/2)^2 dtau^2 + d rho^2, regular for the
     period 4 pi/f'(r_h): T = f'(r_h)/(4 pi) = (1 - Lambda r_h^2)/
     (4 pi r_h), and 1/(8 pi G M) at Lambda = 0.
  F4 (the first law, with one G).  With P3's entropy pi r_h^2/G and
     M(r_h) = r_h (1 - Lambda r_h^2/3)/(2 G), dM = T dS holds
     identically in r_h at fixed Lambda: the G of the horizon's
     entropy is the G of the metric's mass term.
  F5 (where it holds).  At Lambda = 0 the Kretschmann scalar at the
     horizon is 3/(4 G^4 M^4). Its fourth root, 0.931 M_b^2/M, reaches
     the kernel's 10 and 50 percent departures at M = 7.2 and 2.2 M_b.
     Above about seven boundary masses the horizon lies where the
     action matches Einstein's to within a tenth; a horizon of one solar mass sits
     at 1.0e-38 M_b, 8e-38 of the 10 percent scale.
  F6 (what fails).  A linear term gamma r added to f, the potential of
     conformal gravity (Mannheim and Kazanas, Astrophys. J. 342, 635),
     leaves a residual in P1.  So does the screened Poisson potential
     f = 1 - 2 G M exp(-r/xi0)/r with the healing length
     xi0 = l_Pl/d10 (sakharov_normalisation.py, carrier_sakharov.py).

RESULT.  The static spherical vacuum solutions of the traceless
equations are the Schwarzschild de Sitter family, with the mass and
Lambda its two integration constants.  The horizon has the temperature
its Euclidean regularity gives, and with the one list's area term the
first law holds with one G.  The two derivative action that yields the
solution holds at the horizon for masses above about seven boundary
masses.

Scope: static, spherical, vacuum, the two derivative action.  Rotation,
collapse, evaporation, the phases of a horizon and the region below a
few boundary masses lie outside this probe.

Usage: python3 tests/probes/strong_field.py     (about 10 s)
"""

import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import root as r_  # noqa: E402

DEPARTURES = (0.129, 0.430)                # horizon_amount.py A4, in M_b (its test freezes them)
SOLAR_IN_PLANCK = 1.98841e30 / 2.176434e-8  # one solar mass in Planck masses (IAU nominal over CODATA)

t, r, th, ph = sp.symbols("t r theta phi", positive=True)
X = (t, r, th, ph)


def curvature(f, h):
    """Ricci tensor, curvature scalar and Kretschmann scalar of P2."""
    g = sp.diag(-f, 1 / h, r**2, r**2 * sp.sin(th)**2)
    gi = g.inv()
    Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                 for d in range(4)) / 2 for c in range(4)] for b in range(4)] for a in range(4)]

    def riem(a, b, c, d):                                  # R^a_{bcd}
        out = sp.diff(Gam[a][b][d], X[c]) - sp.diff(Gam[a][b][c], X[d])
        return out + sum(Gam[a][c][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][c] for e in range(4))

    Rm = [[[[sp.simplify(riem(a, b, c, d)) for d in range(4)] for c in range(4)] for b in range(4)]
          for a in range(4)]
    Ric = sp.Matrix(4, 4, lambda b, d: sp.simplify(sum(Rm[a][b][a][d] for a in range(4))))
    R = sp.simplify(sum(gi[b, d] * Ric[b, d] for b in range(4) for d in range(4)))
    low = [[[[g[a, a] * Rm[a][b][c][d] for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
    K = sp.simplify(sum(low[a][b][c][d]**2 * gi[a, a] * gi[b, b] * gi[c, c] * gi[d, d]
                        for a in range(4) for b in range(4) for c in range(4) for d in range(4)))
    return g, Ric, R, K


def trace_free(f, h):
    """The diagonal components of P1 (the others vanish for P2)."""
    g, Ric, R, _ = curvature(f, h)
    E = sp.Matrix(4, 4, lambda b, d: sp.simplify(Ric[b, d] - g[b, d] * R / 4))
    assert all(E[b, d] == 0 for b in range(4) for d in range(4) if b != d)
    return [E[i, i] for i in range(4)]


def residual(f, h=None):
    """Whether a given f (and h, f by default) leaves a residual in P1."""
    return any(sp.simplify(e) != 0 for e in trace_free(f, f if h is None else h))


def run(report=print):
    out = {}
    F, H = sp.Function("f", positive=True)(r), sp.Function("h", positive=True)(r)

    # F1: h = f.
    E = trace_free(F, H)
    comb = sp.simplify(H * E[1] + E[0] / F)
    wronskian = F * sp.diff(H, r) - H * sp.diff(F, r)
    factor = sp.simplify(comb / wronskian)
    assert factor != 0 and not factor.has(sp.Derivative)
    k = sp.Symbol("k", positive=True)
    assert sp.simplify(comb.subs(H, k * F).doit()) == 0
    out.update(wronskian_factor=str(factor))
    report(f"F1 h = f: h E_rr + E_tt/f = ({factor}) (f h' - h f'), which vanishes for h = k f alone")

    # F2: one equation, two constants.
    E1 = trace_free(F, F)
    ode = r**2 * sp.diff(F, r, 2) - 2 * F + 2
    ratios = [sp.simplify(e / ode) for e in E1]
    assert all(q != 0 and not q.has(sp.Derivative) for q in ratios)
    A, B = sp.symbols("A B", real=True)
    general = 1 + A / r + B * r**2
    assert sp.simplify(ode.subs(F, general).doit()) == 0
    sol = sp.dsolve(sp.Eq(ode, 0), F).rhs
    consts = sorted(sol.free_symbols - {r}, key=str)
    assert len(consts) == 2 and sp.simplify(ode.subs(F, sol).doit()) == 0
    assert {sp.simplify(sp.diff(sol, c)) for c in consts} == {1 / r, r**2}
    assert sp.simplify(sol.subs({c: 0 for c in consts})) == 1
    _, _, R_general, _ = curvature(general, general)
    assert sp.simplify(R_general + 12 * B) == 0
    G, M, Lam = sp.symbols("G M Lambda", positive=True)
    sds = 1 - 2 * G * M / r - Lam * r**2 / 3
    _, Ric_s, R_s, K_s = curvature(sds, sds)
    assert not residual(sds) and sp.simplify(R_s - 4 * Lam) == 0
    out.update(constants=2, scalar_curvature="4 Lambda")
    report("F2 one equation: every component is proportional to r^2 f'' - 2 f + 2, solved by "
           "f = 1 + A/r + B r^2; R = -12 B = 4 Lambda and A = -2 G M")

    # F3: the horizon's temperature.
    rh, rho = sp.symbols("r_h rho", positive=True)
    f_prime = sp.diff(sds, r)
    M_of = sp.solve(sds.subs(r, rh), M)[0]
    on = lambda expr: expr.subs(M, M_of)
    fp = sp.simplify(on(f_prime).subs(r, rh))
    r_of = rh + fp * rho**2 / 4                             # proper distance rho from the horizon
    f_near = sp.series(on(sds).subs(r, r_of), rho, 0, 4).removeO()
    radial = sp.series(sp.diff(r_of, rho)**2 / on(sds).subs(r, r_of), rho, 0, 2).removeO()
    assert sp.simplify(f_near - (fp * rho / 2)**2) == 0 and sp.simplify(radial - 1) == 0
    T = sp.simplify((f_prime.subs(r, rh) / (4 * sp.pi)).subs(M, M_of))
    assert sp.simplify(T - (1 - Lam * rh**2) / (4 * sp.pi * rh)) == 0
    T0 = sp.simplify(T.subs(Lam, 0).subs(rh, 2 * G * M))
    assert sp.simplify(T0 - 1 / (8 * sp.pi * G * M)) == 0
    out.update(temperature=str(T0))
    report(f"F3 temperature: regular for the period 4 pi/f'(r_h); T = (1 - Lambda r_h^2)/(4 pi r_h), "
           f"and {T0} at Lambda = 0")

    # F4: the first law with one G.
    S = sp.pi * rh**2 / G                                   # A/(4G), A = 4 pi r_h^2
    first_law = sp.simplify(sp.diff(M_of, rh) - T * sp.diff(S, rh))
    Gs = sp.Symbol("G_S", positive=True)                    # an entropy read with another G
    broken = sp.simplify(sp.diff(M_of, rh) - T * sp.diff(S.subs(G, Gs), rh))
    assert first_law == 0 and broken != 0
    out.update(first_law=True)
    report("F4 first law: dM = T dS identically in r_h at fixed Lambda with S = A/(4G). An entropy read "
           "with another G leaves " + str(sp.factor(broken)))

    # F5: where the two derivative action holds.
    K_h = sp.simplify(K_s.subs(Lam, 0).subs(r, 2 * G * M))
    assert sp.simplify(K_h - sp.Rational(3, 4) / (G * M)**4) == 0
    scale = 0.75**0.25                                      # K^(1/4) M / M_b^2 at the horizon, G = M_b^-2
    floors = tuple(scale / q for q in DEPARTURES)
    solar = scale / SOLAR_IN_PLANCK
    assert abs(scale - 0.9306) < 1e-4 and abs(floors[0] - 7.21) < 0.01 and abs(floors[1] - 2.16) < 0.01
    assert 1.0e-38 < solar < 1.1e-38 and solar / DEPARTURES[0] < 1e-37
    out.update(horizon_scale=scale, mass_floors=floors, solar_scale=solar, solar_ratio=solar / DEPARTURES[0])
    report(f"F5 regime: K^(1/4) at the horizon is {scale:.3f} M_b^2/M. It reaches the kernel's 10 and 50 percent "
           f"departures at M = {floors[0]:.1f} and {floors[1]:.1f} M_b. One solar mass sits at {solar:.1e} M_b, "
           f"{solar / DEPARTURES[0]:.0e} of the 10 percent scale")

    # F6: what fails.
    gam = sp.Symbol("gamma", positive=True)
    xi0 = sp.sqrt(G) / r_.d10                              # the healing length l_Pl/d10, hbar = c = 1
    fails = {"linear term gamma r": residual(1 - 2 * G * M / r + gam * r),
             "screened Poisson, xi0 = l_Pl/d10": residual(1 - 2 * G * M * sp.exp(-r / xi0) / r)}
    assert all(fails.values()) and not residual(1 - 2 * G * M / r)
    out.update(fails=fails)
    report("F6 what fails: " + ", ".join(fails) + " each leave a residual; f = 1 - 2 G M/r leaves none")

    report("RESULT: the static spherical vacuum solutions are f = h = 1 - 2 G M/r - Lambda r^2/3, the mass and "
           "Lambda the two integration constants; T = 1/(8 pi G M) and dM = T dS with one G. The action holds "
           f"at the horizon above about {floors[0]:.0f} boundary masses.")
    return out


if __name__ == "__main__":
    run()
