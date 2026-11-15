"""The scope of the vacuum term derivation, reduced to the lane's premises.

lambda_lane.py derives the traceless form within local actions,
algebraic in U, first order in the fields, with the gauge field
Pontryagin case closed by theta bar = 0.  This probe removes those
conditions one by one.

S1 (the axiom on the effective action). The empty field's energy,
   classical or induced by the fluctuations, weighs no history.  The
   induced vacuum term enters with the same density as the classical
   constant (composite_pole.py K1). So the classification below applies
   to the effective action, at every loop order.

S2 (every local order at once). A density whose integral is the same
   for every history with fixed boundary data has identically vanishing
   Euler-Lagrange expressions: it is a null Lagrangian.  If the space of
   field values is star shaped about the zero field, every null
   Lagrangian of finite order is
       mu[phi] = Div K(mu)[phi] + mu[0](x).
   This is a total divergence of a current built from the fields, plus
   a functionally constant term equal to its value at the zero field (arXiv:
   2406.17048 Corollary D.1, Eq. D.15; Olver, Applications of Lie
   Groups to Differential Equations, Thm 4.7).  Globally the leftover
   classes are the de Rham cohomology of the field space (Anderson, The
   Variational Bicomplex). The carrier's C^3 over R^4 is contractible,
   so there are none.  Checked: the carrier's Phi = dJ ^ dJ is
   Div(eps J dJ) with Phi[0] = 0, by finite differences on a smooth
   configuration.

S3 (translation invariance removes the divergence). In any
   translation invariant state <Div K> = Div <K> = 0 for every local
   current built from the carrier, at every order and beyond
   perturbation theory.  So the vacuum value of the density is mu[0],
   its value on the empty field. That is a field independent density,
   the one in the traceless form.  A nonzero divergence needs a field
   profile (a ruler). The empty field vacuum does not carry one
   (gate C), and a ruler would give the graviton a mass (Carone
   arXiv:1911.13166 Sec. III).  On a periodic lattice the average of a
   discrete divergence vanishes by telescoping. That check is an
   identity. Its content is the premise that the vacuum is translation
   invariant.  A ruler profile gives a nonzero average.

S4 (the gauge field case). The Pontryagin density P = tr F F~ is the
   divergence of the Chern-Simons current, which is not gauge
   invariant. The space of connections modulo gauge is not star shaped,
   theta vacua give <P> != 0, and S3 does not apply.  The volume element
   is the determinant of the metric, and the metric is the common mode
   of the one field. A gauge field is a named row, and in this lane
   gravity couples to no named field.  That premise excludes the case.
   theta bar = 0 (interference/couplings.py) excludes it independently.
   Stated.

S5 (nonlocal densities). The nonlocal densities that make the
   constant's term history independent are global normalizations, such
   as sqrt g / int sqrt g, the shape of vacuum energy sequestering
   (arXiv:1309.6562).  They cancel the classical constant but not the
   induced vacuum term, unless the matter scale varies over histories.
   The anchor fixes m_e / M_b locally (candidates L-sequester, FAIL).
   Every other nonlocal density is excluded by the axiom itself:
   "the locality is an axiom of self-presence: a nonlocal self-term
   would be interference with a presence elsewhere, i.e. with a second
   thing, which is coupling, not self-interference"
   (interference/__init__.py).  In a one-field world every term is a
   self-term, so none is nonlocal.  Bell violations do not count
   against this. The framework is local and nonrealist: it rejects the
   absolute outcome and keeps local dynamics (tests/probes/bell.py).
   The conjugate of Lambda, the total 4-volume, is a global property of
   the state. It is a record level quantity like a total charge, not a
   nonlocal term of the action.  Stated.

S6 (an independent 3-form). Henneaux and Teitelboim (Phys. Lett. B
   222 (1989) 195) replace the fixed density by dC, C a 3-form gauge
   field with no local degrees of freedom and one global mode whose
   conjugate is Lambda. DC is gauge invariant, so S4's Chern-Simons
   mechanism does not reach it, and it restores full diffeomorphisms.
   But it is a second field coupled to the volume. That is interference
   of the carrier with something other than itself, which self-presence
   excludes ("a second thing, which is coupling, not
   self-interference", interference/__init__.py).  In a one-field
   world the null Lagrangians of S2 are built from the carrier alone,
   and the fixed density is what remains.  Classically the two have
   the same content. The one-field premise itself decides between
   them.  Stated.

RESULT. The scope of the derivation reduces to premises the lane
already holds and to the axiom it starts from: the empty field vacuum
(translation invariance), the anchor, gravity coupling to the unnamed
common mode, and self-presence, which supplies both the no weight
condition and locality.  The fixed density of the traceless form is
the value of the density on the empty field.

Usage: python3 tests/probes/lambda_scope.py     (a few seconds)
"""

import itertools

import numpy as np

EPS4 = np.zeros((4, 4, 4, 4))
for _p in itertools.permutations(range(4)):
    EPS4[_p] = np.linalg.det(np.eye(4)[list(_p)])


def _field(rng, n_waves=3):
    """A smooth three component complex configuration: sums of plane waves."""
    k = rng.normal(size=(3, n_waves, 4))
    a = rng.normal(size=(3, n_waves)) + 1j * rng.normal(size=(3, n_waves))

    def phi(x):
        return np.array([np.sum(a[c] * np.exp(1j * (k[c] @ x))) for c in range(3)])
    return phi


def _grad(f, x, h=1e-3):
    out = []
    for mu in range(4):
        e = np.zeros(4); e[mu] = h
        out.append((-f(x + 2 * e) + 8 * f(x + e) - 8 * f(x - e) + f(x - 2 * e)) / (12 * h))
    return np.array(out)                      # shape (4, ...)


def s2(rng, n_points=6):
    phi = _field(rng)
    J = lambda x: np.imag(np.conj(phi(x)) @ _grad(phi, x).T)          # J_nu = Im sum conj(phi) d_nu phi
    dJ = lambda x: _grad(J, x)                                         # dJ[mu, nu] = d_mu J_nu
    K = lambda x: np.einsum("mnrs,n,rs->m", EPS4, J(x), dJ(x))       # K^mu = eps J_nu d_rho J_sigma
    worst = 0.0
    for _ in range(n_points):
        x = rng.normal(size=4)
        div = np.trace(_grad(K, x))                                    # d_mu K^mu
        F = dJ(x)
        phi_density = np.einsum("mnrs,mn,rs->", EPS4, F, F)            # eps d_mu J_nu d_rho J_sigma
        worst = max(worst, abs(div - phi_density) / max(1.0, abs(phi_density)))
    zero = np.einsum("mnrs,mn,rs->", EPS4, np.zeros((4, 4)), np.zeros((4, 4)))
    return dict(div_vs_phi=float(worst), phi_at_zero=float(zero))


def s3(rng, L=6):
    """Lattice average of a nonlinear discrete divergence: periodic ensemble vs a ruler."""
    g = rng.normal(size=(3, L, L, L, L)) + 1j * rng.normal(size=(3, L, L, L, L))
    phi = g + 0.3 * g * g                                             # non-Gaussian, translation invariant
    fwd = lambda f, mu: np.roll(f, -1, axis=mu + 1) - f
    J = np.stack([np.imag(np.sum(np.conj(phi) * fwd(phi, mu), axis=0)) for mu in range(4)])
    K = np.stack([np.abs(J[mu])**3 + J[(mu + 1) % 4] * J[mu] for mu in range(4)])   # any local current
    div = sum(np.roll(K[mu], -1, axis=mu) - K[mu] for mu in range(4))
    periodic = float(abs(div.mean()) / np.abs(div).mean())
    xs = np.meshgrid(*[np.arange(L, dtype=float)] * 4, indexing="ij")
    ruler = np.stack([xs[0] + 1j * xs[1], xs[2] + 1j * xs[3], np.zeros_like(xs[0])])
    dr = np.stack([[np.gradient(ruler[c], axis=mu) for mu in range(4)] for c in range(3)])
    Fr = 2 * np.imag(np.einsum("cm...,cn...->mn...", np.conj(dr), dr))
    ruler_phi = float(np.einsum("mnrs,mn...,rs...->...", EPS4, Fr, Fr).mean())
    return dict(periodic_average=periodic, ruler_average=ruler_phi)


def run(report=print, seed=3):
    rng = np.random.default_rng(seed)
    a, b = s2(rng), s3(rng)
    assert a["div_vs_phi"] < 1e-6 and a["phi_at_zero"] == 0.0
    assert b["periodic_average"] < 1e-12 and abs(b["ruler_average"]) > 1.0
    report("S1 the axiom applies to the effective action: the induced vacuum term enters with the "
           "same density as the classical constant.")
    report(f"S2 null Lagrangian = Div K + mu[0]: Phi = d_mu(eps J dJ)^mu to {a['div_vs_phi']:.1e}, "
           f"Phi[0] = {a['phi_at_zero']:.0f}; C^3 over R^4 is contractible, no topological classes.")
    report(f"S3 periodic, non-Gaussian, translation invariant ensemble: |<Div K>| / <|Div K|> = "
           f"{b['periodic_average']:.1e} (telescoping, an identity); a ruler profile gives "
           f"<Phi> = {b['ruler_average']:.1f}.  The vacuum value of the density is mu[0], the "
           f"empty field's.")
    report("S4 the gauge Pontryagin case: excluded because gravity couples to no named field. "
           "Theta bar = 0 "
           "excludes it independently.  S5 nonlocal global normalizations reduce to "
           "sequestering, which the anchor excludes. Every other nonlocal term is excluded by "
           "self-presence.")
    report("S6 an independent 3-form (Henneaux-Teitelboim) is a second field coupled to the volume. "
           "Self-presence excludes it, and the fixed density is what the carrier alone leaves.")
    report("RESULT: the derivation rests on self-presence (the no weight condition and locality) "
           "and premises the lane already holds: the empty field vacuum, the anchor, the unnamed "
           "common mode.")
    return dict(s2=a, s3=b)


if __name__ == "__main__":
    run()
