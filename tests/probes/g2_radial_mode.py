"""G2 radial mode: the colour singlet the G2 face condensate leaves.

THE QUESTION. THE GLUON CONDITION (couplings.py) and scale_ties.py S1
read the 125 GeV boson as the radial mode of the electroweak
condensate. G2 is broken by a condensate of its own, whose vacuum value
the connector fixes at v_break = v_EW/sqrt 2. That condensate has a
radial mode (the fluctuation of its amplitude): a colour singlet with a
gluon coupling. What is its mass?

INPUTS. v_EW and m_t are read from web depth 4 nodes (root.WEB). v_EW is
M_Pl exp(-S_quark), with the Singh response in S_quark, and closes G_F
at 0.02 ppm. m_t rides on m_tau. alpha_G2(mu*) = pi/32 and E[v^2] = 1/2
have one reading each. Loop order below is beta order, which is
transport, not web depth (interference/__init__.py, DEPTH VOCABULARY).
Cited loop mathematics: the one loop vector determinant (Coleman and
Weinberg, PRD 7, 1888 (1973), and Gildener and Weinberg, PRD 13, 3333
(1976)) and the low energy theorem at leading order.

PREMISES (each one is committed elsewhere in the engine).
  Q1  The G2 scalars are not a row at the regulator (boundary.py). The
      7 carries no action of its own (couplings.py, sector
      completeness). Its vacuum value is the connector's per mode
      reading of the electroweak condensate, not a binding of its own.
      Below the electroweak scale it is one real 7, phi, with G2 gauge
      couplings. With no sector it has no Yukawa term.
  Q2  THE GLUON CONDITION: the connector transports the vacuum value,
      not the field.
  Q3  One covariant operator gives a row both its share of b0 and its
      stiffness (boundary.py), so a row cannot run without also
      venting. The coset vectors are rows of that list. Their operator
      is read in the engine as the alpha_s threshold and, by the low
      energy theorem, as the gluon coupling. Its field dependent vacuum
      term sigma^4 ln sigma^2 is a coefficient of the same operator,
      read here as the curvature of the potential of the condensate's
      amplitude sigma.
  Q4  The coset modes are read as nonpropagating below matching
      (couplings.py, WHY THE COSET MODES DON'T APPEAR IN THE
      SPECTRUM): the six vectors and the six Goldstones, which move
      along the coset G2/SU(3) = S^6. The radial mode is the
      fluctuation of |phi|. It is constant on the coset, a colour
      singlet that carries no coset quantum number, so that reading
      does not reach it.

THE MENU (what else the 7's potential can hold).
  F1  A portal term phi^2 H^dagger H. The vacuum value then follows the
      field h, so M_V follows h: reading B of scale_ties.py, with
      kappa_g = -17/4. It fails THE GLUON CONDITION.
  F2  A mass term with a quartic. Two parameters, of which the
      connector fixes one combination. The quartic, and so the mass,
      stay free. It does not close.
  F3  No scale of its own. The quartic transmutes into the vacuum value
      (dimensional transmutation), and the connector fixes that value.
      One parameter, fixed. It closes.
  F3 is the closing form. F1 fails a measured gate, and F2 leaves a
  coupling unfixed.

STATEMENTS (the numbers below check each one).
  R1  G2 built from the octonion structure constants. A unit vector in
      the 7 is fixed by su(3) (8 massless vectors). The six coset
      vectors are degenerate, M_V^2 = g^2 v^2 C2(7)/6. All six
      Goldstones are eaten, so exactly one radial mode remains.
  R2  F3 gives m^2 = 3 sum_a M_a^4 / (8 pi^2 v^2) = g^4 v^2 / (4 pi^2),
      that is m_sigma = 2 alpha_G2 v_break. The code evaluates the
      potential directly, with the quartic set by the minimum, at three
      renormalisation scales, and recovers the closed form.
  R3  alpha_G2 read at mu* (where M_V is read), at v_break, and at M_V
      (the determinant's own scale) gives 34.2, 36.5 and 39.6 GeV.
  R4  At leading order the coset threshold moves 1/alpha_s by 21/(12 pi)
      per unit ln v and the top by -1/(3 pi): a ratio of -21/4. Per
      unit field (v_EW/v_break = sqrt 2) the mode's gluon amplitude is
      -7.42 times the top loop amplitude of a scalar of the same mass,
      a rate 55 times larger. This is the order of the threshold, the
      order scale_ties.py reads (its two loop chain cross check gives
      -4.52). The electroweak vectors take no share of b0 (boundary.py),
      so they are G2 singlets and electric charge commutes with G2. The
      coset vectors are neutral, and the mode has no photon coupling at
      one loop. At tree level it couples to the coset vectors alone.
      Mixing with h, through the quarks in the 7, starts at two loops.
  R5  F3 is not a general rule. Applied to h it gives m_h^2 < 0,
      because the top's term exceeds the W and Z terms
      ((m_t/M_Z)^4 > 3/4 suffices, since M_W < M_Z). The electroweak
      scale is the F4 knot's action S_quark, so h's potential belongs
      to that binding, not to a determinant. m_H stays data
      (higgs_ew.py), and the Coleman and Weinberg cubic coupling
      kappa_lambda = 5/3 is not available to h.

RESULT. The G2 face radial mode has m_sigma = 2 alpha_G2 v_break =
34.2 GeV with alpha_G2 at mu*, and 34 to 40 GeV over the determinant's
scale. It is the closing form on the menu F1 to F3. It is a colour
singlet whose only one loop coupling to Standard Model fields is to
gluons, through the same coset determinant that sets the alpha_s
threshold. No reference datum enters the mass.

EXPOSURE. A gluon fusion resonance decaying to two gluons near 34 to
40 GeV, at about 55 times the top loop rate: roughly 20 nb at 13 TeV
with the point like gg -> S cross section of the LHC Higgs cross section
working group (395 pb at 35 GeV). Searches that tag the event with an
initial state radiation (ISR) photon do not reach it. Gluon initial
states carry no charge, so the photon must leave a quark line
(ATLAS arXiv:2606.03222). Searches tagged by an ISR jet see gluon
initial states but start at 50 GeV (CMS arXiv:1710.00159). A jet ISR
dijet search below 50 GeV is the test.

Usage: python3 tests/probes/g2_radial_mode.py     (under a second)
"""

import contextlib
import io
import itertools
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import masses, root  # noqa: E402
import scale_ties  # noqa: E402

# Fano plane triples: e_i e_j = e_k for each cyclic (i, j, k).
FANO = ((1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5))
M_Z = 91.1876                         # evaluation coordinate (higgs_ew.py)


def g2_generators():
    """Orthonormal real basis of g2 in so(7): -Tr(X_a X_b) = delta_ab."""
    phi = np.zeros((7, 7, 7))
    for t in FANO:
        for p in itertools.permutations(range(3)):
            i, j, k = (t[p[0]] - 1, t[p[1]] - 1, t[p[2]] - 1)
            phi[i, j, k] = np.linalg.det(np.eye(3)[list(p)])
    so7 = []
    for a, b in itertools.combinations(range(7), 2):
        X = np.zeros((7, 7))
        X[a, b], X[b, a] = 1.0, -1.0
        so7.append(X)
    # X preserves phi: X_il phi_ljk + X_jl phi_ilk + X_kl phi_ijl = 0
    rows = [(np.einsum("il,ljk->ijk", X, phi) + np.einsum("jl,ilk->ijk", X, phi)
             + np.einsum("kl,ijl->ijk", X, phi)).ravel() for X in so7]
    _, s, vt = np.linalg.svd(np.array(rows).T)
    null = vt[np.sum(s > 1e-10):]
    gens = [sum(c * X for c, X in zip(vec, so7)) for vec in null]
    basis = []                                      # Gram-Schmidt in -Tr(XY)
    for X in gens:
        for B in basis:
            X = X + np.trace(X @ B) * B
        basis.append(X / math.sqrt(-np.trace(X @ X)))
    return basis


def coset_spectrum(gens):
    """Vector mass matrix per unit g^2 v^2 for a unit vacuum vector."""
    e = np.zeros(7)
    e[6] = 1.0
    cols = np.array([X @ e for X in gens])
    return np.sort(np.linalg.eigvalsh(cols @ cols.T))


def curvature(g, v, mu, k):
    """V''(v) of lambda s^4/4 + one loop vectors, lambda set by V'(v) = 0.

    k: coset eigenvalues of M^2/(g^2 s^2). MS-bar vector constant 5/6."""
    def one_loop(s):
        m2 = g**2 * s**2 * k
        return float(np.sum(3 * m2**2 / (64 * math.pi**2) * (np.log(m2 / mu**2) - 5 / 6)))
    h = 1e-4 * v
    d1 = (one_loop(v + h) - one_loop(v - h)) / (2 * h)
    lam = -d1 / v**3                                # V'(v) = lam v^3 + d1 = 0
    d2 = (one_loop(v + h) - 2 * one_loop(v) + one_loop(v - h)) / h**2
    return 3 * lam * v**2 + d2


def run(report=print):
    with contextlib.redirect_stdout(io.StringIO()):
        m = masses.derive(root.derive())
    v_EW, m_t = root.WEB.state["v_EW_GeV"], root.WEB.state["m_t"] / 1e3
    assert abs(v_EW / m["v_EW_pred_GeV"] - 1) < 1e-12
    v_b = v_EW * math.sqrt(float(root.E_v2))
    a_star = root.alpha_G2_WZW
    g = math.sqrt(4 * math.pi * a_star)
    mu_star = root.M_Pl_GeV * math.exp(-(root.S_lepton - float(root.C2_26)))

    # R1: g2, its su(3) stabiliser, the degenerate coset
    gens = g2_generators()
    spec = coset_spectrum(gens)
    assert len(gens) == 14
    assert np.allclose(spec[:8], 0, atol=1e-10)
    assert np.allclose(spec[8:], float(root.C2_7) / 6, atol=1e-10)
    k = spec[8:]
    M_V = g * v_b * math.sqrt(k[0])
    report(f"R1 g2 dim {len(gens)}, stabiliser {np.sum(spec < 1e-10)} massless (su(3)). "
           f"Coset 6 x {k[0]:.6f} g^2 v^2, M_V = {M_V:.2f} GeV. One radial mode.")

    # R2: closed form against the potential, at three renormalisation scales
    m_closed = 2 * a_star * v_b
    m_sum = math.sqrt(3 * np.sum((g * v_b)**4 * k**2) / (8 * math.pi**2 * v_b**2))
    m_pot = [math.sqrt(curvature(g, v_b, mu, k)) for mu in (50.0, v_b, mu_star)]
    assert abs(m_sum / m_closed - 1) < 1e-12
    assert all(abs(x / m_closed - 1) < 1e-5 for x in m_pot)
    report(f"R2 m_sigma = 2 alpha_G2 v_break = {m_closed:.3f} GeV. Potential curvature at "
           f"mu = 50, v_break, mu*: " + ", ".join(f"{x:.3f}" for x in m_pot))

    # R3: the scale of alpha_G2 (one loop G2 running below mu*)
    def m_at(mu):
        inv = 1 / a_star - float(root.b0_G2) / (2 * math.pi) * math.log(mu_star / mu)
        return 2 * v_b / inv
    band = [m_at(mu) for mu in (mu_star, v_b, M_V)]
    assert 34.0 < band[0] < 34.4 and 36.0 < band[1] < 37.0 and 39.3 < band[2] < 39.9
    report(f"R3 alpha_G2 at mu* = {mu_star:.1f}, v_break, M_V: m_sigma = "
           + ", ".join(f"{x:.1f}" for x in band) + f" GeV. m_sigma/M_V = {m_closed / M_V:.3f}.")

    # R4: gluon coupling at the threshold's order, and F1 (scale_ties reading B)
    ratio_lo = (21 / (12 * math.pi)) / (-1 / (3 * math.pi))
    top1, coset1, kappa1, _, _ = scale_ties.readings(1)
    _, _, kappa2, inv_alpha, v0 = scale_ties.readings(2)
    assert abs(ratio_lo + 21 / 4) < 1e-12 and abs(coset1 / top1 - ratio_lo) < 1e-6
    amp = ratio_lo * v_EW / v_b
    a_s = 1 / inv_alpha(v0, True, True, False, mu=m_closed)
    G_F = 1 / (math.sqrt(2) * v_EW**2)
    width = amp**2 * G_F * a_s**2 * m_closed**3 / (36 * math.sqrt(2) * math.pi**3)
    kg, dkg = scale_ties.KAPPA_G
    assert -7.5 < amp < -7.3
    assert not kg - 2 * dkg <= kappa1["B"] <= kg + 2 * dkg
    report(f"R4 coset/top per ln v = {ratio_lo:+.3f}. Mode amplitude vs top loop {amp:+.2f} "
           f"(rate {amp**2:.0f} x), alpha_s({m_closed:.1f}) = {a_s:.4f}, "
           f"width(gg) = {width * 1e3:.2f} MeV.")
    report(f"   F1 (portal, M_V follows h): kappa_g = {kappa1['B']:+.2f} "
           f"(two loop chain {kappa2['B']:+.2f}), measured {kg} +/- {dkg}.")

    # R5: the same form applied to h
    ratio_h = (m_t / M_Z)**4
    assert ratio_h > 0.75 and 9 * M_Z**4 - 12 * m_t**4 < 0   # M_W < M_Z
    report(f"R5 h with no scale of its own: (m_t/M_Z)^4 = {ratio_h:.1f} > 3/4, so "
           f"m_h^2 < 0. m_H stays data, kappa_lambda = 5/3 is not available to h.")
    report(f"RESULT: the G2 face radial mode sits at {band[0]:.1f} GeV ({min(band):.0f} to "
           f"{max(band):.0f} over the scale), a colour singlet that reaches the Standard "
           f"Model at one loop through gluons alone.")
    return dict(m_sigma=m_closed, band=band, M_V=M_V, amplitude=amp,
                width_GeV=width, alpha_s=a_s, kappa_F1=kappa1["B"])


if __name__ == "__main__":
    run()
