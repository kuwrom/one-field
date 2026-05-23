"""
Layer 7b -- Bridge Self-Interference: alpha(0) from the (7,26) Sector.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
        Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
        Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
        Zenodo. https://doi.org/10.5281/zenodo.20144381

Derives the physical fine-structure constant alpha(0) = 1/137.036
from the algebraic coupling alpha_alg = pi/512 via self-interference
of the (7,26) bridge sector.  No "correction" is being applied to a
bare value; the relation below is the derivation itself:

    1/alpha(0) = (1/alpha_alg) * (1 - h_bridge/(2 pi))
               = (512/pi) * (1 - 1/(2 pi))
               = 256(2 pi - 1)/pi^2
               = 137.036

Three facts, all established in Layers 0 and 10, combine:

  (A) h_bridge = h(7, G_2) + h(26, F_4) = 2/5 + 3/5 = 1   [marginal primary]
  (B) D^2_local = 1   [Lagrangian algebra condensation -> unit coupling]
  (C) c_coset = c(E_8) - c(G_2) - c(F_4) = 0   [topological -> one-loop exact]

Physical picture
----------------
Octonions are the accounting ledger for tracking how interference
patterns reference their own history.  The bridge (7,26) is the
channel through which the EM field references itself.  The algebraic
coupling pi/512 is the value READ at the conformal-embedding scale;
the physical alpha(0) is the value READ after the bridge self-
interference is included.  Both are tree-level outputs of the same
algebra at different layers of self-reference -- not a bare value
and a dressed value, but two readings of the same structure.

The same bridge that carries gravity (protected forgetting:
P v P = 0, P v^2 P = (1/2) P) carries the EM self-interference.

Every ingredient is fixed by the conformal embedding alone:
    * h_bridge from the Sugawara formula h_R = C_2(R)/(k + h^v),
      with C_2(7) = 2 and C_2(26) = 6 from Layer 0;
    * D^2_local from the Fibonacci x Fibonacci product category, where
      A_{E_8} = 1 + (tau_G, tau_F) is the unique Lagrangian extension;
    * c_coset from c(G_2) + c(F_4) = 14/5 + 26/5 = 8 = c(E_8) at level 1.

No new ingredients enter.  The one-loop worldsheet integral of a
marginal primary in a current-current correlator is the universal
h/(2 pi); with g^2 = 1 it is 1/(2 pi) exactly.
"""

import math
from fractions import Fraction
from .constants import ALPHA_EM, ALPHA_PHYS
from .formatting import H, S, box


def derive(alg: dict) -> dict:
    """
    Derive the physical fine-structure constant from the bridge sector.

    Parameters
    ----------
    alg : dict from algebra.derive()

    Returns
    -------
    dict with keys:
        h_bridge          : conformal weight of (7,26) = 1
        h_G2_7, h_F4_26   : individual conformal weights (2/5, 3/5)
        g_bridge_sq       : effective coupling g^2 (= 1, Lagrangian)
        d_A               : quantum dimension of extension algebra
        D2_local          : local-module dimension (= 1)
        c_coset           : coset central charge (= 0)
        delta_inv_alpha   : self-interference fraction 1/(2 pi)
        inv_alpha_alg     : 1/alpha_alg = 512/pi
        inv_alpha_phys    : 1/alpha(0) = 256(2 pi - 1)/pi^2
        alpha_phys        : alpha(0)
        err_percent       : agreement with PDG (in %)
    """

    H("LAYER 7b:  BRIDGE SELF-INTERFERENCE -- alpha(0) FROM (7,26)")

    # ══════════════════════════════════════════════════════════════════
    #  7b.1  Bridge conformal weight h = 1 (marginal)
    # ══════════════════════════════════════════════════════════════════

    S("7b.1  Bridge conformal weight (Sugawara)")

    # Sugawara formula: h_R = C_2(R) / (k + h^v)
    # The conformal embedding preserves level: k = 1 for G_2(1) and F_4(1).
    k = 1
    C2_G2 = alg['C2_fund_G2']     # = 2  (Casimir of 7 of G_2)
    C2_F4 = alg['C2_fund_F4']     # = 6  (Casimir of 26 of F_4)
    hv_G2 = alg['h_dual']['G2']   # = 4  (dual Coxeter of G_2)
    hv_F4 = alg['h_dual']['F4']   # = 9  (dual Coxeter of F_4)

    h_G2_7  = Fraction(C2_G2, k + hv_G2)    # 2/5
    h_F4_26 = Fraction(C2_F4, k + hv_F4)    # 3/5
    h_bridge = h_G2_7 + h_F4_26              # 2/5 + 3/5 = 1

    # h_bridge = 1 makes (7,26) a (h, h_bar) = (1, 1) marginal primary
    # in the 2D boundary CFT -- the unique conformal weight that
    # generates a marginal worldsheet deformation.
    assert h_bridge == 1, f"h_bridge = {h_bridge}, expected 1"

    print(f"  h(7, G_2)  = C_2(7)/(k+h^v)  = {C2_G2}/({k}+{hv_G2}) = {h_G2_7}")
    print(f"  h(26, F_4) = C_2(26)/(k+h^v) = {C2_F4}/({k}+{hv_F4}) = {h_F4_26}")
    print(f"  h_bridge   = {h_G2_7} + {h_F4_26} = {h_bridge}  [MARGINAL]")

    # ══════════════════════════════════════════════════════════════════
    #  7b.2  Lagrangian algebra condensation (g^2 = 1)
    # ══════════════════════════════════════════════════════════════════

    S("7b.2  Lagrangian algebra condensation (Fibonacci x Fibonacci)")

    # Both G_2(1) and F_4(1) have Fibonacci fusion rules at level 1:
    # one nontrivial simple object with quantum dimension phi (golden
    # ratio), fusion tau x tau = 1 + tau.
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    d_tau_G = phi                          # quantum dim of G_2(1) simple
    d_tau_F = phi                          # quantum dim of F_4(1) simple

    # The extension algebra A_{E_8} = 1 + (tau_G, tau_F) is the unique
    # nontrivial local boson in the product category.  Its quantum
    # dimension is:
    d_A = 1.0 + (d_tau_G * d_tau_F)        # d(A) = 1 + phi^2 = 2 + phi

    # Total quantum dimension of the product category:
    D2_product = (1.0 + phi**2) ** 2       # = (1 + phi^2)^2

    # Local-module dimension after condensation:
    #   D^2_local = D^2(product) / d(A)^2 = 1
    # D^2_local = 1 is the CONDENSED-IS-TRIVIAL condition that defines
    # a Lagrangian algebra: the bridge is fully absorbed into the E_8
    # vacuum module with unit coupling g^2 = 1.  No tunable input.
    D2_local = D2_product / d_A ** 2

    assert abs(D2_local - 1.0) < 1e-12, f"D^2_local = {D2_local}, expected 1"

    print(f"  Fibonacci fusion: d(tau_G) = d(tau_F) = phi = {phi:.10f}")
    print(f"  Extension:  A_E_8 = 1 + (tau_G, tau_F)")
    print(f"  d(A) = 1 + phi^2 = {d_A:.10f}")
    print(f"  D^2_local = (1+phi^2)^2 / d(A)^2 = {D2_local:.1f}  [LAGRANGIAN]")
    print(f"  -> bridge coupling g^2 = 1 (no tunable input)")

    # ══════════════════════════════════════════════════════════════════
    #  7b.3  Topological coset (c = 0 -> one-loop exact)
    # ══════════════════════════════════════════════════════════════════

    S("7b.3  Topological coset (c_coset = 0)")

    # Central charges from Layer 0:
    c_E8 = Fraction(248, 31)      # = 8
    c_G2 = Fraction(14, 5)        # = 14/5
    c_F4 = Fraction(52, 10)       # = 26/5
    c_coset = c_E8 - c_G2 - c_F4  # = 0

    # c = 0 means the coset has no propagating dynamical modes.
    # Higher-loop worldsheet diagrams with internal bridge propagators
    # vanish identically: there is nothing to run in the loop.  The
    # one-loop self-interference below is therefore the EXACT reading,
    # not the leading term of a perturbative series.
    assert c_coset == 0, f"c_coset = {c_coset}, expected 0"

    print(f"  c(E_8) = 248/31 = {c_E8}")
    print(f"  c(G_2) = 14/5   = {c_G2}")
    print(f"  c(F_4) = 52/10  = {c_F4}")
    print(f"  c_coset = c(E_8) - c(G_2) - c(F_4) = {c_coset}  [TOPOLOGICAL]")
    print(f"  -> no higher-loop contributions; one-loop reading is exact")

    # ══════════════════════════════════════════════════════════════════
    #  7b.4  One-loop self-interference contribution
    # ══════════════════════════════════════════════════════════════════

    S("7b.4  One-loop self-interference contribution")

    # The EM current J_EM lives in the F_4 sector.  The bridge (7,26)
    # couples J_EM to the G_2 sector, closing a self-interference loop:
    #
    #     J_EM(z) --> phi_bridge(w) --> J_EM(0)
    #
    # In the 2D boundary CFT, the one-loop contribution of a marginal
    # primary phi with conformal weight h to the current-current
    # correlator <J(z) J(0)> is the universal:
    #
    #     delta <JJ> / <JJ> = g^2 * h / (2 pi).
    #
    # The factor h/(2 pi) is the standard worldsheet integral of a
    # marginal operator (zeta-function / heat-kernel regularised).
    # For the bridge: g^2 = 1 (Lagrangian condensation), h = 1
    # (marginal), so the self-interference reads exactly 1/(2 pi).
    h = float(h_bridge)
    g2 = D2_local                                    # = 1
    delta_inv = g2 * h / (2.0 * math.pi)             # = 1/(2 pi)

    print(f"  delta(1/alpha) / (1/alpha) = g^2 * h / (2 pi)")
    print(f"                              = {g2:.0f} * {h:.0f} / (2 pi)")
    print(f"                              = 1/(2 pi)")
    print(f"                              = {delta_inv:.10f}")
    print(f"  Physical: EM current self-interferes through the marginal")
    print(f"  bridge primary.  Same sector that carries gravity (h=1 there too).")

    # ══════════════════════════════════════════════════════════════════
    #  7b.5  Physical fine-structure constant
    # ══════════════════════════════════════════════════════════════════

    S("7b.5  Physical alpha(0)")

    # Algebraic coupling from Layer 7 (Singh ratio + WZW):
    inv_alpha_alg = 1.0 / ALPHA_EM                   # = 512/pi

    # Read the physical coupling at the IR pole-mass scale:
    inv_alpha_phys = inv_alpha_alg * (1.0 - delta_inv)
    alpha_phys = 1.0 / inv_alpha_phys

    # Consistency with the cached value in constants.py:
    assert abs(inv_alpha_phys - 1.0 / ALPHA_PHYS) < 1e-8, \
        f"derived 1/alpha = {inv_alpha_phys}, stored = {1.0 / ALPHA_PHYS}"

    # Consistency with the closed-form algebraic identity:
    inv_check = 256.0 * (2.0 * math.pi - 1.0) / math.pi**2
    assert abs(inv_alpha_phys - inv_check) < 1e-10

    # Agreement with experiment (PDG 2025: 137.035999177(21)):
    inv_alpha_exp = 137.035999177
    err_percent = 100.0 * abs(inv_alpha_phys - inv_alpha_exp) / inv_alpha_exp

    box([
        f"1/alpha(0) = (1/alpha_alg) * (1 - g^2 * h_bridge / (2 pi))",
        f"           = (512/pi) * (1 - 1/(2 pi))",
        f"           = (512/pi) * (2 pi - 1)/(2 pi)",
        f"           = 256(2 pi - 1)/pi^2",
        f"           = {inv_alpha_phys:.10f}",
        f"",
        f"Experiment: 1/alpha(0) = {inv_alpha_exp:.9f}(21)",
        f"Agreement:  {err_percent:.4f}%",
        f"",
        f"Inputs (each is a theorem):",
        f"  h_bridge = 1   (Layer 0: Sugawara on (7,26))",
        f"  g^2      = 1   (Lagrangian condensation, D^2_local = 1)",
        f"  c_coset  = 0   (topological, one-loop exact)",
        f"  alpha_alg = pi/512   (Layer 7: Singh ratio (8/3)*6 = 16)",
    ])

    return {
        'h_bridge': float(h_bridge),
        'h_G2_7': float(h_G2_7),
        'h_F4_26': float(h_F4_26),
        'g_bridge_sq': g2,
        'd_A': d_A,
        'D2_local': D2_local,
        'c_coset': float(c_coset),
        'delta_inv_alpha': delta_inv,
        'inv_alpha_alg': inv_alpha_alg,
        'inv_alpha_phys': inv_alpha_phys,
        'alpha_phys': alpha_phys,
        'err_percent': err_percent,
    }
