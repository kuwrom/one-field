"""
Layer 2 -- Lepton Masses from the G₂ sector.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

The Brannen Z₃ formula with B/A = √2, θ = 2/9:

    √mₖ = A [1 + √2 cos(2/9 + 2πk/3)]

yields all three charged-lepton pole masses from the
lepton confinement scale Λ_conf = ½ M_Pl exp(−9π²/2),
with the QED pole-mass factor (1 − α(0)/(2π)) where
α(0) = π²/(256(2π−1)) is the physical fine-structure constant,
an algebraic identity from the conformal embedding.
"""

import math
from .constants import (
    M_Pl_MeV, QED_FACTOR, PDG_MASSES,
    THETA_BRANNEN, BA_RATIO,
)
from .formatting import H, S, pct


def derive(scale: dict):
    """
    Derive the three charged-lepton masses.

    Parameters
    ----------
    scale : dict from scale.derive()

    Returns
    -------
    dict with keys:
        lepton_tree   : {e, mu, tau} tree-level masses in MeV
        lepton_pred   : {e, mu, tau} pole masses in MeV
        m_e, m_mu, m_tau : shorthand MeV values
    """

    S_0 = scale['S_0']

    H("LAYER 2:  LEPTON MASSES -- G₂ SECTOR")

    S("2.1  Brannen formula with B/A = √2, θ = 2/9")

    # The Brannen functional form (Brannen 2006) was originally an empirical
    # numerological ansatz.  In this framework B/A and theta are NOT fitted.
    # They are derived from the SU(3)_3 modular tensor category and the
    # octonionic G_2/SU(3) Clebsch-Gordan data:
    #
    #   theta = h(fund) = C_2(fund,SU(3))/(k+h^v) = (4/3)/(3+3) = 2/9
    #     (Sugawara: conformal weight of the lightest SU(3)_3 primary)
    #
    #   B/A = 2|C_3bar|/(|C_1|*d(fund)) = 2*sqrt(2)/(1*2) = sqrt(2)
    #     (octonionic Clebsch-Gordan: |C_3bar|=sqrt(2), |C_1|=1, d(fund)=2;
    #      see Innocent Lepton paper, Sec. 4)
    #
    # The Koide value Q=2/3 follows as a corollary of |B/A|=sqrt(2)
    # on the positive branch (BdG identity, Kidane 2026).
    ratio = BA_RATIO                 # = sqrt(2), derived (NOT fit)
    theta_L = THETA_BRANNEN          # = 2/9, derived (NOT fit)

    lepton_tree = {}
    lepton_labels = [(1, "e"), (2, "mu"), (0, "tau")]
    for k, name in lepton_labels:
        gap = 1.0 + ratio * math.cos(theta_L + 2.0 * math.pi * k / 3.0)
        lepton_tree[name] = 0.5 * M_Pl_MeV * math.exp(-S_0) * gap * gap

    lepton_pred = {n: m * QED_FACTOR for n, m in lepton_tree.items()}

    print(f"  √mₖ = A [1 + √2 cos(2/9 + 2πk/3)]")
    print(f"  QED correction: multiply by 1 − α/(2π) = {QED_FACTOR:.9f}")
    print()
    for name in ['e', 'mu', 'tau']:
        tree = lepton_tree[name]
        pred = lepton_pred[name]
        pdg = PDG_MASSES[name]
        err = pct(pred, pdg)
        print(f"  m_{name:<3s} = {pred:12.6f} MeV   (tree: {tree:12.6f}, PDG: {pdg:12.6f},  {err:+.3f}%)")

    return {
        'lepton_tree': lepton_tree,
        'lepton_pred': lepton_pred,
        'm_e': lepton_pred['e'],
        'm_mu': lepton_pred['mu'],
        'm_tau': lepton_pred['tau'],
    }
