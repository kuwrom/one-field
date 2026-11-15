"""A horizon on the one list: how much of its entropy is entanglement.

THE QUESTION. With G anchored a horizon of area A carries the entropy
A/(4G) = A Sigma/(8 pi t0), Sigma the one list's stiffness
(boundary.py). The rows that induce the stiffness are also entangled
across the horizon. Row by row, is the entropy the horizon carries the
entanglement of that row, and what is left over?

PREMISES. P1 the one list and its heat kernel weights per state:
scalar 1/6 - xi, Weyl 1/12, vector 1/6 - c/2 for vector contact c
(boundary.a1). P2 the entropy of a row across a flat horizon from its
heat trace on the cone of angle 2 pi alpha, under the proper time
cutoff of FACT 1:
    S = (alpha d/d alpha - 1) ln Z at alpha = 1,
    ln Z = +-(1/2) int ds/s rho(s) Tr exp(-s D),
upper sign for bosons. P3 on the cone of angle 2 pi / N the trace is a
sum over images: the k-th image of a point rotated by 2 pi k / N
contributes 1/(4 sin^2(pi k / N)) per unit transverse area over
4 pi s, times cos(pi k / N) (-1)^k per spinor component, N odd (He,
Numasawa, Takayanagi, Watanabe arXiv:1412.5606 Eqs. 2.6, 2.7, 2.15 to
2.18).

STATEMENTS (the numbers below check each one).

  H1 (scalar). (1/N) sum_k 1/(4 sin^2(pi k/N)) = (N^2 - 1)/(12 N) for
     N = 2 to 12: the cone coefficient is C(alpha) = (1/12)(1/alpha -
     alpha), and the entropy weight -C'(1) is 1/6, the heat kernel
     weight of a minimal scalar. A row at xi = 0 carries
     S = A/(48 pi t0), all of it entanglement (Solodukhin
     arXiv:1104.3712 Eqs. 26 to 28).
  H2 (spinor). With the spinor's rotation the sum is
     -(N^2 - 1)/(24 N) per component for odd N from 3 to 13: weight
     1/12 per state, the Weyl row's heat kernel weight. "For fields of
     spin zero and one half the black hole entropy is identical to the
     entropy of entanglement" (Kabat arXiv:hep-th/9503016).
  H3 (what the image sums do not hold). A coupling xi enters the cone
     through the curvature at its tip and takes xi from the weight
     without changing the trace over images: "there is no known way to
     give a statistical meaning to this entropy. Moreover, (142) does
     not have a definite sign and may become negative if xi is
     positive" (Solodukhin Sec. 3.16). A vector's two polarizations
     carry the scalar weight each, 1/6 per state (his Eq. 126), and its
     stiffness has an extra -1 per field (Eq. 148), Kabat's "contact
     interaction with the horizon which appears in the black hole
     entropy but not in the entropy of entanglement". That term "is
     the entanglement entropy of edge modes", and "the sign of the
     leading order divergence is not universal under a change of
     regulator scheme" (Donnelly and Wall arXiv:1412.1895).  The
     contact coefficient c of the list measures that spread.
  H4 (the one list). The entanglement of the list is
     182/6 + 90/12 + 36/6 = 43.8333, the stiffness at contact c = 0.
     As a fraction of A/(4G), set by the stiffness, it is 1.70, 1.00
     and 0.71 at the contacts +1, 0, -1.  At the covariant contact the
     vectors' contact term, -18, is the whole difference.
  H5 (the amount). Minimally coupled, the amount's 30.33 is
     entanglement entirely, 69 percent of the list's. With the label's
     xi = 1/(48 pi) its weight would contain 1.21 that no image sum
     supplies. The nameable rows alone carry entanglement 13.5 against
     a stiffness of -4.5 at the covariant contact: a horizon of
     negative entropy, which the amount makes positive.

RESULT. On the one list a horizon's entropy is the entanglement of the
rows across it plus the vectors' edge term.  Every scalar and Weyl row
contributes to the stiffness exactly what it entangles.  The amount is
included because it is minimally coupled.  The vector contact
coefficient, the edge modes' regulator dependent share, is the one
place the two counts differ. The entropy's value is A/(4G) by the
anchor. No value of G is used here.

Scope: flat horizon, free rows, the leading area term. The
continuation from integer N to alpha is the replica continuation.

Usage: python3 tests/probes/horizon_entropy.py     (instant)
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import boundary as B  # noqa: E402
from interference import root as r  # noqa: E402

XI = float(r.xi_bridge)                          # the label's coupling
ENTANGLE = {"vector": 1 / 6, "weyl": 1 / 12, "scalar": 1 / 6, "bridge": 1 / 6}   # per state, from H1-H3


def image_sum(N, spinor=False):
    """Cone coefficient at angle 2 pi / N per component: (1/N) sum over the N - 1 images."""
    k = np.arange(1, N)
    weight = np.cos(np.pi * k / N) * (-1.0)**k if spinor else 1.0
    return float(np.sum(weight / (4 * np.sin(np.pi * k / N)**2)) / N)


def cone_weight(Ns, spinor=False):
    """Fit N C(N) = p N^2 + q N + r, read C(alpha) = p/alpha + q + r alpha, return the entropy weight.

    Bosons: -C'(1) = p - r. A fermion loop enters ln Z with the other sign.
    """
    data = [N * image_sum(N, spinor) for N in Ns]
    p, q, c = np.polyfit(Ns, data, 2)
    resid = max(abs(p * N * N + q * N + c - d) for N, d in zip(Ns, data))
    sign = -1.0 if spinor else 1.0
    return sign * (p - c), (p, q, c), resid


def run(report=print):
    out = {}

    # H1, H2: the image sums and the weights they continue to.
    scalar_N, spinor_N = list(range(2, 13)), list(range(3, 14, 2))
    assert max(abs(image_sum(N) - (N * N - 1) / (12 * N)) for N in scalar_N) < 1e-12
    assert max(abs(image_sum(N, True) + (N * N - 1) / (24 * N)) for N in spinor_N) < 1e-12
    w0, fit0, res0 = cone_weight(scalar_N)
    wh, fith, resh = cone_weight(spinor_N, spinor=True)
    assert abs(w0 - 1 / 6) < 1e-10 and abs(wh - 1 / 12) < 1e-10 and max(res0, resh) < 1e-10
    assert abs(w0 - B.a1("scalar", 0.0, 0.0)) < 1e-10 and abs(wh - B.a1("weyl", 0.0, 0.0)) < 1e-10
    out.update(scalar_weight=w0, spinor_weight=wh)
    report(f"H1 scalar: image sums equal (N^2 - 1)/(12 N). Entropy weight {w0:.6f} = 1/6.  "
           f"H2 spinor: -(N^2 - 1)/(24 N) per component. Weight {wh:.6f} = 1/12 per state")

    # H3, H4: the list's entanglement against the stiffness the horizon carries.
    one = B.rows(r.d10, r.d11, r.n7, r.n26)
    live = [(name, n, kind) for name, n, kind, _, reading in one if reading in ("light", "amount")]
    ent = sum(n * ENTANGLE[kind] for _, n, kind in live)
    sig = {c: B.sigma(one, 0.0, c, XI) for c in B.CONTACTS}
    share = {c: ent / s for c, s in sig.items()}
    vectors = sum(n for _, n, kind in live if kind == "vector")
    assert abs(ent - sig[0.0]) < 1e-12 and abs(ent - 43.8333333) < 1e-6
    assert abs((ent - sig[1.0]) - vectors / 2) < 1e-12 and vectors / 2 == 18
    out.update(entanglement=ent, share=share, contact_term=-(ent - sig[1.0]))
    report(f"H4 one list: entanglement {ent:.4f} = Sigma at contact 0. Over the horizon's Sigma it is "
           + ", ".join(f"{share[c]:.3f}" for c in B.CONTACTS) + " at contacts +1, 0, -1. The vectors' contact "
           f"term is {out['contact_term']:+.0f} at the covariant contact")

    # H5: the amount, and the nameable rows alone.
    amount = sum(n * ENTANGLE[kind] for _, n, kind in live if kind == "bridge")
    amount_sigma = sum(n * B.a1(kind, 1.0, B.xi_of("amount", XI)) for _, n, kind in live if kind == "bridge")
    label_gap = 182 * XI
    named_ent, named_sig = ent - amount, sig[1.0] - amount_sigma
    assert abs(amount - amount_sigma) < 1e-12 and abs(label_gap - 1.2069) < 1e-3
    assert named_ent > 0 > named_sig and abs(named_ent - 13.5) < 1e-12 and abs(named_sig + 4.5) < 1e-12
    out.update(amount=amount, amount_share=amount / ent, label_gap=label_gap, named=(named_ent, named_sig))
    report(f"H5 amount: {amount:.4f}, all entanglement, {100 * amount / ent:.0f}% of the list's. With the label's "
           f"xi its weight would contain {label_gap:.3f} that no image sum supplies. The nameable rows alone entangle "
           f"{named_ent:.1f} against a stiffness of {named_sig:+.1f}")
    report("RESULT: a horizon's entropy on the one list is the entanglement of its rows plus the vectors' edge "
           "term. Scalar and Weyl rows, the amount included, weigh what they entangle.")
    return out


if __name__ == "__main__":
    run()
