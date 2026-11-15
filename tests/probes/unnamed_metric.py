"""What an unnamed sector can build a metric from.

THE QUESTION.  The metric is the common mode of the field, and the
sector past the searchlight's horizon is amount without identity: its
channels carry no names, so any construction built from them is
invariant under every relabelling of the channels, odd ones included.
With no background metric, which first derivative fields can such a
sector use to build a symmetric rank-2 tensor?

THE CONTRACTIONS.  Without a metric the only index tool is the
Levi-Civita density.  From p-form field strengths in D = 4:
  scalars (p = 1):  h = sum_a dphi_a (x) dphi_a, needing only delta_ab.
  one forms (p = 2): chains F *G F *G ... F of odd length, *G the dual
    of a 2-form. Even chains carry one upper and one lower index.
  two form potentials (p = 3): H_mu.. eps H_nu.., the wedge of the
    2-forms H_mu and H_nu.
  three form potentials (p = 4): a density, no free index.

U1 (2-forms).  Every factor of a chain is an antisymmetric matrix, so
   the transpose of an odd chain is minus the reversed chain.  A label
   coefficient even under reversal gives a purely antisymmetric
   tensor. The symmetric part needs a coefficient odd under reversal,
   an orientation of the labels.  Repeated labels collapse at once,
   F^i *F^i = -(1/4)(F^i.*F^i) delta.  The surviving construction is
   Urbantke's, eps_ijk on an oriented triple: a named structure.
   Checked numerically for unlabelled sums at chain length three and
   five, and for eps_ijk at length three.
U2 (3-forms).  H_mu.. eps H_nu.. vanishes identically, for one field
   and for every pair.
U3 (scalars).  h from six unlabelled real components is a full rank
   symmetric tensor.
U4 (spinors).  Wetterich's metric free fermion form needs Clifford
   matrices, an internal frame: a named structure.  Stated, not
   computed.

RESULT.  A sector with no names builds a metric free symmetric tensor
from scalars alone.  This is the field level reason (the territory,
not the ledger's map) that the amount enters the stiffness with the
scalar weight: any metric built without names is scalar.  An oriented
triple can build one from 2-forms (the chiral connection route). But
an oriented triple is a name, so that structure belongs to the
nameable side, not to the amount.

Scope: first derivative fields, contractions with the Levi-Civita
density alone, chain length up to five for 2-forms.  The numbers are
random configuration checks, to 1e-12.

Usage: python3 tests/probes/unnamed_metric.py     (under a second)
"""

import itertools

import numpy as np

EPS = np.zeros((4, 4, 4, 4))
for _p in itertools.permutations(range(4)):
    EPS[_p] = np.linalg.det(np.eye(4)[list(_p)])
TOL = 1e-12


def _forms(rng, n):
    a = rng.normal(size=(n, 4, 4))
    return a - a.transpose(0, 2, 1)


def _dual(F):
    return 0.5 * np.einsum("abgd,...bg->...ad", EPS, F)


def _chain(*ms):
    out = ms[0]
    for m in ms[1:]:
        out = out @ m
    return out


def _parts(t):
    s, a = 0.5 * (t + t.T), 0.5 * (t - t.T)
    scale = max(np.abs(t).max(), 1.0)
    return np.abs(s).max() / scale, np.abs(a).max() / scale


def two_forms(seed=1, n=5):
    rng = np.random.default_rng(seed)
    F = _forms(rng, n)
    D = _dual(F)
    Fb = F.sum(0)
    rng_n = range(n)
    unlabelled = {
        "F^i *F^j F^k, all labels": sum(_chain(F[i], D[j], F[k])
                                        for i in rng_n for j in rng_n for k in rng_n),
        "F^i *F^i Fbar": sum(_chain(F[i], D[i], Fb) for i in rng_n),
        "F^i *F^j F^i *F^j Fbar": sum(_chain(F[i], D[j], F[i], D[j], Fb)
                                      for i in rng_n for j in rng_n),
        "F^i *F^j F^k *F^j F^i": sum(_chain(F[i], D[j], F[k], D[j], F[i])
                                     for i in rng_n for j in rng_n for k in rng_n),
        "distinct labels, length five": sum(_chain(F[i], D[j], F[k], D[l], F[m])
                                            for i, j, k, l, m in itertools.permutations(rng_n, 5)),
    }
    out = {k: _parts(v) for k, v in unlabelled.items()}
    e3 = np.zeros((3, 3, 3))
    for p in itertools.permutations(range(3)):
        e3[p] = np.linalg.det(np.eye(3)[list(p)])
    urb = np.einsum("ijk,abgd,ima,jbg,kdn->mn", e3, EPS, F[:3], F[:3], F[:3])
    out["Urbantke eps_ijk"] = _parts(urb)
    collapse = _chain(F[0], D[0]) + 0.25 * np.einsum("ab,ab", F[0], D[0]) * np.eye(4)
    return out, float(np.abs(collapse).max())


def three_forms(seed=3, n=3):
    rng = np.random.default_rng(seed)
    Hs = []
    for _ in range(n):
        A = rng.normal(size=(4, 4, 4))
        H = np.zeros_like(A)
        for p in itertools.permutations(range(3)):
            H += np.linalg.det(np.eye(3)[list(p)]) * A.transpose(p)
        Hs.append(H)
    X = sum(np.einsum("mab,abgd,ngd->mn", Hs[i], EPS, Hs[j]) for i in range(n) for j in range(n))
    return float(np.abs(X).max())


def scalars(seed=4, n=6):
    d = np.random.default_rng(seed).normal(size=(n, 4))
    return int(np.linalg.matrix_rank(d.T @ d))


def run(report=print):
    tf, collapse = two_forms()
    out = dict(two_forms=tf, collapse=collapse, three_forms=three_forms(), scalar_rank=scalars())
    for k, (s, a) in tf.items():
        if k == "Urbantke eps_ijk":
            assert s > 1e-2 and a < TOL, k
        else:
            assert s < TOL, k
    assert collapse < TOL
    assert out["three_forms"] < TOL
    assert out["scalar_rank"] == 4
    report("U1 2-forms: symmetric part over the unlabelled chains "
           + ", ".join(f"{k} {v[0]:.1e}" for k, v in tf.items() if k != "Urbantke eps_ijk")
           + f". The repeated label collapse holds to {collapse:.1e}; Urbantke's oriented "
           f"triple is symmetric ({tf['Urbantke eps_ijk'][0]:.2f}, antisymmetric part "
           f"{tf['Urbantke eps_ijk'][1]:.1e}).")
    report(f"U2 3-forms: H eps H over every pair = {out['three_forms']:.1e}.")
    report(f"U3 scalars: six unlabelled components give rank {out['scalar_rank']}.")
    report("RESULT: with no names, a metric free symmetric tensor is built from scalars alone; "
           "the 2-form route needs an oriented triple, a name.")
    return out


if __name__ == "__main__":
    run()
