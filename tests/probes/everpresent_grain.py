"""Lambda option D3 against the data: a record count per grain.

THE MODEL. Everpresent Lambda (Ahmed, Dodelson, Greene, Sorkin,
astro-ph/0209274; Zwane, Afshordi, Sorkin arXiv:1703.06265 Eqs. 5-8):
the vacuum term fluctuates with the Poisson count of elements,
rho_Lambda(t+1) = [S(t) + beta xi sqrt(N(t+1) - N(t))] / V(t), with
N = V / l_p^4 in units 8 pi G = hbar = c = 1.  The data need
0.01 < beta < 0.02 (1703.06265, citing Ahmed et al.). If beta is
larger, the total density turns negative before the universe reaches
its size. If it is smaller, the universe is younger than its oldest
galaxies.  The
heuristic causal set value is of order one per element.

THE FRAMEWORK'S GRAIN. With G anchored the regulator is
t0 = Sigma/(2 pi) l_Pl^2 = 4.112 l_Pl^2 on the one list
(boundary.py), l_Pl^2 = G hbar.  Take one record per grain 4-volume
v_g = k t0^2, where k is a shape factor (1 for t0^2, pi^2/2 for the
4-ball of radius sqrt(t0)). The count is then N = V / v_g. With an
order one action per record, the amplitude in the paper's units is
beta = 1 / sqrt(v_g / l_p^4).

E1. In the paper's stated units (l_p^2 = 8 pi G hbar) the grain is
    v_g = k (4.112 / 8 pi)^2 l_p^4 = 0.0268 k l_p^4, and beta is 6.1,
    2.8 and 1.0 for k = 1, pi^2/2 and 40: above the window by a factor
    of 50 to 600.
E2. With the unreduced length (l_p^2 = G hbar), v_g = 16.9 k and
    beta = 0.24 and 0.11 for k = 1 and pi^2/2: above the window by 5 to
    25.
E3. The window needs a grain 4-volume of 2.5e3 to 1e4 l_p^4, or an
    action per record near 1e-2: a new small number either way.

RESULT. With the anchor's Planckian grain and an order one action per
record, the record count fluctuates too strongly for the data in
either convention. The density turns negative before today.  D3 in
this form is excluded: it needs a small factor the state does not supply.  The
estimate carries the shape factor k and the action per record as
order one conventions. No verdict here rests on better than a factor
of a few.

Usage: python3 tests/probes/everpresent_grain.py     (instant)
"""

import math

T0_LPL2 = 4.112                      # anchored regulator, units of G hbar (boundary.py)
WINDOW = (0.01, 0.02)                # beta allowed by the data (1703.06265 Model 1)
SHAPES = {"t0^2": 1.0, "4-ball": math.pi**2 / 2, "k = 40": 40.0}


def beta(k, reduced=True):
    t0 = T0_LPL2 / (8 * math.pi) if reduced else T0_LPL2
    return 1.0 / math.sqrt(k * t0 * t0)


def run(report=print):
    out = {"reduced": {s: beta(k) for s, k in SHAPES.items()},
           "unreduced": {s: beta(k, reduced=False) for s, k in SHAPES.items() if k < 10}}
    need = (1 / WINDOW[1]**2, 1 / WINDOW[0]**2)          # grain 4-volume in l_p^4
    out["grain_needed"] = need
    for conv in ("reduced", "unreduced"):
        for s, b in out[conv].items():
            assert b > WINDOW[1], (conv, s)
    report("E1 reduced units (l_p^2 = 8 pi G): beta = "
           + ", ".join(f"{b:.2f} ({s})" for s, b in out["reduced"].items())
           + f". The data allow {WINDOW[0]} to {WINDOW[1]}.")
    report("E2 unreduced (l_p^2 = G): beta = "
           + ", ".join(f"{b:.2f} ({s})" for s, b in out["unreduced"].items()) + ".")
    report(f"E3 the window needs a grain 4-volume of {need[0]:.0f} to {need[1]:.0f} l_p^4, against "
           f"{(T0_LPL2 / (8 * math.pi))**2:.4f} k from the anchor.")
    report("RESULT: the anchor's grain makes the record count fluctuate too strongly. D3 in the "
           "everpresent form needs an unexplained small factor.")
    return out


if __name__ == "__main__":
    run()
