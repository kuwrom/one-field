"""Formation on the self-bound sign: a random cloud in the empty field.

lump_in_empty_field.py and lump_stability.py solve the lump and its
second variation. This probe runs the formation direction on the same
sign. A random three component cloud is released in the empty field
under the full real time dynamics, and the probe measures what is left.

PREMISES. The nonrelativistic limit of the carrier, hbar = m = 1,
    i d_t psi_a = -(1/2) psi_a'' + (M[psi] psi)_a + kappa R^2 psi_a,
M the interaction matrix of tests/substrate/nls_soliton.py, R the total
density, on the self-bound sign g0 = -1, g1 = 1/sqrt 2, kappa = 1. One
space dimension, a periodic box of length 100 with an absorbing layer
beyond |x| = 36, so radiated waves leave the box. Outside is the empty
field.
Initial data: a Gaussian envelope of width 4 times smooth band limited
random complex amplitudes in all three components. The data are made
even in x (zero momentum, no recoil, so a single lump stays at the
centre) and carry three times the charge N1 = 1.26359 of the reference
lump. The integrator is fourth order Runge-Kutta on the interaction,
with the kinetic flow treated exactly (integrating factor). In a closed box it holds N to 1e-7 and E to 1e-6
over t = 40.

MEASURES. In the window |x - x_peak| < 8: the charge N_w, the energy
E_w, the peak density, the Z3 Fourier shares of
psi_hat_q = sum_a conj(w)^(q a) psi_a / 3, and the lock
    lam = -6 int sum_k rho_k Re(conj(psi_{k+1}) psi_{k+2}) / int R^2,
the phase lock energy in units of its clock value. lam = 1 on either
clock orientation (q = 1 or q = 2), and lam = -2 on the symmetric state.

STATEMENTS (the numbers below check each one).

  F1 (a lump forms). By t = 100 the cloud has shed about two fifths of
     its charge through the layer and a bound lump sits at the centre:
     N_w = 2.1, E_w = -0.22 and -0.26, peak density 0.54 and 0.59
     (seeds 1 and 2), against 0.61 for the closed form clock lump at
     that charge.
  F2 (the G2 term turns the lump toward the clock). The lock rises
     from the cloud's -0.64 and -0.50 to +0.30 and +0.80 over
     100 <= t <= 300 and to +0.41 and +0.83 over 300 <= t <= 500. At
     g1 = 0, same data, a lump still forms (the common term binds) and
     its Z3 shares and lock are frozen: the two windows agree to 1e-2.
  F3 (the clock does not lock in real time). lam stays below one, and
     the lump stays E_w - E_clock(N_w) = +0.06 and +0.02 above the
     clock lump at its own charge. The excess sits in relative modes
     bound inside the lump, below the continuum edge
     (lump_stability.py S1). Such modes cannot radiate at linear order.
     The minimum at fixed N is the clock lump, and gradient flow
     reaches it (lump_stability.py S4). Conservative real time flow
     approaches it on this time horizon but does not arrive. The lump's
     soft relative mode reaches the continuum only at its seventh
     harmonic and does not radiate. A large excess drains, and the
     draining slows (lump_internal_modes.py). The dynamics is far from
     the integrable Manakov point.
  F4 (the gas sign disperses). With g0 = +1 and no sextic, same data,
     the cloud spreads into the layer: late peak density below 0.02
     and no bound object.
  F5 (not every cloud gives one lump). With the same settings, seed 4
     sends its charge into the layer as travelling lumps.

RECORD, eight seeds at T = 600 with these settings (not rerun here).
Seven give a centred lump, seed 4 does not. Lock of the seven, as
cloud -> 100..300 -> 300..600: seed 1 -0.64, +0.29, +0.41. Seed 2 -0.50,
+0.80, +0.83. Seed 3 -0.13, -0.39, -0.30. Seed 5 +0.36, +0.52, +0.52.
Seed 6 -0.04, +0.07, +0.20. Seed 7 -0.75, +0.37, +0.59. Seed 8 -0.53,
+0.23, +0.32. The lock does not fall between the two late windows in
any of the seven and ends above its cloud value in six. The same
seven at g1 = 0 are frozen and take either sign. Excess over the clock
lump: +0.02 to +0.06.

Three dimensions are not run here. RECORD (Kaggle kernel
scratchpad/kaggle/formation3d, 80^3, complex128, box 40, absorbing
beyond r = 14, cloud of charge 300, t = 2000): both seeds leave a
centred bound lump (N_w = 150, 191, E_w < 0, peak 0.98, 0.96). The
lock rises from the cloud's -0.25 and +0.28 to +0.66 and +0.48 and then
holds. The excess over the lattice clock lump at the same charge stays
at +3.07 -> +3.05 and +6.18 -> +6.05 from t = 500 to 2000, the Z3 shares
settle near (0.25, 0.37, 0.37), and the g1 = 0 control is frozen. The
registered prediction, a soft mode excess draining as (1 + Gamma t)^(-1/3)
(the 3D soft mode sits at 0.28 of the continuum edge), is not seen on
this horizon. A 48^3 reproduction of seed 1 ends on one clock orientation
(q = 1 in every interior cell, interior lock 0.90): the content of the
excess is resolution dependent and is not settled by these runs.

RECORD, a fluctuating field in a closed box (Kaggle kernel
scratchpad/kaggle/universe3d, 128^3, box 80, Gaussian random field in
all three components at mean density 0.08, no absorber, nothing held,
t = 2000; charge conserved to 6e-7 and energy to 6e-6). The field clumps
into 92 lumps by t = 125, which merge to 20 by t = 2000 (largest
N = 4077) while 46% of the charge sits in lumps and the rest stays a
dilute background. The charge weighted lock rises 0.52, 0.64, 0.74,
0.83, 0.87, 0.91 at t = 200, 400, 600, 1000, 1600, 2000. So in a
populated box the lumps keep turning toward the clock as they merge,
whereas the isolated lump stalls. Each large lump sits on one clock orientation
(0.93 to 0.95 of its charge in q = 1 or in q = 2) and the two
orientations come in equal numbers (41/46 at t = 125, 9/11 at
t = 2000). One lump (N = 2220) holds both, 0.63 and 0.33. No lump
spins: |L|/N stays below 0.25 and ends below 0.06 for every lump above
N = 300.

RECORD, the same field at seed 2 run to t = 4000 with orientation
domains and vortex plaquettes counted inside each lump (kernel
kuwrom/z3nls-universe3d-2; charge conserved to 5e-6, energy to 1e-4).
The lumps coarsen from 87 to 8 while their share of the charge grows to
62%, and the charge weighted lock reaches 0.94. The orientations stay
balanced in number (41:43 at t = 400, 3:5 at t = 4000) but not in
charge. Mergers pass the lumps to the larger orientation, so the charge
in the two goes 0.57:0.41, 0.22:0.78, 0.17:0.83 at t = 400, 2000, 4000.
Seed 1 drifts the other way (0.66:0.34 at t = 2000). Lumps holding both
orientations number 25 at t = 400 and none from t = 2800. Vortices
appear inside large lumps during mergers (544 plaquettes in a lump of
charge 4532 at t = 2000) and are gone by t = 2800. No lump keeps spin
(|L|/N at most 0.16, 0.04 at t = 4000). Every lump at t = 4000 sits on
one orientation with lock 0.93 to 0.95.

Usage: python3 tests/probes/lump_formation.py     (about fifteen seconds)
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from nls_soliton import G1_OVER_G0, interaction_matrix  # noqa: E402
import lump_in_empty_field as lf  # noqa: E402
import lump_stability as ls  # noqa: E402

W3 = lf.W3
N_REF = 1.26359                       # charge of the 1D lump at mu = 0.6 mu0
BOX, POINTS, DT = 100.0, 512, 0.04
WINDOW = 8.0


def force(psi, g0, g1, kappa):
    """M[psi] psi + kappa R^2 psi, written out component by component."""
    rho = abs(psi)**2
    R = rho.sum(axis=0)
    out = np.empty_like(psi)
    for m in range(3):
        a, b = (m + 1) % 3, (m + 2) % 3
        out[m] = ((g0 * R + kappa * R * R + g1 * np.real(np.conj(psi[a]) * psi[b])) * psi[m]
                  + 0.5 * g1 * (rho[b] * psi[a] + rho[a] * psi[b]))
    return out


def grid():
    dx = BOX / POINTS
    x = (np.arange(POINTS) - POINTS // 2) * dx
    return x, 2 * np.pi * np.fft.fftfreq(POINTS, d=dx), dx


def cloud(N, seed, width=4.0, kcut=1.5):
    """Smooth random three component cloud, even in x, total charge N."""
    x, k, dx = grid()
    rng = np.random.default_rng(seed)
    white = rng.normal(size=(3, POINTS)) + 1j * rng.normal(size=(3, POINTS))
    psi = np.exp(-(x / width)**2) * np.fft.ifft(np.fft.fft(white, axis=1) * (abs(k) < kcut),
                                                axis=1)
    psi = psi + np.roll(psi[:, ::-1], 1, axis=1)              # psi(x) = psi(-x)
    return psi * math.sqrt(N / (np.sum(abs(psi)**2) * dx))


def measure(psi, g0, g1, kappa):
    x, k, dx = grid()
    rho = abs(psi)**2
    R = rho.sum(axis=0)
    x_peak = x[np.argmax(R)]
    win = abs(x - x_peak) < WINDOW
    px = np.fft.ifft(1j * k * np.fft.fft(psi, axis=1), axis=1)
    exchange = sum(rho[j] * np.real(np.conj(psi[(j + 1) % 3]) * psi[(j + 2) % 3])
                   for j in range(3))
    e = 0.5 * np.sum(abs(px)**2, axis=0) + 0.5 * g0 * R**2 + kappa * R**3 / 3 + g1 * exchange
    modes = np.array([np.sum(np.conj(W3)[:, None]**q * psi, axis=0) / 3 for q in range(3)])
    weight = np.sum(abs(modes[:, win])**2, axis=1)
    return dict(N=R.sum() * dx, E=e.sum() * dx, N_w=R[win].sum() * dx, E_w=e[win].sum() * dx,
                peak=R.max(), x_peak=x_peak, shares=weight / weight.sum(),
                lock=-6 * exchange[win].sum() / (R[win]**2).sum())


def evolve(g0, g1, kappa, seed, T, absorb=True, samples=30):
    """Integrating factor RK4 from the random cloud. Samples at equal times."""
    x, k, dx = grid()
    psi = cloud(3 * N_REF, seed)
    half = np.exp(-0.25j * k * k * DT)                         # exact kinetic half step
    edge = np.clip((abs(x) - 0.36 * BOX) / (0.14 * BOX), 0, 1)
    mask = np.exp(-DT * edge**2) if absorb else 1.0

    def flow(u):
        return np.fft.ifft(np.fft.fft(u, axis=1) * half, axis=1)

    def rate(u):
        return -1j * force(u, g0, g1, kappa)
    steps = samples * int(round(T / DT / samples))
    out = [dict(t=0.0, **measure(psi, g0, g1, kappa))]
    for i in range(steps):
        k1 = rate(psi)
        a = flow(psi)
        k2 = rate(a + 0.5 * DT * flow(k1))
        k3 = rate(a + 0.5 * DT * k2)
        b = flow(a)
        k4 = rate(b + DT * flow(k3))
        psi = (b + DT / 6 * (flow(flow(k1)) + 2 * flow(k2) + 2 * flow(k3) + k4)) * mask
        if (i + 1) % (steps // samples) == 0:
            out.append(dict(t=(i + 1) * DT, **measure(psi, g0, g1, kappa)))
    return out


def window(out, t_lo, t_hi):
    """Averages over the samples with t_lo < t <= t_hi."""
    sel = [o for o in out if t_lo < o["t"] <= t_hi]
    return dict(N_w=float(np.mean([o["N_w"] for o in sel])),
                E_w=float(np.mean([o["E_w"] for o in sel])),
                peak=float(np.mean([o["peak"] for o in sel])),
                lock=float(np.mean([o["lock"] for o in sel])),
                shares=np.mean([o["shares"] for o in sel], axis=0),
                off_centre=float(max(abs(o["x_peak"]) for o in sel)),
                kept=sel[-1]["N"] / out[0]["N"])


def clock_lump(N):
    """Energy and peak density of the closed form 1D clock lump at charge N."""
    mu, E = ls.lump_1d(N)
    return E, float(lf.profile_1d(mu)[0](0.0))**2


def run(report=print):
    rng = np.random.default_rng(0)
    p = rng.normal(size=(3, 9)) + 1j * rng.normal(size=(3, 9))
    R = np.sum(abs(p)**2, axis=0)
    force_err = np.max(abs(force(p, -1.0, G1_OVER_G0, 1.0) - R * R * p
                           - np.einsum("xab,bx->ax", interaction_matrix(p, -1.0, G1_OVER_G0), p)))
    assert force_err < 1e-12

    closed = evolve(-1.0, G1_OVER_G0, 1.0, 1, T=40.0, absorb=False, samples=4)
    dN = abs(closed[-1]["N"] / closed[0]["N"] - 1)
    dE = abs(closed[-1]["E"] / closed[0]["E"] - 1)
    assert dN < 1e-6 and dE < 1e-5

    # F1-F3: the self-bound sign, seeds 1 and 2.
    main = {}
    for seed in (1, 2):
        out = evolve(-1.0, G1_OVER_G0, 1.0, seed, T=500.0, samples=100)
        early, late = window(out, 100.0, 300.0), window(out, 300.0, 500.0)
        E_clock, peak_clock = clock_lump(late["N_w"])
        a = dict(start=out[0]["lock"], early=early, late=late, E_clock=E_clock,
                 peak_clock=peak_clock, excess=late["E_w"] - E_clock)
        assert early["off_centre"] < 1.0 and late["off_centre"] < 1.0
        assert late["N_w"] > N_REF and abs(late["N_w"] / early["N_w"] - 1) < 0.05
        assert late["E_w"] < -0.1 and 0.4 < late["peak"] < peak_clock
        assert early["lock"] > a["start"] + 0.5 and late["lock"] > early["lock"] - 0.05
        assert late["lock"] < 0.95 and 0.005 < a["excess"] < 0.1
        main[seed] = a

    # F2: the phase lock switched off.
    out = evolve(-1.0, 0.0, 1.0, 2, T=300.0, samples=60)
    free_a, free_b = window(out, 100.0, 200.0), window(out, 200.0, 300.0)
    assert free_b["E_w"] < -0.1 and free_b["peak"] > 0.4 and free_b["off_centre"] < 1.0
    assert abs(free_b["lock"] - free_a["lock"]) < 0.02
    assert np.max(abs(free_b["shares"] - free_a["shares"])) < 0.01

    # F4: the gas sign.
    gas = window(evolve(1.0, G1_OVER_G0, 0.0, 2, T=150.0, samples=15), 100.0, 150.0)
    assert gas["peak"] < 0.02 and gas["E_w"] > 0 and gas["kept"] < 0.3

    # F5: a cloud that does not leave one lump.
    split = window(evolve(-1.0, G1_OVER_G0, 1.0, 4, T=160.0, samples=16), 120.0, 160.0)
    assert split["off_centre"] > 10.0

    report(f"Force against the interaction matrix: {force_err:.1e}. Closed box, t = 40: "
           f"dN/N = {dN:.1e}, dE/E = {dE:.1e}.")
    for seed, a in main.items():
        e, l = a["early"], a["late"]
        report(f"F1 seed {seed}: N_w = {l['N_w']:.3f}, E_w = {l['E_w']:+.4f}, peak "
               f"{l['peak']:.3f} ({100 * l['kept']:.0f}% of the charge still in the box); "
               f"clock lump at that N: E = {a['E_clock']:+.4f}, peak {a['peak_clock']:.3f}.")
        report(f"F2 seed {seed}: lock {a['start']:+.2f} (cloud) -> {e['lock']:+.2f} "
               f"(100..300) -> {l['lock']:+.2f} (300..500). Shares (q = 0, 1, 2) "
               f"{np.round(l['shares'], 2)}.")
        report(f"F3 seed {seed}: E_w - E_clock = {a['excess']:+.4f}.")
    report(f"F2 g1 = 0: N_w = {free_b['N_w']:.3f}, E_w = {free_b['E_w']:+.4f}, peak "
           f"{free_b['peak']:.3f}. Lock {free_a['lock']:+.3f} -> {free_b['lock']:+.3f}, "
           f"shares {np.round(free_a['shares'], 3)} -> {np.round(free_b['shares'], 3)}.")
    report(f"F4 gas sign: late peak {gas['peak']:.4f}, N_w = {gas['N_w']:.3f}, E_w = "
           f"{gas['E_w']:+.4f}, {100 * gas['kept']:.0f}% of the charge still in the box.")
    report(f"F5 seed 4: the peak sits {split['off_centre']:.1f} from the centre by t = 160.")
    report("RESULT: on the self-bound sign a random cloud in the empty field leaves a "
           "bound lump. The G2 term turns it toward the clock: the lock rises from the "
           "cloud's value and does not fall back. At g1 = 0 the lump still binds, but its "
           "internal state stays frozen. Conservative real time flow does not complete "
           "the lock on this time horizon. The excess over the clock lump stays in relative "
           "modes bound inside it. The gas sign disperses.")
    return dict(force_error=force_err, closed_dN=dN, closed_dE=dE, main=main,
                free=(free_a, free_b), gas=gas, split=split)


if __name__ == "__main__":
    run()
