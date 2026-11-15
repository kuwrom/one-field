"""Gravitationally induced entanglement: what the framework predicts for two masses in superposition.

THE QUESTION. Two masses, each in a superposition of two places and
interacting only through gravity, pick up a phase that depends on both
places. If gravity's mediator can carry the masses' joint state, the
masses entangle. If it records which place each mass is in, the
superposition decoheres. If it is a classical field, no entanglement
forms (Bose et al., arXiv:1707.06050; Marletto and Vedral,
arXiv:1707.06036). What does the framework say, and what ends it?

PREMISES. P1 records are primitive: a phase between two branches of
one closure is a record, a ratio of two clocks (bell.py,
gates_as_records.py). P2 the web holds no stochastic term: every edge
is a registered backreaction, and none reads a branch (root.py, THE
KERNEL). P3 the metric is the second order
common mode bilinear of the carrier, a quantum field, and its kernel is
the massless spin-2 pole (composite_pole.py, gate T). P4 the engine's
G (root.py). P5 the setup of Bose et al.: m = 1e-14 kg, superposition
Delta x = 250 um, trap separation d = 450 um, interaction time
tau = 2.5 s, the branch state mapped to two spins.

STATEMENTS (the numbers below check each one).

  Q1 (the phase). Each branch pair (i, j) acquires
     phi_ij = G m^2 tau/(hbar d_ij), with d_ij = d, d - Delta x,
     d + Delta x: 0.352 at d, 0.791 at the closest approach and 0.226
     at the farthest. The nonseparable part,
     Phi = phi_LL + phi_RR - phi_LR - phi_RL, is -0.314.
  Q2 (the reading). The mediator is the common mode, a composite of the
     quantum carrier (P3). The field each branch configuration produces
     is a function of the branch positions.  No quanta are radiated at
     this order, so the field returns to the vacuum when the branches
     recombine.  It keeps no record of the branch, and the branches keep
     their relative phases (P1). The masses entangle: with the engine's
     G the concurrence |sin(Phi/2)| is 0.156 of a maximum of 1.  The 55
     ppm between the engine's G and the recommended value moves it by
     under 1e-5. The concurrence is basis free. The relative phases the
     paper quotes, -0.2 and 0.7, exceed the 0.44 and -0.13 that its
     stated parameters give here.  Its witness value 1.16 rests on
     them.
  Q3 (the alternatives). A classical mediator acts by local operations
     and classical communication, which leave the product state
     separable: concurrence 0, separable by construction. A
     gravitational collapse at the Diosi-Penrose rate
     Gamma = E_Delta/hbar, with E_Delta = (12/5) G m^2/R for two uniform
     spheres of radius 1 um displaced by more than their size, gives
     tau_DP = 6.6 ms, under a three hundredth of tau: the coherence is
     gone before the phase forms. The web carries no such term (P2).

RESULT. The framework predicts gravitationally induced entanglement at
the Newtonian phase, and no gravitational decoherence.  The reason is
that the mediator is a composite of the quantum carrier, and the web
holds no stochastic term. Linearized quantized gravity makes the same
prediction.  The experiment therefore tests the framework against
classical mediators and collapse models, not against other quantum
readings of gravity.  The reading fails if a BMV class experiment, with
its environmental decoherence controlled, finds the spins separable at
the predicted phase.  It also fails on any measured gravity sourced
decoherence rate.

Scope: the static exchange between the branches at the engine's G.
Radiation from the branches' accelerations is not computed. The
Diosi-Penrose rate is the uniform sphere estimate.

Usage: python3 tests/probes/gravity_entanglement.py
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import root as r  # noqa: E402

HBAR = 1.054571817e-34
SETUP = dict(m=1e-14, dx=250e-6, d=450e-6, tau=2.5, R=1e-6)


def phases(G, m, dx, d, tau, **_):
    """phi_ij for (mass 1, mass 2) in branches L, R. Mass 1 at -/+ dx/2, mass 2 at d -/+ dx/2."""
    x1 = {"L": -dx / 2, "R": dx / 2}
    x2 = {"L": d - dx / 2, "R": d + dx / 2}
    return {(a, b): G * m * m * tau / (HBAR * abs(x2[b] - x1[a])) for a in "LR" for b in "LR"}


def spin_state(phi):
    """(|LL> e^{i phi_LL} + ...)/2 as a two qubit vector, L = up."""
    order = [("L", "L"), ("L", "R"), ("R", "L"), ("R", "R")]
    return np.array([np.exp(1j * phi[k]) for k in order]) / 2


def concurrence(psi):
    a, b, c, d = psi
    return 2 * abs(a * d - b * c)


def run(report=print):
    out = {}
    s = SETUP

    # Q1: the phases.
    phi = phases(r.G_PRED, **s)
    Phi = phi[("L", "L")] + phi[("R", "R")] - phi[("L", "R")] - phi[("R", "L")]
    assert abs(phi[("L", "L")] - 0.352) < 1e-3 and abs(Phi + 0.314) < 1e-3
    out.update(phi={"".join(k): v for k, v in phi.items()}, Phi=Phi)
    report("Q1 branch phases: " + ", ".join(f"{''.join(k)} {v:.4f}" for k, v in phi.items())
           + f". Nonseparable part {Phi:+.4f}")

    # Q2: entanglement at the engine's G, and its sensitivity to G.
    psi = spin_state(phi)
    C = concurrence(psi)
    C_one = concurrence(spin_state(phases(r.G_CODATA, **s)))
    rel = phases(r.G_PRED, **s)
    rel_RL = rel[("R", "L")] - rel[("L", "L")]
    rel_LR = rel[("L", "R")] - rel[("L", "L")]
    assert abs(C - abs(math.sin(Phi / 2))) < 1e-12 and abs(C - 0.156) < 1e-3 and abs(C - C_one) < 1e-5
    assert abs(rel_RL - 0.44) < 0.01 and abs(rel_LR + 0.126) < 0.01
    out.update(concurrence=C, C_shift=C - C_one, rel=(rel_RL, rel_LR))
    report(f"Q2 concurrence {C:.4f}. The recommended G in place of the engine's moves it by {C - C_one:.1e}; "
           f"relative phases "
           f"{rel_RL:+.3f}, {rel_LR:+.3f}")

    # Q3: the failing cases.
    separable = concurrence(spin_state({k: phi[("L", "L")] for k in phi}))     # no joint phase
    E_delta = 12 / 5 * r.G_PRED * s["m"]**2 / s["R"]
    tau_dp = HBAR / E_delta
    assert separable < 1e-12
    assert abs(tau_dp - 6.6e-3) < 1e-4 and tau_dp < s["tau"] / 300
    out.update(separable=separable, tau_dp=tau_dp)
    report(f"Q3 a separable (recording or classical) mediator: concurrence {separable:.0f}; "
           f"Diosi-Penrose tau = {1e3 * tau_dp:.2f} ms against tau = {s['tau']} s")
    return out


if __name__ == "__main__":
    run()
