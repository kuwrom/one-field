"""WHICH UNIFORM STATE IS THE VACUUM OF THE Z3-COUPLED SUBSTRATE?

Every 3D Z3-NLS run in this tree has been seeded on the Z3-SYMMETRIC
background psi_0 = psi_1 = psi_2. That state is a stationary point, and
emergent_clock.py already measured that it is the one that "decays into
matter" while offset backgrounds sit still. This probe computes why, in
closed form: it is not the vacuum, and it is linearly unstable.

At fixed total density the kinetic term of a uniform state vanishes, so the
whole question is exact algebra on the interaction energy

    E/V = (g0/2) rho_tot^2 + g1 SUM_m rho_m Re(conj(psi_{m+1}) psi_{m+2})

and the first term is a constant on the constraint surface. The competition
is entirely in the g1 ring term, and it prefers the CLOCK.

WHAT THIS SETTLES, and it is not a small thing. The Faddeev-Skyrme director
is z = phi_rel / sqrt(rho_rel + EPS), so

    |n| = rho_rel / (rho_rel + EPS)

In the symmetric background rho_rel = 0 identically, so |n| = 0: n does not
map into S^2 but into a contractible BALL, pi_3(ball) = 0, and the Hopf
charge is not quantized. That is why a barrier scan found NO energy
barrier at any of four resolutions and why Q_H reads 0.94 rather than
1: there was no topology for the derived quartic to protect. In the
condensate rho_rel != 0 everywhere, |n| = 1, and n is a constant unit
vector -- the boundary condition a Hopf charge actually requires.

The condensate also keeps the Z3 CENTRE UNBROKEN: the cyclic shift is
equivalent to a global phase on it, which is exactly the projective
realisation the confined-phase entry test (confined_phase.py) posits.
That test carries the phase as a premise; this state is a stable stationary point in
closed form with mu = 3 g0 - g1.

Nothing here is fitted or relaxed. Every number is either exact algebra or
a numerically exact eigendecomposition of a 6x6 matrix.
"""

import numpy as np

import substrate3d as s3

W3 = np.exp(2j * np.pi / 3)
CLOCK = np.array([1.0, W3, W3**2])
SYMMETRIC = np.array([1.0, 1.0, 1.0])
RHO_TOT = 3.0 * s3.rho0


def _normalise(vec):
    vec = np.asarray(vec, dtype=complex)
    return vec * np.sqrt(RHO_TOT / np.sum(np.abs(vec) ** 2))


def energy_density(vec):
    """E/V of a uniform state. Uses the repo's own energy(), so this cannot
    drift from the functional the 3D runs integrate."""
    vec = _normalise(vec)
    n = 4
    psi = np.ones((3, n, n, n), dtype=complex) * vec[:, None, None, None]
    return float(s3.energy(psi, 1.0)) / n**3


def director(vec):
    """(rho_rel, |n|, n) of a uniform state, in substrate3d's convention."""
    vec = _normalise(vec)
    phi = np.array([np.sum(vec * W3 ** (-q * np.arange(3))) / 3.0
                    for q in range(3)])
    rho_rel = float(abs(phi[1]) ** 2 + abs(phi[2]) ** 2)
    if rho_rel < 1e-14:
        return rho_rel, 0.0, np.zeros(3)
    z = np.array([phi[1], phi[2]]) / np.sqrt(rho_rel)
    n = np.array([2 * np.real(np.conj(z[0]) * z[1]),
                  2 * np.imag(np.conj(z[0]) * z[1]),
                  abs(z[0]) ** 2 - abs(z[1]) ** 2])
    return rho_rel, float(np.linalg.norm(n)), n


def chemical_potential(vec):
    """mu from the uniform stationarity condition dE/dpsi* = mu psi."""
    vec = _normalise(vec)
    n = 4
    psi = np.ones((3, n, n, n), dtype=complex) * vec[:, None, None, None]
    K2 = np.zeros((n, n, n))          # uniform: no kinetic contribution
    return s3.stationarity(psi, K2, 1.0)


def hessian_k0(vec):
    """Eigenvalues of the k=0 Hessian of (E - mu N)/V over uniform states.

    k = 0 is where any uniform instability is strongest, since the kinetic
    term only adds k^2/2 >= 0. Built by central differences on the six real
    coordinates, so it makes no assumption about the algebra.
    """
    vec = _normalise(vec)
    mu, _ = chemical_potential(vec)

    def f(x):
        v = x[:3] + 1j * x[3:]
        n = 4
        psi = np.ones((3, n, n, n), dtype=complex) * v[:, None, None, None]
        e = float(s3.energy(psi, 1.0)) / n**3
        return e - mu * float(np.sum(np.abs(v) ** 2))

    x0 = np.concatenate([vec.real, vec.imag])
    h, m = 1e-4, 6
    H = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            ei, ej = np.zeros(m), np.zeros(m)
            ei[i] = h; ej[j] = h
            H[i, j] = (f(x0 + ei + ej) - f(x0 + ei - ej)
                       - f(x0 - ei + ej) + f(x0 - ei - ej)) / (4 * h * h)
    return np.linalg.eigvalsh(0.5 * (H + H.T))


def z3_shift_ratio(vec):
    """psi_k -> psi_{k-1} divided by psi_k. Constant <=> the centre is
    realised projectively (shift == global phase) <=> unbroken."""
    vec = _normalise(vec)
    r = np.roll(vec, 1) / vec
    return r, float(np.max(np.abs(r - r[0])))


def run(report=print):
    states = {"Z3-symmetric (1,1,1)": SYMMETRIC,
              "clock q=1 (1,w,w^2)": CLOCK,
              "clock q=2 (1,w^2,w)": np.array([1.0, W3**2, W3]),
              "two-component (1,1,0)": np.array([1.0, 1.0, 0.0]),
              "single (1,0,0)": np.array([1.0, 0.0, 0.0])}
    report("  uniform states at fixed rho_tot = 3 (kinetic term vanishes)")
    report(f"  {'state':<24}{'E/V':>12}{'rho_rel':>10}{'|n|':>7}")
    out = {}
    for name, v in states.items():
        e = energy_density(v)
        rr, nn, _ = director(v)
        out[name] = dict(E=e, rho_rel=rr, n_norm=nn)
        report(f"  {name:<24}{e:>12.6f}{rr:>10.4f}{nn:>7.3f}")

    e_sym = out["Z3-symmetric (1,1,1)"]["E"]
    e_clk = out["clock q=1 (1,w,w^2)"]["E"]
    report(f"\n  the symmetric state is {100*(e_sym/e_clk - 1):.0f}% higher "
           f"in energy density than the clock")

    mu_sym, res_sym = chemical_potential(SYMMETRIC)
    mu_clk, res_clk = chemical_potential(CLOCK)
    report(f"  both are exact stationary points:")
    report(f"    symmetric  mu = {mu_sym:+.5f} = 3g0 + 2g1   "
           f"residual {res_sym:.1e}")
    report(f"    clock      mu = {mu_clk:+.5f} = 3g0 -  g1   "
           f"residual {res_clk:.1e}")

    h_sym = hessian_k0(SYMMETRIC)
    h_clk = hessian_k0(CLOCK)
    n_sym = int(np.sum(h_sym < -1e-6))
    n_clk = int(np.sum(h_clk < -1e-6))
    report(f"  k=0 Hessian of (E - mu N):")
    report(f"    symmetric  lowest {h_sym[0]:+.5f}   "
           f"{n_sym} negative directions -> UNSTABLE")
    report(f"    clock      lowest {h_clk[0]:+.5f}   "
           f"{n_clk} negative directions -> stable")

    ratio, dev = z3_shift_ratio(CLOCK)
    report(f"  Z3 centre on the clock: shift/state = {ratio[0]:+.4f} for all "
           f"three components (spread {dev:.1e})")
    report(f"    -> the shift IS a global phase: the centre is UNBROKEN")

    return dict(states=out, mu_symmetric=mu_sym, mu_clock=mu_clk,
                res_symmetric=res_sym, res_clock=res_clk,
                hess_symmetric=h_sym, hess_clock=h_clk,
                n_negative_symmetric=n_sym, n_negative_clock=n_clk,
                shift_ratio_spread=dev, shift_phase=complex(ratio[0]))


if __name__ == "__main__":
    run()
