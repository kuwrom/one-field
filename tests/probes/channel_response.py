"""The clock medium's static backreaction to two held biases, channel by channel.

THE QUESTION. substrate_coupling.py defines the substrate's coupling
as the force the medium carries between two of its quanta, in units
of hbar c_s, and finds 1/(8 pi rho xi^3) for a one component medium
with any equation of state. Nothing is exchanged and nothing is
dressed. The force is the medium's static backreaction to two held
biases. It is the inverse of the medium's stiffness at second order,
and the e^(-r/xi)/r shape is the Green's function of that stiffness,
with xi its screening length. The committed medium, the lump interior,
is the three component clock psi_k = sqrt(rho0) w^k of
lump_in_empty_field.py (g0 < 0 < g1, sextic (kappa/3) R^3, zero
pressure). What does the same definition give there, in the common
channel and in the two relative channels? And in which convention
should the depth it implies be quoted?

THE COMPUTATION. Six real fields at the clock point. The static
Hessian of E - mu N is H(k^2) = k^2 I + H0, H0 = W - 2 mu I, with W
the Hessian of the interaction density w. Every eigenvector of H0 is
an eigenvector of H(k^2), so the backreaction energy between two static
biases S and S' is an exact sum of screened 1/r terms,
    V(r) = - sum_i (S.v_i)(S'.v_i) e^{-sqrt(lam_i) r} / (4 pi r),
alpha_i = (S.v_i)(S'.v_i)/(4 pi c_s) in units of hbar c_s. A bias of
species k is a held quantum coupled through the same w by its
density, Omega_k = dw/drho_k at fixed phases = g0 R + g1 Re(psibar_
{k+1} psi_{k+2}) + kappa R^2. The Z3 Fourier transform over (k, k')
separates the channels q = 0, 1, 2 exactly. Hessian and gradients are
7-point stencils, exact for the degree-6 polynomial w.

STATEMENTS (the numbers below check each one).
  X1 (common channel; two readings of "a quantum"). The q = 0 channel
     is one screened term at 2 c_s, i.e. the range xi = 1/(2 c_s). Its
     strength depends on what the quantum couples by. Coupled by the
     chemical potential d mu/dR = c_s^2/R0 it is exactly 1/(8 pi R0 xi^3),
     substrate_coupling.py's law with R0 the total density. Coupled by
     its own Hartree potential, g0 - g1/6 + 2 kappa R0, it is
     (1 + g1 R0/(6 c_s^2))^2 = 1.41787 times that at g1/g0 = 1/sqrt 2.
     The two coincide for one component. The interspecies exchange
     term g1 rho_{k+1} Re(psibar_{k+2} psi_k) separates them. The law
     "for any equation of state" is a one component statement.
  X2 (relative channels). q = 1 and q = 2 are each two screened terms
     at the static relative stiffnesses k^2 = (5 -+ sqrt 13) r/2, r = g1 rho0, the
     zero frequency roots of gravity_tensor.clock_frequencies: ranges
     3.878 xi and 1.561 xi. Summed over both channels the same species
     relative attraction is 0.12830 of the common Hartree channel.
     Between different species it is repulsive and half as strong.
     The contact strengths are g1 (2/3) per relative channel.
  X3 (depth conventions). On the clock mu0 = c_s^2/2 (the sextic),
     so the FACT 5 length xi^2 = 1/(4|mu0|) is 2 xi_range^2, and the
     O4 depth sqrt(kappa)/(2 sqrt3 a^2) (saturation_origin.py) is
     2^(3/2) = 2.828 times the range convention depth R0/(8 c_s^3).
     Depth 3 in the range convention is depth 8.485 in O4's and
     kappa = 125.9. Kappa = 15.7 is depth 3 in O4's convention, which
     is depth 1.061 in the range convention.
  X4 (what alpha_G2 = 1/(24 pi) fixes). In the range convention the
     depth is 3 under the chemical potential reading and 4.254 under
     the Hartree reading: the normalization is a binary choice between
     two readings of "a quantum of the medium", not a continuum.

RESULT. The relative channel backreaction is computed. It is a pure
function of g1/g0, and no engine number uses it. The depth implied by
alpha_G2 carries one discrete reading (3 or 4.25) and must be quoted
in the range convention. The O4 coefficient at depth 3 is 126.

Usage: python3 tests/probes/channel_response.py     (about a second)
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from nls_soliton import G1_OVER_G0  # noqa: E402

G0, KAPPA = -1.0, 1.0
W3 = np.exp(2j * np.pi * np.arange(3) / 3)
_STENCIL1 = np.array([-1, 9, -45, 0, 45, -9, 1]) / 60.0          # f'(0), exact to degree 6
_STENCIL2 = np.array([2, -27, 270, -490, 270, -27, 2]) / 180.0   # f''(0), exact to degree 7
_H = 0.05


def medium(g1):
    """Clock medium at zero pressure for g0 = -1, sextic kappa = 1."""
    a = G0 / 2 - g1 / 6
    R0 = 3 * abs(a) / (2 * KAPPA)
    rho0 = R0 / 3
    mu = 2 * a * R0 + KAPPA * R0**2
    cs2 = R0 * (2 * a + 2 * KAPPA * R0)
    point = np.concatenate([np.sqrt(rho0) * W3.real, np.sqrt(rho0) * W3.imag])
    return dict(g1=g1, a=a, R0=R0, rho0=rho0, mu=mu, cs2=cs2, cs=math.sqrt(cs2), point=point)


def _fields(z):
    psi = z[:3] + 1j * z[3:]
    rho = np.abs(psi)**2
    O = np.array([(np.conj(psi[(k + 1) % 3]) * psi[(k + 2) % 3]).real for k in range(3)])
    return psi, rho, O


def w_density(z, g1):
    psi, rho, O = _fields(z)
    R = rho.sum()
    return 0.5 * G0 * R * R + g1 * float(rho @ O) + KAPPA * R**3 / 3


def omega(z, g1, k):
    """Potential felt by a held quantum of species k."""
    psi, rho, O = _fields(z)
    R = rho.sum()
    return G0 * R + g1 * O[k] + KAPPA * R * R


def _d1(f, z, u):
    return sum(c * f(z + j * _H * u) for j, c in zip(range(-3, 4), _STENCIL1)) / _H


def _d2(f, z, u):
    return sum(c * f(z + j * _H * u) for j, c in zip(range(-3, 4), _STENCIL2)) / _H**2


def hessian(f, z):
    n = z.size
    E = np.eye(n)
    H = np.zeros((n, n))
    for i in range(n):
        H[i, i] = _d2(f, z, E[i])
    for i in range(n):
        for j in range(i + 1, n):
            H[i, j] = H[j, i] = (_d2(f, z, E[i] + E[j]) - _d2(f, z, E[i] - E[j])) / 4
    return H


def gradient(f, z):
    E = np.eye(z.size)
    return np.array([_d1(f, z, E[i]) for i in range(z.size)])


def static_problem(g1):
    m = medium(g1)
    z = m["point"]
    Wh = hessian(lambda q: w_density(q, g1), z)
    H0 = Wh - 2 * m["mu"] * np.eye(6)
    lam, V = np.linalg.eigh(H0)
    S = np.array([gradient(lambda q, k=k: omega(q, g1, k), z) for k in range(3)])
    grad_R = gradient(lambda q: float(np.sum(q * q)), z)
    m.update(H0=H0, lam=lam, V=V, S=S, grad_R=grad_R)
    return m


def channel_amplitudes(m):
    """A[q, i] = |S_q . v_i|^2 with S_q = sum_s w^{-qs} S_s / sqrt 3.
    V_kk'(r) = (1/3) sum_q w^{q(k - k')} sum_i (-A[q, i]) e^{-sqrt(lam_i) r}/(4 pi r)."""
    A = np.zeros((3, 6))
    for q in range(3):
        Sq = sum(np.exp(-2j * np.pi * q * s / 3) * m["S"][s] for s in range(3)) / math.sqrt(3)
        A[q] = np.abs(Sq @ m["V"])**2
    return A


def run(report=print):
    out = {}
    g1 = G1_OVER_G0
    m = static_problem(g1)
    cs, R0, lam = m["cs"], m["R0"], m["lam"]
    xi = 1 / (2 * cs)
    r = g1 * m["rho0"]
    report("CHANNEL BACKREACTION of the clock medium (g0 = -1, g1/g0 = 1/sqrt 2, kappa = 1)")
    report(f"  R0 = {R0:.5f}, mu0 = {m['mu']:.5f} = c_s^2/2 ({m['cs2'] / 2:.5f}), c_s = {cs:.5f}, xi = 1/(2 c_s) = {xi:.5f}")

    # structure of the static Hessian: Goldstone, common density at 4 c_s^2, two relative pairs
    expected = np.sort(np.array([0.0, 4 * m["cs2"], (5 - math.sqrt(13)) * r / 2, (5 - math.sqrt(13)) * r / 2,
                                 (5 + math.sqrt(13)) * r / 2, (5 + math.sqrt(13)) * r / 2]))
    out["eigen_error"] = float(np.max(abs(lam - expected)))
    assert out["eigen_error"] < 1e-8, (lam, expected)
    # the radial direction is the 4 c_s^2 eigenvector
    vR = m["grad_R"] / np.linalg.norm(m["grad_R"])
    out["radial_residual"] = float(np.linalg.norm(m["H0"] @ vR - 4 * m["cs2"] * vR))
    assert out["radial_residual"] < 1e-8
    report(f"  H0 eigenvalues k^2: {np.round(lam, 6)} = 0, (5 -+ sqrt13) r/2 (x2 each), 4 c_s^2 "
           f"(error {out['eigen_error']:.1e}). Grad R is the 4 c_s^2 eigenvector (residual {out['radial_residual']:.1e})")

    A = channel_amplitudes(m)
    alpha = A / (4 * math.pi * cs) / 3                     # per pair of quanta, units hbar c_s
    # X1: common channel
    alpha_law = cs**3 / (math.pi * R0)                      # 1/(8 pi R0 xi^3)
    assert abs(alpha_law - 1 / (8 * math.pi * R0 * xi**3)) < 1e-14
    g_hartree = G0 - g1 / 6 + 2 * KAPPA * R0
    g_mu = m["cs2"] / R0
    assert abs(g_mu - (G0 - g1 / 3 + 2 * KAPPA * R0)) < 1e-12
    ratio_sq = (g_hartree / g_mu)**2
    common = alpha[0]
    assert np.sum(common > 1e-12) == 1 and abs(lam[np.argmax(common)] - 4 * m["cs2"]) < 1e-8
    out["alpha_common_hartree"] = float(common.max())
    out["alpha_common_mu"] = alpha_law
    out["hartree_over_mu"] = float(common.max() / alpha_law)
    assert abs(out["hartree_over_mu"] - ratio_sq) < 1e-8
    assert abs(ratio_sq - (1 + g1 * R0 / (6 * m["cs2"]))**2) < 1e-12
    report(f"X1 common channel: one screened term at 2 c_s. mu coupled alpha = 1/(8 pi R0 xi^3) = {alpha_law:.6f}. "
           f"Hartree coupled alpha = {common.max():.6f} = {ratio_sq:.5f} x that "
           f"(couplings {g_mu:.5f} vs {g_hartree:.5f}).")

    # X2: relative channels
    rel = alpha[1]
    assert np.allclose(alpha[1], alpha[2])
    keep = rel > 1e-12
    masses = np.sqrt(lam[keep])
    ranges = (1 / masses) / xi
    out["relative_ranges_xi"] = sorted(set(np.round(ranges, 3)))
    same_species_rel = 2 * rel.sum()                        # (1/3)(V_1 + V_2) with V_q = 3 alpha_q
    cross_species_rel = -rel.sum()
    out["relative_over_common"] = float(same_species_rel / common.max())
    out["contact_relative"] = float(np.sum(A[1] / np.where(lam > 1e-12, lam, np.inf)) / 3)
    assert abs(out["contact_relative"] - 2 * g1 / 3) < 1e-8
    assert abs(out["relative_over_common"] - 0.12830) < 2e-5
    assert abs(sorted(set(np.round(ranges, 3)))[0] - 1.561) < 2e-3 and abs(sorted(set(np.round(ranges, 3)))[1] - 3.878) < 2e-3
    report(f"X2 relative channels q = 1, 2: screened terms at ranges {np.round(sorted(set(np.round(ranges, 3))), 3)} xi, "
           f"amplitudes {np.round(rel[keep], 6)} each. Same species sum {same_species_rel:.6f} = "
           f"{out['relative_over_common']:.5f} of the common Hartree channel. Cross species {cross_species_rel:.6f}. "
           f"Contact {out['contact_relative']:.5f} = 2 g1/3.")

    # X3: conventions
    depth_range = R0 * xi**3
    xi_o4 = 1 / math.sqrt(2 * abs(m["a"]) * R0)
    depth_o4 = R0 * xi_o4**3
    out["depth_range"], out["depth_o4"] = depth_range, depth_o4
    assert abs(depth_o4 / depth_range - 2**1.5) < 1e-12
    assert abs(depth_o4 - math.sqrt(KAPPA) / (2 * math.sqrt(3) * m["a"]**2)) < 1e-12
    kappa_depth3 = (3 * 2**1.5 * 2 * math.sqrt(3) * m["a"]**2)**2       # depth 3 in the range convention
    out["kappa_at_range_depth_3"] = kappa_depth3
    out["o4_depth_at_kappa_15_7"] = 3 / 2**1.5
    assert abs(kappa_depth3 - 125.9) < 0.1
    report(f"X3 conventions: depth {depth_range:.5f} (range, xi = 1/2c_s) = {depth_o4:.5f} (O4, xi^2 = 1/(4|mu0|)) / 2^(3/2). "
           f"Range depth 3 = O4 depth {3 * 2**1.5:.3f}, kappa = {kappa_depth3:.1f}. Kappa = 15.7 is range depth {3 / 2**1.5:.3f}.")

    # X4: the depth alpha_G2 fixes under the two readings
    alpha_G2 = 1 / (24 * math.pi)
    out["depth_mu"] = 1 / (8 * math.pi * alpha_G2)
    out["depth_hartree"] = out["depth_mu"] * ratio_sq
    assert abs(out["depth_mu"] - 3) < 1e-12 and abs(out["depth_hartree"] - 4.254) < 2e-3
    report(f"X4 alpha_G2 = 1/(24 pi): range convention depth {out['depth_mu']:.3f} (mu coupled) or "
           f"{out['depth_hartree']:.3f} (Hartree coupled). A binary reading, not a continuum.")

    # the ratio across g1/g0, for the record
    scan = {}
    for g in (0.3, 0.5, 0.9, 1.2):
        mm = static_problem(g)
        AA = channel_amplitudes(mm)
        al = AA / (4 * math.pi * mm["cs"]) / 3
        scan[g] = dict(hartree_over_mu=float(al[0].max() / (mm["cs"]**3 / (math.pi * mm["R0"]))),
                       relative_over_common=float(2 * al[1].sum() / al[0].max()))
    out["scan"] = scan
    report("    scan g1/g0 -> (Hartree/mu, relative/common): " +
           ", ".join(f"{g}: ({v['hartree_over_mu']:.4f}, {v['relative_over_common']:.4f})" for g, v in scan.items()))
    report("RESULT: the relative backreaction is computed, and no engine number uses it. The depth must be quoted in the "
           "range convention and carries one discrete reading (3 or 4.25).")
    return out


if __name__ == "__main__":
    run()
