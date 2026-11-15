"""The committed medium's own record area, and the smallest closed event.

THE QUESTION. record_area.py computed the substrate's record
(entanglement entropy) per area on a one component gas sign medium and
read xi0/l_Pl = 1.918 from it (R4). The committed medium is the three
component clock at zero pressure (lump_in_empty_field.py). This probe
asks three things. What record does that medium leave per healing
area, sector by sector? How many quanta and how much area does the
smallest closed event (the minimal lump) hold? Does any reading of the
record match the identification xi0 = l_Pl/d10?

The record is a lump interior quantity throughout. The vacuum of the
lane is the empty field, and its Fock vacuum carries no record. A
medium filling space would need xi > 2.7e6 l_Pl to keep c_T = c on the
gas sign (carrier_sakharov.py, record_area.py R5). On the committed
sign the empty field is the global minimum (lump_stability.py S3). So
nothing below normalizes G. That normalization is the boundary
identification of root.py. What is fixed here is what a lump's
interior records.

THE COMPUTATION. Per momentum the quadratic Hamiltonian is H2 = (1/2)
Z^T K Z, Z = (X, P) canonical, K(k) = H(k^2)/2 with H the static
Hessian of channel_response.py. The ground state covariance is
    gamma = (1/2) K^{-1/2} (K^{1/2} J^T K J K^{1/2})^{1/2} K^{-1/2}.
The corotating Z3 sectors q = 0 and q = {1, 2} decouple at every k
and the transform is site local, so the half space entropy is a sum
over sectors. The half space entropy follows record_area.py: a chain
along x with lattice spacing a and N sites, k_perp integrated, the
symplectic eigenvalues of the half chain covariance; a periodic
half chain has two cuts, so the per boundary record is half of it.

STATEMENTS (the numbers below check each one).
  K1 (sectors). The corotating transform block diagonalises H to
     1e-13: q = 0 is diag(0, 4 c_s^2), the common density and phase;
     q = {1, 2} has the static stiffnesses (5 -+ sqrt13) r/2 twice.
  K2 (common sector = record_area). The q = 0 covariance is X_k =
     eps/(2w), P_k = w/(2 eps) with mu -> c_s^2, and its record per
     area at a = 0.4 xi and 0.2 xi is record_area's 0.02108 and
     0.02217 in units of c_s^2 (there, of mu). kappa_0 = 0.0057 per
     xi^2 with xi = 1/(2 c_s), depth free, so R4's closure value
     xi/l_Pl = 1.918 (range convention) is the same number on the
     committed medium. It is the screening length at which the lump
     interior's record per area would equal the Einstein coefficient.
     By R5 it is not a vacuum regulator and does not fix G.
  K3 (relative sector). The gapped sector q = {1, 2} leaves
     kappa_rel = 0.00085 per xi^2 per boundary, 0.15 of kappa_0,
     converged in the lattice and in k_perp. It is gapped, so under
     the per state reading of R4 it regulates no massless state and
     changes no closure. It is the record held by the clock's lock
     itself.
  K4 (the smallest closed event). The minimal lump of
     lump_in_empty_field.py (mu/mu0 = 0.1362, N = 43.87 in the
     probe's units) has half density radius 2.755 xi and area
     95.4 xi^2 there, and holds N_min = 164.1 x depth quanta: 492 at
     the range convention depth 3 and 698 at 4.25. The values of
     xi/l_Pl are: the per state Bogoliubov record from the interior,
     1.918. One record per quantum, sqrt(4 N_min xi^2 / A) = 4.54
     (depth 3) or 5.41. The identification l_Pl/d10 = 0.5. No reading
     of the interior's record matches the identification. None of the
     interior readings is a candidate for G's normalization. Each is a
     property of the smallest closed event, not of the empty field,
     and each fixes a Planckian carrier mass that no gate reads
     (closed_event_scale.py).

RESULT. The committed medium's record is the one component record in
its common sector plus a 15% gapped relative record. Both are records
of a lump's interior. The closure value 1.918 l_Pl is the same number
on this medium, and no reading reproduces the identification 0.5.
Neither fixes G. G's normalization is the boundary identification of
root.py (record_area.py R5). The smallest closed event is a lump of
about five hundred quanta. Its record per area is a derived count of
what the lump holds, not of what the vacuum resolves.

Usage: python3 tests/probes/clock_record.py     (about fifteen seconds)
"""

import math
import os
import sys

import numpy as np
from scipy.linalg import sqrtm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import channel_response as C  # noqa: E402
import lump_in_empty_field as L  # noqa: E402

IDX0, IDX12 = [0, 3], [1, 2, 4, 5]


def sector_transform():
    """Orthogonal 6x6 map (x_k, y_k) <- corotating sector fields (Re phi_q, Im phi_q):
    eta_k = w^k sum_q w^{qk} phi_q / sqrt 3."""
    T = np.zeros((6, 6))
    for k in range(3):
        for q in range(3):
            c = np.exp(2j * np.pi * k * (1 + q) / 3) / math.sqrt(3)
            T[k, q] += c.real
            T[k, 3 + q] += -c.imag
            T[3 + k, q] += c.imag
            T[3 + k, 3 + q] += c.real
    return T


def covariance(K):
    n = K.shape[0] // 2
    J = np.block([[np.zeros((n, n)), np.eye(n)], [-np.eye(n), np.zeros((n, n))]])
    Kh = np.real(sqrtm(K))
    Kmh = np.linalg.inv(Kh)
    return 0.5 * Kmh @ np.real(sqrtm(Kh @ J.T @ K @ J @ Kh)) @ Kmh


def chain_entropy(kperp, a, N, Hblock):
    """Half chain entropy per boundary of a sector with Hblock (2n x 2n static Hessian)."""
    n = Hblock.shape[0]
    kx = 2 * np.pi * np.fft.fftfreq(N, d=a)
    gam = np.zeros((N, n, n))
    for i, k in enumerate(kx):
        k2 = max(k * k + kperp * kperp, 1e-12)
        gam[i] = covariance((k2 * np.eye(n) + Hblock) / 2)
    Gr = np.real(np.fft.ifft(gam, axis=0))
    M = N // 2
    sep = np.subtract.outer(np.arange(M), np.arange(M)) % N
    G = Gr[sep].transpose(0, 2, 1, 3).reshape(M * n, M * n)
    h = n // 2
    Jb = np.block([[np.zeros((h, h)), np.eye(h)], [-np.eye(h), np.zeros((h, h))]])
    ev = np.linalg.eigvals(np.kron(np.eye(M), Jb) @ G)
    nu = np.clip(np.sort(np.abs(ev.imag))[::2], 0.5, None)
    nu = nu[nu > 0.5 + 1e-10]
    return float(np.sum((nu + 0.5) * np.log(nu + 0.5) - (nu - 0.5) * np.log(nu - 0.5))) / 2


def area_density(a, N, Hblock, kmax, nk):
    kp = np.linspace(1e-3, kmax, nk)
    s = np.array([chain_entropy(x, a, N, Hblock) for x in kp])
    return float(np.trapezoid(kp * s, kp) / (2 * np.pi)), kp, s


def run(report=print):
    out = {}
    m = C.static_problem(C.G1_OVER_G0)
    cs2, cs = m["cs2"], m["cs"]
    xi = 1 / (2 * cs)
    r = C.G1_OVER_G0 * m["rho0"]
    T = sector_transform()
    assert np.allclose(T @ T.T, np.eye(6), atol=1e-13)
    Hs = T.T @ m["H0"] @ T
    out["decoupling"] = float(np.max(abs(Hs[np.ix_(IDX0, IDX12)])))
    H00, Hrel = Hs[np.ix_(IDX0, IDX0)], Hs[np.ix_(IDX12, IDX12)]
    e0, erel = np.linalg.eigvalsh(H00), np.linalg.eigvalsh(Hrel)
    assert out["decoupling"] < 1e-12 and abs(e0[0]) < 1e-12 and abs(e0[1] - 4 * cs2) < 1e-12
    exp_rel = np.sort([(5 - math.sqrt(13)) * r / 2] * 2 + [(5 + math.sqrt(13)) * r / 2] * 2)
    assert np.max(abs(erel - exp_rel)) < 1e-12
    report(f"K1 sectors: decoupling {out['decoupling']:.1e}; q = 0 stiffnesses {np.round(e0, 6)} = (0, 4 c_s^2); "
           f"q = 1,2 {np.round(erel, 6)} = (5 -+ sqrt13) r/2.")

    # K2: common sector covariance and record
    for k2 in (0.1, 1.0, 5.0):
        g = covariance((k2 * np.eye(2) + H00) / 2)
        eps = k2 / 2
        w = math.sqrt(eps * (eps + 2 * cs2))
        ev = np.sort(np.linalg.eigvalsh(g))
        assert abs(ev[0] * ev[1] - 0.25) < 1e-10 and abs(ev[1] - max(eps / (2 * w), w / (2 * eps))) < 1e-9
    common = {}
    for frac, N in ((0.4, 400), (0.2, 800)):
        v, _, _ = area_density(frac * xi, N, H00, kmax=12 * cs, nk=60)
        common[frac] = v / cs2
    out["common_S_over_A"] = common
    assert abs(common[0.4] - 0.02108) < 3e-4 and abs(common[0.2] - 0.02217) < 3e-4
    limit = common[0.2] + (common[0.2] - common[0.4]) / 2
    out["kappa_0"] = limit * cs2 * xi**2                     # = limit / 4
    assert abs(out["kappa_0"] - 0.0057) < 1e-4
    report(f"K2 common sector: S/A = {common[0.4]:.5f}, {common[0.2]:.5f} c_s^2 at a = 0.4, 0.2 xi "
           f"(record_area 0.02108, 0.02217 mu); kappa_0 = {out['kappa_0']:.5f} per xi^2, depth free, lump interior.")

    # K3: relative sector record
    rel = {}
    for frac, N in ((0.3, 200), (0.15, 400)):
        v, kp, s = area_density(frac * xi, N, Hrel, kmax=30 * cs, nk=60)
        rel[frac] = v * xi**2
        assert (kp * s)[-1] < 2e-2 * max(kp * s)                 # k_perp integrand converged (gapped)
    out["kappa_rel_grid"] = rel
    out["kappa_rel"] = rel[0.15] + (rel[0.15] - rel[0.3]) / 2
    out["rel_over_common"] = out["kappa_rel"] / out["kappa_0"]
    assert abs(out["kappa_rel"] - 0.00085) < 5e-5 and 0.12 < out["rel_over_common"] < 0.18
    report(f"K3 relative sector: kappa_rel = {rel[0.3]:.5f}, {rel[0.15]:.5f} per xi^2 at a = 0.3, 0.15 xi; "
           f"limit {out['kappa_rel']:.5f} = {out['rel_over_common']:.3f} kappa_0. Gapped, regulates no massless state.")

    # K4: the smallest closed event
    d = L.lump(0.13624 * L.MU0)
    assert abs(d["N"] - 43.8675) < 2e-3
    f0 = d["f0"]
    prof = d["profile"]
    from scipy.optimize import brentq
    r_half = brentq(lambda x: prof(x)**2 - 0.5 * f0**2, 1e-3, 40.0)
    xi_L = 1 / (2 * math.sqrt(-L.A * L.R0))
    D_L = L.R0 * xi_L**3
    area = 4 * math.pi * (r_half / xi_L)**2
    per_depth = d["N"] / D_L
    out.update(r_half_xi=r_half / xi_L, area_xi2=area, quanta_per_depth=per_depth)
    assert abs(r_half / xi_L - 2.755) < 3e-3 and abs(area - 95.4) < 0.3 and abs(per_depth - 164.1) < 0.2
    bridge = 182 / 6                                       # the amount, minimally coupled
    sigma = 7 / 6 + 45 / 6 + 36 * (1 / 6 - 1 / 2) + bridge
    xi_state = math.sqrt(24 * out["kappa_0"] * sigma)
    readings = {"per state Bogoliubov": xi_state, "identification l_Pl/d10": 0.5}
    for D in (3.0, 4.254):
        readings[f"one record per quantum, depth {D}"] = math.sqrt(4 * per_depth * D / area)
    out["readings"] = readings
    assert abs(xi_state - 1.918) < 0.03 and abs(readings["one record per quantum, depth 3.0"] - 4.54) < 0.02
    assert min(abs(v - 0.5) for k, v in readings.items() if "identification" not in k) > 1.0
    report(f"K4 smallest closed event: r_half = {r_half / xi_L:.3f} xi, area {area:.1f} xi^2, N_min = {per_depth:.1f} x depth "
           f"(492 at 3, 698 at 4.25). Interior readings of xi/l_Pl: " + ", ".join(f"{k} {v:.3f}" for k, v in readings.items())
           + ". None normalizes G.")
    report("RESULT: lump interior records only (record_area.py R5). The closure value 1.918 l_Pl is the same "
           "number on this medium. No reading reproduces the identification 0.5. The relative sector adds "
           "a 15% record that no state reads. None of it normalizes G.")
    return out


if __name__ == "__main__":
    run()
