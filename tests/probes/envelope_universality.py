"""
envelope_universality.py: the NLS form derived as the universal
envelope normal form of weak local self-interference, up to
locality.

THE ARGUMENT.  Taking the substrate to be an NLS medium looks like
an assumption.  The
multiple-scales theorem behind this probe: for ANY local, dispersive,
weakly nonlinear wave field, the slow modulation of a carrier obeys
the nonlinear Schrodinger equation at leading order, with
coefficients read off the microdynamics (kinetic = half the
dispersion curvature, cubic = the leading nonlinear response).  The
NLS form is therefore not a modeling choice about the substrate; it
is the normal form of self-interference at lowest order, whatever the
microdynamics.  The Z3 coupling structure passes to the envelope
linearly (the cubic coefficient matrix projects channel by channel),
so the Z3-NLS of tests/substrate/ is the normal form of any Z3-coupled
local medium.  What remains assumed is locality alone.

THE COMPUTATION.  The theorem is run on a microdynamics that is NOT
NLS: the phi^4 Klein-Gordon field
u_tt = u_xx - u - beta u^3 (second order in time, relativistic
dispersion).  A modulated carrier is evolved under the full KG
dynamics; its envelope is extracted by exact demodulation
(A = (u + i u_t/w0) e^{-i theta} / 2, low-passed); and the extracted
envelope is compared against the NLS evolution of the initial
envelope with coefficients DERIVED from the KG data, no fitting:
group velocity vg = k0/w0, kinetic w'' = 1/w0^3, cubic
lambda = 3 beta / (2 w0).

MEASURED (deterministic, frozen):
  eps = 0.12: NLS envelope error 0.040; linear envelope error 0.66
  eps = 0.06: NLS envelope error 0.010; linear envelope error 0.17
  convergence: the NLS error scales as eps^2 (ratio 0.24 at eps
  halving), the normal-form signature; the cubic term is load-bearing
  (NLS beats the linear envelope by a factor > 15).

Usage: python3 envelope_universality.py     (~20 s)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401

L, N = 200.0, 1024
X = np.linspace(0, L, N, endpoint=False)
DX = L / N
K = 2 * np.pi * np.fft.fftfreq(N, d=DX)
K0 = 1.0
W0 = np.sqrt(1 + K0 ** 2)
VG = K0 / W0                      # group velocity
WPP = 1.0 / W0 ** 3               # dispersion curvature
BETA = 1.0
LAM = 3 * BETA / (2 * W0)         # cubic coefficient, derived not fitted
T = 60.0


def _kg_run(eps, t_end, dt=0.01):
    """The full (non-NLS) microdynamics: phi^4 Klein-Gordon."""
    a0 = eps * np.exp(-((X - L / 2) / 18.0) ** 2)
    th = K0 * X
    u = 2 * a0 * np.cos(th)
    ut = 2 * a0 * W0 * np.sin(th)
    lap = lambda f: np.fft.ifft(-(K ** 2) * np.fft.fft(f)).real
    acc = lap(u) - u - BETA * u ** 3
    for _ in range(int(t_end / dt)):
        ut2 = ut + 0.5 * dt * acc
        u = u + dt * ut2
        acc = lap(u) - u - BETA * u ** 3
        ut = ut2 + 0.5 * dt * acc
    return u, ut


def _extract(u, ut, t):
    """Exact demodulation: A = (u + i u_t / w0) e^{-i theta} / 2,
    low-passed to the envelope band."""
    th = K0 * X - W0 * t
    a = 0.5 * (u + 1j * ut / W0) * np.exp(-1j * th)
    ak = np.fft.fft(a)
    ak[np.abs(K) > 0.5 * K0] = 0
    return np.fft.ifft(ak)


def _envelope_run(a, t_end, lam, dt=0.005):
    """The derived envelope equation in the lab frame:
    i A_t = -i vg A_x - (w''/2) A_xx + lam |A|^2 A."""
    kin = np.exp(-1j * dt * (VG * K + 0.5 * WPP * K ** 2))
    for _ in range(int(t_end / dt)):
        a = np.fft.ifft(kin * np.fft.fft(a))
        a = a * np.exp(-1j * dt * lam * np.abs(a) ** 2)
    return a


def run(report=print):
    report("ENVELOPE UNIVERSALITY: the NLS normal form, run on a "
           "non-NLS microdynamics")
    report("=" * 70)
    report(f"  microdynamics: u_tt = u_xx - u - u^3 (KG); carrier "
           f"k0 = {K0}, w0 = {W0:.4f}")
    report(f"  derived envelope coefficients: vg = {VG:.4f}, "
           f"w''/2 = {WPP / 2:.4f}, lambda = {LAM:.4f} (no fitting)")

    errs = {}
    for eps in (0.12, 0.06):
        u, ut = _kg_run(eps, T)
        a_true = _extract(u, ut, T)
        u0, ut0 = _kg_run(eps, 0.0)
        a_0 = _extract(u0, ut0, 0.0)
        a_nls = _envelope_run(a_0.copy(), T, LAM)
        a_lin = _envelope_run(a_0.copy(), T, 0.0)
        e_nls = float(np.linalg.norm(a_true - a_nls)
                      / np.linalg.norm(a_true))
        e_lin = float(np.linalg.norm(a_true - a_lin)
                      / np.linalg.norm(a_true))
        errs[eps] = (e_nls, e_lin)
        report(f"  eps = {eps}: NLS envelope err {e_nls:.4f}, "
               f"linear envelope err {e_lin:.4f} "
               f"(cubic term load-bearing x{e_lin / e_nls:.0f})")

    # frozen assertions (measured 0.040/0.010; linear 0.66/0.17)
    assert errs[0.12][0] < 0.08, "NLS must track the KG envelope"
    assert errs[0.06][0] < 0.025
    ratio = errs[0.06][0] / errs[0.12][0]
    assert ratio < 0.45, \
        f"normal-form eps^2 convergence broken (ratio {ratio:.2f})"
    assert errs[0.12][1] / errs[0.12][0] > 5.0, \
        "the cubic term must be load-bearing vs the linear envelope"
    assert errs[0.06][1] / errs[0.06][0] > 5.0

    report("-" * 70)
    report("  VERDICT: the envelope of a non-NLS microdynamics obeys "
           "the DERIVED NLS to")
    report(f"  normal-form accuracy (error ~ eps^2: ratio "
           f"{ratio:.2f} at eps halving).  The NLS form is")
    report("  the universal envelope of weak local self-interference; "
           "the residue is locality.")
    return dict(errors=errs, ratio=ratio,
                coefficients=dict(vg=VG, wpp=WPP, lam=LAM))


if __name__ == "__main__":
    run()
