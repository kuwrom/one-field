"""The substrate's own record area: vacuum entanglement of the Bogoliubov condensate.

THE QUESTION. The gravity lane's coefficient is an area density of
records, 1/(4G). Every regulator used for it is a declared choice
(sakharov_normalisation.py FACT 7). A condensate carries a physical
regulator: its Bogoliubov vacuum becomes a product state above the
healing scale, so its entanglement across a surface needs no cutoff. What
record area does the substrate fix, and what does it give for t0?

PREMISES. One complex field of mass 1 (hbar = 1) at chemical potential
mu, its Bogoliubov vacuum in atomic field quadratures,
X_k = eps/(2w), P_k = w/(2 eps), eps = k^2/2, w = sqrt(eps(eps + 2mu)),
both tending to 1/2 at large k. The healing length is the repo's
screening length, xi^2 = 1/(4 mu), the static response of the q = 0
dispersion w^2 = mu p^2 + p^4/4 (sakharov_normalisation.py FACT 5);
xi0 = l_Pl/d10 (root.py C_EFF). A minimal scalar's area law under the
proper time cutoff is S/A = 1/(48 pi t0) (the induced gravity and
entanglement coefficients coincide for minimal scalars, Jacobson
gr-qc/9404039, Solodukhin 1104.3712 eqs. 318-320).

STATEMENTS (the numbers below check each one).
  R1  The half space entropy per area converges with no cutoff:
      S/A = 0.02108, 0.02217, 0.02253 mu at lattice spacing 0.2, 0.1,
      0.05, independent of the box, tending to about 0.0227 mu.
  R2  In healing units kappa = (S/A) xi^2 = 0.0057 per phonon branch.
  R3  Read as a minimal scalar's area law, t0 = xi^2/(48 pi kappa)
      = 1.17 xi^2, which at xi0 = l_Pl/2 is 0.29 l_Pl^2 against the
      4.297 l_Pl^2 the counted ledger needs (4.112 on the engine's one
      list, a factor 14.1). It is short by a factor near 15, and by 28
      and 56 for the other healing conventions. Closing would
      need xi0 = 1.92 l_Pl.

  R4  Read the other way, the closure fixes the resolution. With the
      substrate's regulator under every ledger state, 1/(4G) = 6 Sigma
      kappa / xi^2 for Sigma in units of one scalar = 1/6, so
      xi0/l_Pl = sqrt(24 kappa Sigma) = 1.918 for the counted ledger's
      covariant contact term, 2.48 and 2.93 for the contact coefficients 0 and
      -1, against the identification l_Pl/d10 = 0.5.
  R5  As a vacuum regulator the closure fails gate C. A condensate
      filling space shifts the tensor speed by order one unless
      xi > 2.7e6 l_Pl (GW170817, carrier_sakharov.py), so 1.918 l_Pl
      misses by six orders of magnitude. The vacuum of the lane is the
      empty field, whose Fock vacuum carries no record. The Bogoliubov record is a
      lump interior quantity, and a lump interior does not set a
      universal G. The closure therefore does not constrain G's
      normalization, which the boundary identification carries
      (root.py: t0 = Sigma/(2 pi) l_Pl^2 in proper time form). The bound is
      computed on the gas sign uniform clock, this probe's own medium.
      On the committed sign a liquid filling space is excluded by
      energy, the empty field being the global minimum
      (lump_stability.py S3).

RESULT. The substrate fixes a finite record area. It is a lump interior
record that acts as a regulator only for a Planckian carrier
(closed_event_scale.py). With the Planck identification xi0 = l_Pl/2
it gives an induced stiffness about 14 times too strong.
The identification and this reading of the record do not both hold.
Taken as the closure, the substrate's healing length is 1.9 to 2.9
Planck lengths, the spread being the vector contact convention. The
premise is that every ledger state is an excitation of the substrate
(the gauge as substrate program). This record area is then the only
quantity that depends on xi0. On the committed clock medium the common
sector is this same problem with mu -> c_s^2, so the value holds there
too. In the range convention xi = 1/(2 c_s) it is at least half the
carrier's Compton length (closed_event_scale.py S1). The gapped
relative sector adds a 15% record that no ledger state reads
(clock_record.py).

Usage: python3 tests/probes/record_area.py     (about a minute)
"""

import math

import numpy as np


def chain_entropy(kperp, a, N, mu=1.0):
    k = 2 * np.pi * np.fft.fftfreq(N, d=a)
    eps = (k**2 + kperp**2) / 2
    w = np.sqrt(eps * (eps + 2 * mu))
    with np.errstate(divide="ignore", invalid="ignore"):
        Xk = np.where(eps > 0, eps / (2 * w), 0.0)
        Pk = np.where(eps > 0, w / (2 * eps), 0.0)
    Xr, Pr = np.real(np.fft.ifft(Xk)), np.real(np.fft.ifft(Pk))
    M = N // 2
    idx = np.abs(np.subtract.outer(np.arange(M), np.arange(M)))
    nu = np.sqrt(np.clip(np.real(np.linalg.eigvals(Xr[idx] @ Pr[idx])), 0.25, None))
    nu = nu[nu > 0.5 + 1e-12]
    return np.sum((nu + 0.5) * np.log(nu + 0.5) - (nu - 0.5) * np.log(nu - 0.5)) / 2


def area_density(a, N, kmax=12.0, nk=60):
    kp = np.linspace(1e-3, kmax, nk)
    s = np.array([chain_entropy(x, a, N) for x in kp])
    return float(np.trapezoid(kp * s, kp) / (2 * np.pi))


def run(report=print, fine=False):
    grid = ((0.2, 400), (0.1, 800)) + (((0.05, 1600),) if fine else ())
    vals = {a: area_density(a, N) for a, N in grid}
    a1, a2 = 0.2, 0.1
    limit = vals[a2] + (vals[a2] - vals[a1]) / 2          # first order Richardson in a
    kappa = limit / 4.0                                    # xi^2 = 1/(4 mu)
    t0_xi2 = 1 / (48 * math.pi * kappa)
    t0_lpl2 = t0_xi2 / 4.0                                 # xi0 = l_Pl/2
    need = 4.297                                           # the counted ledger's t0 (FACT 7's reference)
    assert abs(vals[0.2] - 0.02108) < 2e-4 and abs(vals[0.1] - 0.02217) < 2e-4
    assert 0.0222 < limit < 0.0232 and 10 < need / t0_lpl2 < 18
    report("R1 S/A = " + ", ".join(f"{v:.5f} (a = {a})" for a, v in vals.items())
           + f" mu. Limit ~ {limit:.4f} mu")
    report(f"R2 kappa = {kappa:.5f} per xi^2 (xi^2 = 1/(4 mu))")
    report(f"R3 t0 = {t0_xi2:.3f} xi^2 = {t0_lpl2:.3f} l_Pl^2 at xi0 = l_Pl/2 against {need} "
           f"(short by {need / t0_lpl2:.1f}). Closing needs xi0 = {math.sqrt(need / t0_xi2):.2f} l_Pl")
    bridge = 182 / 6                                       # the amount, minimally coupled
    xi_closure = {c: math.sqrt(24 * kappa * (7 / 6 + 45 / 6 + 36 * (1 / 6 - c / 2) + bridge))
                  for c in (1.0, 0.0, -1.0)}
    assert abs(xi_closure[1.0] - 1.918) < 0.02
    report("R4 closure: xi0/l_Pl = sqrt(24 kappa Sigma) = " + ", ".join(
        f"{v:.3f} (contact {k:+.0f})" for k, v in xi_closure.items()) + ". Identification 0.5")
    return dict(xi_closure=xi_closure, S_over_A=vals, limit=limit, kappa=kappa, t0_xi2=t0_xi2, t0_lpl2=t0_lpl2,
                shortfall=need / t0_lpl2)


if __name__ == "__main__":
    run(fine=True)
