#!/usr/bin/env python3
"""
E8 -- Modular end-to-end derivation for the three-paper series:
    1. "The Innocent Lepton" (lepton branch)
    2. "One Substrate, Three Generations" (quark branch)
    3. "The Echo of Standing Waves" (gravity branch)
       Zenodo. https://doi.org/10.5281/zenodo.20144381

Reproduces every numerical prediction from the conformal embedding
E_8(1) > G_2(1) x F_4(1), using only the Planck mass M_Pl
(dimensional, sets the unit system) as input.  All dimensionless
couplings are derived; zero free parameters.

Layer 10 (gravity): Sakharov induced gravity from the (7,26) bridge
sector, deriving G_N from the heat-kernel sum.

References (in dependency order):
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
        Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
        Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
        Zenodo. https://doi.org/10.5281/zenodo.20144381

Usage:
    python -m E8              # full layer-by-layer derivation
    python -m E8 --summary    # final scorecard + emergence tree only

Each derivation layer lives in its own module.  This runner calls
them in dependency order and collects all results.
"""

import argparse
import contextlib
import io
import sys

from . import algebra
from . import scale
from . import leptons
from . import quarks
from . import wzw
from . import ckm
from . import pmns
from . import neutrinos
from . import alpha_s
from . import higgs
from . import gravity
from . import baryogenesis
from . import proofs
from . import scorecard


def _silent_if(quiet: bool):
    """Return a context manager that mutes stdout when quiet is True."""
    if quiet:
        return contextlib.redirect_stdout(io.StringIO())
    return contextlib.nullcontext()


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="python -m E8",
        description="End-to-end E_8(1) > G_2(1) x F_4(1) derivation of the Standard Model.",
    )
    parser.add_argument(
        "-s", "--summary",
        action="store_true",
        help="suppress per-layer output; print only the final scorecard and emergence tree",
    )
    args = parser.parse_args(argv)
    quiet = args.summary

    with _silent_if(quiet):
        # Layer 0 -- Algebraic root: E_8 > G_2 x F_4
        alg = algebra.derive()

        # Layer 1 -- Scale emergence: v_EW from M_Pl
        scl = scale.derive(alg)

        # Layer 2 -- Lepton masses (G_2 sector)
        lep = leptons.derive(scl)

        # Layer 3 -- Quark masses (F_4 sector)
        qrk = quarks.derive(alg, scl, lep)

        # Layer 4 -- SU(3)_3 WZW data (CKM/PMNS engine)
        wzw_data = wzw.derive()

        # Layer 5 -- CKM matrix (four WZW structures)
        ckm_data = ckm.derive(wzw_data)

        # Layer 6 -- PMNS mixing (conjugation invariant)
        pmns_data = pmns.derive(wzw_data, lep, ckm_data)

        # Layer 6b -- Neutrino masses (F_4 neutral sector)
        nu = neutrinos.derive(alg)

        # Layer 7 -- Strong coupling (embedding index chain)
        as_data = alpha_s.derive(alg, scl, qrk)

        # Layer 8 -- Higgs mass from F_4(1) fusion + bridge threshold
        hig = higgs.derive(scl, qrk, as_data, lep, alg)

        # Layer 10 -- Gravity: Sakharov induced gravity from the bridge sector
        grav = gravity.derive(alg, scl)

        # Layer 11 -- Baryogenesis: eta_B from emergence
        baryo = baryogenesis.derive(alg, pmns_data)

        # Layer 12 -- Proof certificates
        proof_data = proofs.derive(alg, scl, lep, qrk, wzw_data, ckm_data, hig, pmns_data, as_data, grav)

    # Global scorecard -- always printed
    scorecard.derive(lep, qrk, ckm_data, hig, pmns_data, as_data, proof_data, grav, nu, baryo)


if __name__ == "__main__":
    main(sys.argv[1:])
