"""
casimir.py — Spectral Casimir Effect (QVG prediction)
======================================================
FALSIFIABLE PREDICTION P5.

The spectral Planck law of the vacuum modifies the Casimir pressure:

    P_Cas^QVG = -π²ℏc/(240 d⁴) × (1 + ℏc/(d E₀))

The characteristic length scale is:
    d* = ℏc/E₀ ≈ 62 μm  (for E₀ = 3.2 meV)

Physical regimes:
  d << d* : QVG correction >> 1 (non-perturbative QVG regime)
  d ~  d* : transition (correction ~ 1)
  d >> d* : QFT standard regime (correction → 0)

NOTE ON THE BOOK: The value "0.3% at 65 μm" in Chapter 20 contains
an error. The correct value is δP/P ≈ 95% at d = 65 μm, meaning
the QVG correction dominates the standard Casimir pressure near d*.
This is a stronger prediction: the QVG introduces a NEW REGIME in
Casimir physics at d ~ 62 μm, not just a small correction.
This will be corrected in the revised edition.

Experimental implication:
  - For d < d*: QVG Casimir differs fundamentally from QFT prediction
  - For d > d*: QFT standard limit recovered
  - The transition at d ~ 62 μm is the observable signature

Reference: Chapter 20 (revised edition pending).
"""

import numpy as np

HBAR_C_EV_M = 197.3e-9  # eV·m  (ℏc in natural units)
E0_EV       = 3.2e-3    # eV    (E₀ = 3.2 meV)
D_STAR      = HBAR_C_EV_M / E0_EV  # ≈ 62 μm


def casimir_standard(d_meters):
    """Standard Casimir pressure P = -π²ℏc/(240 d⁴) in Pa."""
    hbar_c_SI = 3.162e-26  # J·m
    return -np.pi**2 * hbar_c_SI / (240 * d_meters**4)


def casimir_qvg(d_meters, E0_eV=E0_EV):
    """
    QVG-corrected Casimir pressure.

    P^QVG = P_std × (1 + ℏc/(d E₀))

    Parameters
    ----------
    d_meters : float or array  Plate separation in meters
    E0_eV : float  E₀ in eV (default: 3.2 meV)

    Returns
    -------
    P_qvg    : float  Pressure in Pa
    delta    : float  Relative correction δP/P = ℏc/(d E₀)
    d_star   : float  Characteristic scale ℏc/E₀ in meters
    """
    d_star_m  = HBAR_C_EV_M / E0_eV
    delta     = d_star_m / d_meters
    P_std     = casimir_standard(d_meters)
    P_qvg     = P_std * (1 + delta)
    return P_qvg, delta, d_star_m


def print_casimir_table():
    """Print correction table for various separations."""
    separations = [0.1e-6, 1e-6, 10e-6, D_STAR, 100e-6, 1e-3]
    labels = ['0.1 μm', '1 μm', '10 μm', f'{D_STAR*1e6:.0f} μm (d*)',
              '100 μm', '1 mm']

    print(f"\nSPECTRAL CASIMIR EFFECT — QVG PREDICTION")
    print(f"Characteristic scale: d* = ℏc/E₀ = {D_STAR*1e6:.1f} μm")
    print("─"*70)
    print(f"{'d':15s} {'P_std (Pa)':15s} {'δP/P':10s} {'Regime':20s}")
    print("─"*70)

    for d, label in zip(separations, labels):
        P_std = casimir_standard(d)
        _, delta, _ = casimir_qvg(d)
        if delta > 10:
            regime = "QVG dominant"
        elif delta > 0.5:
            regime = "TRANSITION ← observe here"
        elif delta > 0.05:
            regime = "QVG correction visible"
        else:
            regime = "QFT standard"
        print(f"{label:15s} {P_std:15.3e} {delta:10.3f} {regime:20s}")

    print("─"*70)
    print(f"\nNOTE: The QVG predicts a fundamental change of regime at d ~ {D_STAR*1e6:.0f} μm.")
    print(f"For d < d*, the standard Casimir formula is no longer valid.")
    print(f"This is a stronger prediction than a small correction.")


if __name__ == '__main__':
    print_casimir_table()
