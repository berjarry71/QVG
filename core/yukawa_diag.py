"""
yukawa_diag.py — Yukawa diagonalization and physical mass extraction
===================================================================
Phase C of the QVG algorithm: extract physical observables from Y*.

Physical masses via biunitary diagonalization:
    V_L(f)† Y*_f V_R(f) = diag(y_f1, y_f2, y_f3)
    m_fi = y_fi × v / sqrt(2),  v = 246.22 GeV
"""

import numpy as np


HIGGS_VEV = 246.22  # GeV

# PDG 2024 reference values (GeV)
PDG_MASSES = {
    'top':    172.69,
    'bottom':   4.18,
    'charm':    1.274,
    'strange':  0.0934,
    'tau':      1.7769,
    'muon':     0.10566,
    'electron': 5.11e-4,
}


def biunitary_diag(Y):
    """
    Biunitary diagonalization of a 3×3 complex matrix Y.
    Returns (VL, VR, singular_values) such that VL† Y VR = diag(sv).
    """
    U, sv, Vh = np.linalg.svd(Y)
    VL = U
    VR = Vh.conj().T
    return VL, VR, np.sort(sv)[::-1]  # descending order


def yukawa_to_masses(Y, vev=HIGGS_VEV):
    """
    Convert Yukawa matrix to physical masses.
    m_i = y_i × vev / sqrt(2)
    """
    _, _, sv = biunitary_diag(Y)
    masses = sv * vev / np.sqrt(2)
    return masses


def mass_table(Y_u, Y_d, Y_e, vev=HIGGS_VEV):
    """
    Print comparison table: predicted masses vs PDG 2024.
    """
    m_u = yukawa_to_masses(Y_u, vev)
    m_d = yukawa_to_masses(Y_d, vev)
    m_e = yukawa_to_masses(Y_e, vev)

    quark_names   = ['top', 'charm', 'up']
    down_names    = ['bottom', 'strange', 'down']
    lepton_names  = ['tau', 'muon', 'electron']

    print("\n{'='*55}")
    print(f"{'Particle':12s} {'QVG (GeV)':12s} {'PDG (GeV)':12s} {'Δ':8s}")
    print("─"*55)

    for name, m_pred in zip(quark_names, m_u):
        if name in PDG_MASSES:
            delta = abs(m_pred - PDG_MASSES[name]) / PDG_MASSES[name] * 100
            print(f"{name:12s} {m_pred:12.4f} {PDG_MASSES[name]:12.4f} {delta:6.1f}%")

    for name, m_pred in zip(down_names, m_d):
        if name in PDG_MASSES:
            delta = abs(m_pred - PDG_MASSES[name]) / PDG_MASSES[name] * 100
            print(f"{name:12s} {m_pred:12.4f} {PDG_MASSES[name]:12.4f} {delta:6.1f}%")

    for name, m_pred in zip(lepton_names, m_e):
        if name in PDG_MASSES:
            delta = abs(m_pred - PDG_MASSES[name]) / PDG_MASSES[name] * 100
            print(f"{name:12s} {m_pred:12.4g} {PDG_MASSES[name]:12.4g} {delta:6.1f}%")

    print("─"*55)
