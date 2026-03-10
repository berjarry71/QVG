"""
ckm_pmns.py — CKM and PMNS matrix extraction
=============================================
Extract mixing angles and CP phases from Yukawa matrices at the fixed point.

CKM: V_CKM = V_L(u)† V_L(d)
PMNS: U_PMNS from neutrino see-saw mechanism
"""

import numpy as np


def extract_ckm(Y_u, Y_d):
    """
    Compute V_CKM = V_L(u)† V_L(d) from Yukawa matrices.

    Returns
    -------
    V_CKM : ndarray, shape (3,3), complex
    angles : dict with lambda, A, rho, eta (Wolfenstein)
    Jarlskog : float
    """
    # Biunitary diagonalization
    U_u, _, _ = np.linalg.svd(Y_u)
    U_d, _, _ = np.linalg.svd(Y_d)
    V_CKM = U_u.conj().T @ U_d

    # Wolfenstein parameters
    lam = abs(V_CKM[0, 1])                              # λ = |V_us|
    A   = abs(V_CKM[1, 2]) / lam**2                     # A
    rho_eta = V_CKM[0, 2].conj() / (A * lam**3)
    rho = rho_eta.real
    eta = -rho_eta.imag

    # Jarlskog invariant
    J = abs(np.imag(
        V_CKM[0,0] * V_CKM[1,1].conj() *
        V_CKM[0,1].conj() * V_CKM[1,0]
    ))

    # Unitarity triangle angles
    # γ = arg(-V_ud V_ub* / V_cd V_cb*)
    gamma_rad = np.angle(
        -V_CKM[0,0] * V_CKM[0,2].conj() /
        (V_CKM[1,0] * V_CKM[1,2].conj())
    )

    return {
        'V_CKM'  : V_CKM,
        'lambda' : lam,
        'A'      : A,
        'rho'    : rho,
        'eta'    : eta,
        'Jarlskog': J,
        'gamma_deg': np.degrees(gamma_rad) % 180,
    }


def print_ckm_summary(ckm):
    """Print CKM results with PDG comparison."""
    PDG = {'lambda': 0.2250, 'A': 0.826, 'Jarlskog': 3.08e-5,
           'gamma_deg': 63.8}

    print("\nCKM PARAMETERS")
    print("─"*50)
    for key in ['lambda', 'A', 'Jarlskog', 'gamma_deg']:
        pred = ckm[key]
        ref  = PDG.get(key, None)
        if ref:
            delta = abs(pred - ref) / abs(ref) * 100
            print(f"  {key:12s}: {pred:.4g}  (PDG: {ref:.4g}, Δ={delta:.1f}%)")
        else:
            print(f"  {key:12s}: {pred:.4g}")
    print(f"\n  QVG prediction: γ = {ckm['gamma_deg']:.1f}°")
    print(f"  Belle II target: γ = 65.4°  (test by 2028)")
