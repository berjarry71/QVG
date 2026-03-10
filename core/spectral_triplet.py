"""
spectral_triplet.py — Internal spectral triplet (A_F, H_F, D_F)
================================================================
N = 96 = 2 × 3 generations × 16 modes/generation.
(particle + antiparticle via real structure J)

ATTRIBUTION: A_F = C ⊕ H ⊕ M_3(C), gauge group SU(3)×SU(2)×U(1),
and sin²θ_W = 3/8 are established results of Connes-Chamseddine (1997).

Key relation verified here:
    Tr_HF(Y²) = (5/3) Tr_HF(T₃²)
This GUT normalization of U(1) in A_F gives sin²θ_W = 3/8.

References:
  [1] A. Chamseddine, A. Connes, Commun. Math. Phys. 186, 731 (1997).
  [2] A. Chamseddine, A. Connes, M. Marcolli, Adv. Theor. Math. Phys. 11 (2007).
"""

import numpy as np


def build_fermion_modes():
    """
    Build the 96-dimensional fermion mode list for H_F.
    N = 96 = 2 × 48 (particles + antiparticles via J).

    Returns list of dicts: name, gen, Y, T3, color, chirality, ptype
    """
    modes = []
    for ptype, sign in [('particle', +1), ('antiparticle', -1)]:
        for gen in range(1, 4):
            # Leptons (4 modes per generation)
            modes.append({'name': f'nuL_g{gen}_{ptype}', 'gen': gen,
                           'Y': sign*(-0.5), 'T3': sign*(+0.5),
                           'color': 0, 'chirality': 'L', 'ptype': ptype})
            modes.append({'name': f'eL_g{gen}_{ptype}', 'gen': gen,
                           'Y': sign*(-0.5), 'T3': sign*(-0.5),
                           'color': 0, 'chirality': 'L', 'ptype': ptype})
            modes.append({'name': f'eR_g{gen}_{ptype}', 'gen': gen,
                           'Y': sign*(-1.0), 'T3': 0.0,
                           'color': 0, 'chirality': 'R', 'ptype': ptype})
            modes.append({'name': f'nuR_g{gen}_{ptype}', 'gen': gen,
                           'Y': 0.0, 'T3': 0.0,
                           'color': 0, 'chirality': 'R', 'ptype': ptype})
            # Quarks (12 modes per generation: 4 types × 3 colors)
            for col in range(1, 4):
                modes.append({'name': f'uL_g{gen}_c{col}_{ptype}', 'gen': gen,
                               'Y': sign*(+1/6), 'T3': sign*(+0.5),
                               'color': col, 'chirality': 'L', 'ptype': ptype})
                modes.append({'name': f'dL_g{gen}_c{col}_{ptype}', 'gen': gen,
                               'Y': sign*(+1/6), 'T3': sign*(-0.5),
                               'color': col, 'chirality': 'L', 'ptype': ptype})
                modes.append({'name': f'uR_g{gen}_c{col}_{ptype}', 'gen': gen,
                               'Y': sign*(+2/3), 'T3': 0.0,
                               'color': col, 'chirality': 'R', 'ptype': ptype})
                modes.append({'name': f'dR_g{gen}_c{col}_{ptype}', 'gen': gen,
                               'Y': sign*(-1/3), 'T3': 0.0,
                               'color': col, 'chirality': 'R', 'ptype': ptype})

    assert len(modes) == 96, f"Expected 96 modes, got {len(modes)}"
    return modes


def compute_sin2_thetaW(rho, modes):
    """
    Compute sin²θ_W using the Connes-Chamseddine GUT normalization.

    The key identity Tr_HF(Y²) = (5/3) Tr_HF(T₃²) gives:
        sin²θ_W = Tr_ρ(T₃²) / [Tr_ρ(T₃²) + Tr_ρ(Y²)]
    which yields exactly 3/8 for uniform ρ (Connes-Chamseddine 1997).

    Parameters
    ----------
    rho : ndarray shape (96,)
    modes : list of 96 dicts

    Returns
    -------
    sin2_w : float
    """
    Y2   = np.array([m['Y']**2  for m in modes])
    T3_2 = np.array([m['T3']**2 for m in modes])
    tr_Y2   = np.dot(rho, Y2)
    tr_T3_2 = np.dot(rho, T3_2)
    denom = tr_T3_2 + tr_Y2
    if denom == 0:
        return 0.0
    return tr_T3_2 / denom


def verify_connes_identity(modes):
    """
    Verify Tr_HF(Y²) = (5/3) Tr_HF(T₃²) — the Connes-Chamseddine key relation.
    This is a structural property of A_F = C ⊕ H ⊕ M_3(C).
    """
    Y2   = np.array([m['Y']**2  for m in modes])
    T3_2 = np.array([m['T3']**2 for m in modes])
    tr_Y2   = Y2.sum()
    tr_T3_2 = T3_2.sum()
    ratio   = tr_Y2 / tr_T3_2 if tr_T3_2 > 0 else 0
    expected = 5/3
    ok = abs(ratio - expected) < 1e-6
    return ratio, expected, ok


if __name__ == '__main__':
    modes = build_fermion_modes()
    print(f"N = {len(modes)} fermionic modes  ✓")

    # Connes-Chamseddine identity
    ratio, expected, ok = verify_connes_identity(modes)
    status = '✓' if ok else '✗'
    print(f"Tr(Y²)/Tr(T₃²) = {ratio:.6f}  [expected 5/3 = {expected:.6f}]  {status}")

    # sin²θ_W at uniform ρ
    rho_unif = np.ones(96) / 96
    s = compute_sin2_thetaW(rho_unif, modes)
    status2 = '✓' if abs(s - 3/8) < 1e-6 else '✗'
    print(f"sin²θ_W(uniform ρ) = {s:.6f}  [Connes: 3/8 = {3/8:.6f}]  {status2}")
