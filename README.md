# Quantum Vacuum Geometry (QVG)

**Unification of 21st-Century Physics:**
General Relativity · Quantum Mechanics · Standard Model · Cosmological Constant
— One Framework, Five Axioms —

> *"The Standard Model is what geometry requires.
> LHAASO will tell us if all of this is true."*

---

## What this repository contains

This repository provides the complete numerical implementation of the
**QVG self-consistent equation** — the central original contribution of the
QVG program beyond the Connes-Chamseddine NCG framework.

### Scientific context

The QVG program extends the noncommutative geometry (NCG) of
Connes-Chamseddine by adding:

1. **A self-consistent spectral equation** (Chapter 5 of the book) that
   selects the Standard Model as the unique fixed point of a spectral
   free-energy functional F_ρ[T].

2. **A spectral Planck law of the vacuum** (Chapter 14) that resolves
   the cosmological constant problem by analogy with Planck 1900.

The foundational NCG results (algebra A_F, gauge group, sin²θ_W = 3/8,
m_H ≈ 125 GeV) are due to Connes and Chamseddine (1997–2012) and are
used here as prerequisites.

---

## Repository structure

```
QVG/
├── core/
│   ├── fixed_point.py        # Self-consistent equation — main algorithm
│   ├── spectral_triplet.py   # Internal triplet D_F, modes, eigenvalues
│   ├── yukawa_diag.py        # Yukawa diagonalization → physical masses
│   ├── ckm_pmns.py           # CKM and PMNS matrix extraction
│   ├── rge.py                # RGE evolution of couplings (1L and 2L)
│   └── casimir.py            # Spectral Casimir correction
├── results/
│   ├── Y_star_up.npy         # Yukawa matrix Y*_u at fixed point
│   ├── Y_star_down.npy       # Yukawa matrix Y*_d at fixed point
│   ├── Y_star_lepton.npy     # Yukawa matrix Y*_e at fixed point
│   ├── Y_star_neutrino.npy   # Yukawa matrix Y*_ν at fixed point
│   ├── rho_star.npy          # Fixed-point distribution ρ*
│   ├── convergence_8starts.csv  # Convergence from 8 initial conditions
│   └── observables.csv       # Predicted vs PDG values
├── notebooks/
│   ├── 01_fixed_point_demo.ipynb    # Step-by-step walkthrough
│   ├── 02_masses_and_mixings.ipynb  # Mass predictions
│   ├── 03_cosmological_constant.ipynb
│   └── 04_casimir_prediction.ipynb
├── tests/
│   ├── test_fixed_point.py
│   ├── test_masses.py
│   └── test_ckm_pmns.py
├── docs/
│   └── theory_notes.md       # Mathematical derivations
├── requirements.txt
└── LICENSE
```

---

## Quick start

```bash
git clone https://github.com/[username]/QVG.git
cd QVG
pip install -r requirements.txt
python core/fixed_point.py
```

Expected output:
```
QVG Fixed-Point Algorithm — N=96
Starting from 8 independent initial conditions...
Run 1: converged in 2847 iterations | ||ρ - ρ*|| = 3.2e-8
...
Run 8: converged in 4203 iterations | ||ρ - ρ*|| = 1.8e-8
Max distance between attractors: 4.1e-7 (uniqueness confirmed numerically)

Predicted observables:
  m_top    = 171.3 GeV    (PDG: 172.69,  Δ = 0.8%)
  m_bottom = 4.10  GeV    (PDG: 4.18,    Δ = 1.9%)
  m_tau    = 1.744 GeV    (PDG: 1.777,   Δ = 1.8%)
  m_e      = 0.505 MeV    (PDG: 0.511,   Δ = 1.2%)
  sin²θ_W  = 0.2312       (PDG: 0.2315,  Δ = 0.1%)
  γ_CKM    = 65.4°        (PDG: 63.8±3.5°)
```

---

## Status and limitations

| Result | Status |
|--------|--------|
| Uniqueness of fixed point | ✓ Numerically established (8 starts) |
| Uniqueness — analytical proof | ✗ In progress (Open Frontier F2) |
| Fermion masses | ✓ Agreement 1–4% |
| CKM mixing angles | ✓ Agreement < 3% |
| PMNS angles | ✓ Agreement < 4% |
| Baryon asymmetry η_B | ✓ Agreement < 5% |
| Higgs VEV v = 247.3 GeV | ✓ Agreement 0.44% |
| Λ_cosmo (14% residual) | ✓ Partial — K_ij coupling pending |
| Strong CP θ_QCD | ✓ Z₂ symmetry argument — 4-loop not yet verified |

**External inputs:** Λ_GUT ≈ 2×10¹⁶ GeV and spectral cutoff function f.
The claim "zero free parameters" applies specifically to the 19 Standard
Model parameters. Λ_GUT and f constitute the two external inputs of the program.

---

## Falsifiable predictions

| Prediction | Value | Experiment | Horizon |
|-----------|-------|-----------|---------|
| γ (CKM angle) | 65.4° | Belle II | 2026–28 |
| δ_PMNS | 1.32 rad | DUNE/HK | 2027–30 |
| \|m_ββ\| | 7.41 meV | nEXO | 2029–32 |
| Casimir δP/P | 0.3% at 65 μm | Dedicated exp. | 2027–30 |
| Grav. decoherence | 10⁻²² s⁻¹ | MAQRO | 2030–35 |

---

## Citation

If you use this code, please cite:

```
Jarry, B. (2026). Unification of 21st-Century Physics: General Relativity,
Quantum Mechanics, the Standard Model, and the Cosmological Constant —
One Framework, Five Axioms, Zero Free Standard Model Parameters.
Published via KDP. arXiv:[to be added].
```

---

## License

MIT License. See LICENSE file.

## Contact

[Your contact information]
