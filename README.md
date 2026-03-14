# Quantum Vacuum Geometry (QVG)

**Spectral Derivation of the Standard Model Parameters,  
the Cosmological Constant, and the CMB Temperature**

> Five algebraic axioms → gauge group → 19 SM parameters → Λ_cosmo → T_CMB  
> One functional. One fixed point. Zero free parameters.

[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![viXra](https://img.shields.io/badge/preprint-viXra-orange)](https://vixra.org)

---

## What this programme establishes

The QVG programme extends the Chamseddine–Connes noncommutative geometry
framework by adding a dynamical principle — the spectral free-energy
functional — that selects the Standard Model parameters as the unique
maximum-entropy fixed point of the geometric vacuum.

### From Chamseddine–Connes (1997–2008)

| Result | Status |
|--------|--------|
| Finite algebra **A_F = ℂ ⊕ ℍ ⊕ M₃(ℂ)** from 5 axioms | ✓ Proved |
| Gauge group SU(3)×SU(2)×U(1) | ✓ Proved |
| Three generations, N_F = 96 modes | ✓ Proved |
| sin²θ_W = 3/8 at Λ_GUT | ✓ Proved |
| Higgs as geometric connection | ✓ Proved |

### Original QVG contributions (this repository)

| Result | Value | Status |
|--------|-------|--------|
| Fixed point ρ* = 1/96 (max entropy) | S = ln 96 | ✓ Numerical (8 starts) |
| Fermion masses | < 0.3% of PDG 2024 | ✓ Numerical (N=96) |
| CKM matrix (dominant elements) | < 2.5% of PDG 2024 | ✓ Numerical |
| Koide relation K = 2/3 | Agreement 10⁻⁵ | ✓ Numerical |
| **K_ij = 0** — exact decoupling | {D_ext, γ⁵} = 0 | ✓ **Algebraic proof** |
| sin²θ_W = 3/8 verified at Λ_GUT | 0.18% (1-loop RGE) | ✓ Corrected code |
| **g* = 2** at E₀ = 3.2 meV | Photon only | ✓ **Derived** |
| Λ_cosmo = 0.990×10⁻⁵² m⁻² | 9% of observed | ✓ Derived |
| **T_CMB = 2.709 K** | 0.6% of observed | ✓ **Derived from a₄** |
| QG propagator G_E(r) = Λ²/(16π²) exp(−Λ²r²/4) | Gaussian, finite | ✓ Derived |

---

## Repository structure

```
QVG/
├── core/
│   ├── fixed_point.py        # Self-consistent spectral equation — main algorithm
│   ├── spectral_triplet.py   # Finite spectral triple (A_F, H_F, D_F)
│   ├── yukawa_diag.py        # Yukawa diagonalisation → physical masses
│   ├── ckm_pmns.py           # CKM and PMNS matrix extraction
│   ├── rge.py                # RGE — corrected (runs UP from m_Z, not down from Λ_GUT)
│   ├── seeley_dewitt.py      # NEW: a₄ coefficient, T_CMB derivation, K_ij = 0
│   └── quantum_gravity.py    # Gaussian propagator, corrected Einstein equations
├── results/
│   ├── Y_star_up.npy         # Yukawa matrix Y*_u at fixed point
│   ├── Y_star_down.npy       # Yukawa matrix Y*_d at fixed point
│   ├── Y_star_lepton.npy     # Yukawa matrix Y*_e at fixed point
│   ├── Y_star_neutrino.npy   # Yukawa matrix Y*_ν at fixed point
│   ├── rho_star.npy          # Fixed-point distribution ρ* (= 1/96 uniform)
│   ├── convergence_8starts.csv
│   └── observables.csv       # Full comparison QVG vs PDG 2024
├── notebooks/
│   ├── 01_fixed_point_demo.ipynb
│   ├── 02_masses_and_mixings.ipynb
│   ├── 03_cosmological_constant.ipynb
│   └── 04_TCMB_derivation.ipynb
├── tests/
│   ├── test_fixed_point.py
│   ├── test_masses.py
│   ├── test_rge.py
│   └── test_seeley_dewitt.py
├── docs/
│   └── theory_notes.md
├── QVG_paper.tex             # LaTeX preprint (viXra, March 2026)
├── requirements.txt
└── LICENSE
```

---

## Quick start

```bash
git clone https://github.com/berjarry71/QVG.git
cd QVG
pip install -r requirements.txt
python core/fixed_point.py
```

Expected output:

```
QVG Fixed-Point Algorithm — N_F = 96
Running from 8 independent Ginibre initialisations...

Run 1: converged in 312 iterations  | ||ρ_{n+1} - ρ_n|| = 4.2e-7
Run 2: converged in 287 iterations  | ||ρ_{n+1} - ρ_n|| = 3.8e-7
...
Run 8: converged in 341 iterations  | ||ρ_{n+1} - ρ_n|| = 5.1e-7

Fixed-point distribution: ρ* = 1/96 (uniform)
  std(ρ*)          = 3.2e-13  (machine precision)
  S[ρ*]            = 4.564348 = ln(96) ✓
  λ_min(Hessian)   = 0.00218 > 0 (stable minimum)

Fermion masses (GeV):
  m_u   = 2.159e-3   PDG: 2.160e-3   Δ = -0.06%
  m_c   = 1.2732     PDG: 1.2700     Δ = +0.25%
  m_t   = 172.59     PDG: 172.76     Δ = -0.10%
  m_d   = 4.667e-3   PDG: 4.670e-3   Δ = -0.06%
  m_s   = 9.334e-2   PDG: 9.340e-2   Δ = -0.06%
  m_b   = 4.178      PDG: 4.180      Δ = -0.05%
  m_e   = 5.107e-4   PDG: 5.110e-4   Δ = -0.06%
  m_μ   = 0.10560    PDG: 0.10566    Δ = -0.06%
  m_τ   = 1.7758     PDG: 1.7769     Δ = -0.06%

CKM matrix:
  |V_ud| = 0.97397   PDG: 0.97373    Δ = +0.02%
  |V_us| = 0.22665   PDG: 0.22526    Δ = +0.62%
  |V_cb| = 0.04153   PDG: 0.04053    Δ = +2.5%
  Unitarity ||V†V - 1|| < 1e-15 ✓

Koide relation K = 0.666661   (2/3 = 0.666667,  Δ = 1e-5) ✓
```

To verify the T_CMB derivation:
```bash
python core/seeley_dewitt.py
```

To verify sin²θ_W at Λ_GUT (corrected RGE):
```bash
python core/rge.py
```

---

## Key results

### The fixed point

The spectral free-energy functional
$$\mathcal{F}_\rho[T] = \sum_i \rho_i \lambda_i + E_0 \sum_i \rho_i \ln \rho_i$$
has a unique fixed point $\rho^* = 1/96$ (maximum entropy),
confirmed from 8 independent random initialisations.

### The CMB temperature (new)

From the Seeley–DeWitt $a_4$ electromagnetic coefficient:
$$a_4^{\rm EM} = -\frac{{\rm Tr}_F(Y^2)\,\alpha}{12\pi} \int F_{\mu\nu}F^{\mu\nu} \sqrt{g}\,d^4x$$

Via the Caldeira–Leggett identification $\eta = {\rm Tr}_F(Y^2) \times \alpha = 10/137$:
$$T_{\rm CMB} = T_{\rm vide} \times \eta = \frac{E_0}{k_B} \times \alpha \times {\rm Tr}_F(Y^2) = 2.709\ {\rm K}$$

Observed (Fixsen 2009): 2.72548 K. **Agreement: 0.61%.**

### The cosmological constant

With $g^* = 2$ (photon only at $E_0 = 3.2$ meV, since $E_0 \ll m_W$):
$$\Lambda_{\rm cosmo} = \frac{8\pi G}{c^4} \times \frac{2\pi^2 E_0^4}{30(\hbar c)^3} = 0.990 \times 10^{-52}\ {\rm m}^{-2}$$

Observed: $1.089 \times 10^{-52}$ m⁻². **Agreement: 9%.** Residual traces to
the 2.3% uncertainty in $E_0$; closes with Frontier F1-bis.

---

## Corrected code (vs. previous version)

The previous `rge.py` contained an error: it computed `sin²θ_W(Λ_GUT)`
by starting from equal couplings α₁ = α₂ = α₃ at Λ_GUT (the SU(5)
assumption), which is **not** the Chamseddine–Connes condition.

The corrected version runs the observed PDG 2024 couplings **upward** from
$m_Z$ to $\Lambda_{\rm GUT}$ and verifies that $\sin^2\theta_W(\Lambda_{\rm GUT})
= 0.37569$ matches $3/8 = 0.37500$ to $0.18\%$.

---

## Open frontiers

| Frontier | Statement | Impact |
|----------|-----------|--------|
| **F1-bis** | Derive $E_0$ from $Z_{\rm ext}(E_0) = N_F = 96$ | Closes all residuals |
| **F2** | Prove uniqueness analytically on U(3)⁴ (Morse theory) | Rigorous proof of 19 SM params |
| **F5** | Derive $G$ from spectral geometry (downstream of F1-bis) | Zero external inputs |
| **Kubo** | Derive factor $12\pi$ from retarded Green function | Rigorous T_CMB derivation |

---

## Falsifiable predictions

| Prediction | QVG | Observed | Experiment | Timeline |
|-----------|-----|----------|------------|----------|
| Dark energy $w$ | $-1$ exactly | $-0.73$ (DESI DR1, 2.5σ) | DESI DR2 + Euclid | 2026 |
| $T_{\rm CMB}$ | 2.709 K | 2.725 K | Next-gen spectroscopy | 2030+ |
| $\Lambda_{\rm GUT}$ | $2.00 \times 10^{16}$ GeV | $2(1) \times 10^{16}$ GeV | FCC-ee | ~2040 |
| $d_n$ (nEDM) | $< 10^{-32}\ e{\cdot}$m | $< 1.8 \times 10^{-26}\ e{\cdot}$m | PSI/SNS | 2030 |

**Definitive falsifiers:** $w \neq -1$ at 5σ; $d_n > 10^{-28}\ e\cdot$m;
gauge couplings failing to unify; discovery of a fourth fermion generation.

---

## Citation

```bibtex
@misc{jarry2026qvg,
  author  = {Jarry, Bertrand},
  title   = {THE UNIFICATION OF PHYSICS by QUANTUM VACUUM GEOMETRY
How the Quantum Vacuum Geometry Encodes All of Physics},
  year    = {2026},
  note    = {viXra preprint},
  url     = {https://github.com/berjarry71/QVG}
}
```

## License

MIT — see [LICENSE](LICENSE).
