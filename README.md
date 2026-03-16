# Quantum Vacuum Geometry (QVG)

**Author:** Bertrand Jarry  
**Preprint:** *Quantum Vacuum Geometry: Spectral Derivation of the Standard Model Parameters, the Cosmological Constant, and the CMB Temperature* (2026)  
**Repository:** `github.com/berjarry71/QVG`

---

## Overview

QVG derives the physical constants of the Standard Model from five algebraic axioms (C1–C5) applied to a non-commutative spectral triple. The finite algebra **A_F = ℂ ⊕ ℍ ⊕ M₃(ℂ)** forces N_F = 96 fermion modes; the spectral free energy functional F_ρ[T] has a unique maximum-entropy fixed point ρ* = 1/96 from which all SM parameters emerge.

---

## Results: Three Levels of Certainty

### Level 1 — Algebraically proven (exact, from axioms C1–C5)

| Result | Value | Derivation |
|---|---|---|
| Fermion modes N_F | 96 | dim(H_F) from A_F |
| Hypercharge invariant Tr_F(Y²) | 10 | Trace over H_F |
| Weinberg angle sin²θ_W at Λ_GUT | 3/8 | Tr_F(T₃²)/Tr_F(Y²)_GUT |
| Off-diagonal curvature K_ij | 0 | Algebraic proof |
| Loop closure condition | g* = N_F/e at T = 149 MeV | From 4π³g*/(45N_F) = 1 |

### Level 2 — Numerically established (computed, reproducible)

| Result | QVG | Observed | Deviation |
|---|---|---|---|
| ρ* = 1/96 (fixed point) | 1/96 exact | — | ✓ Machine precision |
| T_CMB | 2.709 K | 2.7255 K | 0.61% |
| Λ_cosmo | 0.990×10⁻⁵² m⁻² | 1.089×10⁻⁵² m⁻² | 9% |
| sin²θ_W (numerical RGE) | 0.37569 | 3/8 = 0.37500 | 0.18% |
| Koide K (emergent) | 0.666661 | 2/3 = 0.666667 | 4×10⁻⁶ |

### Level 3 — Demonstration models (illustrate conjectures, not full proofs)

| Script | What it does | Key limitation |
|---|---|---|
| `fixed_point.py` | Demonstrates ρ → 1/N_F convergence | M-step simplified; masses depend on random seed |
| `true_qvg_solver.py` | FSS O(1/N) scaling toy model | PDG masses injected as targets — tautological |
| `rge_2loop_solver_corrected.py` | 2-loop RGE Runge-Kutta | GUT-normalised α₁ has Landau pole; results unreliable |
| `spectral_anomaly.py` | Trace anomaly h(T) peak at 155 MeV | m_gap = ms√π injected, not derived |
| `mcmc_lattice.py` | Lattice NCG MCMC prototype | Structural prototype; does not produce SM masses |

---

## The Single Open Frontier

One calculation remains to close the theoretical loop:

**Derive ⟨a₄(D_ext, T)⟩_KMS from M₃(ℂ) during confinement — without injecting any phenomenological mass parameter.**

If proven, this derives E₀ = 3.2 meV from the axioms alone and closes:

**R_H ↔ E₀ ↔ Λ_cosmo ↔ R_H**

---

## Falsifiable Prediction

QVG predicts **w = −1** exactly. DESI DR3 (2026–2027) is the experimental verdict.  
Current status: DESI DR2 (March 2025) shows 2.8–4.2σ preference for w ≠ −1.

---

## Repository Structure

```
core/
  fixed_point.py                    # Demonstrates ρ*=1/96 (M-step simplified)
  rge.py                            # 1-loop RGE running
  seeley_dewitt.py                  # a₄, Caldeira-Leggett, T_CMB
  observables.csv                   # QVG vs PDG 2024

research/
  spectral_anomaly.py               # h(T) proof-of-concept (m_gap injected)
  rge_2loop_solver_corrected.py     # 2-loop prototype (see limitations)
  mcmc_lattice.py                   # Lattice NCG prototype
  true_qvg_solver.py                # FSS toy model (PDG targets injected)
  theory_notes.md                   # Full theoretical derivations
```

---

## Running

```bash
pip install numpy scipy matplotlib
python core/fixed_point.py        # ρ* = 1/96 convergence
python core/seeley_dewitt.py      # T_CMB = 2.709 K
python research/spectral_anomaly.py  # h(T) peak at 155 MeV
```

---

## Citation

```
Jarry, B. (2026). Quantum Vacuum Geometry. viXra. github.com/berjarry71/QVG
```
