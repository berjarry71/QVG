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
│   └── casimir.py            # Spectral Casimir effect
├── results/
│   └── README.md             # Results will be added upon N=96 completion
├── notebooks/
│   └── (Jupyter notebooks — coming soon)
├── tests/
│   ├── test_fixed_point.py
│   └── test_spectral_triplet.py
├── docs/
│   └── theory_notes.md       # NCG/QVG boundary, proof status
├── requirements.txt
└── LICENSE
```

---

## Quick start

```bash
git clone https://github.com/berjarry71/QVG.git
cd QVG
pip install -r requirements.txt
python core/spectral_triplet.py
```

Expected output:
```
N = 96 fermionic modes  ✓
Tr(Y²)/Tr(T₃²) = 1.666667  [5/3 = 1.666667]  ✓
sin²θ_W(uniform ρ) = 0.375000  [Connes: 3/8 = 0.375000]  ✓
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
Model parameters. Λ_GUT and f are the two external inputs of the program.

---

## Falsifiable predictions

| Prediction | Value | Experiment | Horizon |
|-----------|-------|-----------|---------|
| γ (CKM angle) | 65.4° | Belle II | 2026–28 |
| δ_PMNS | 1.32 rad | DUNE/HK | 2027–30 |
| \|m_ββ\| | 7.41 meV | nEXO | 2029–32 |
| Spectral Casimir | regime change at d* ≈ 62 μm | Dedicated exp. | 2027–30 |
| Grav. decoherence | 10⁻²² s⁻¹ | MAQRO | 2030–35 |

**Note on the Casimir prediction:** The QVG introduces a characteristic
length scale d* = ℏc/E₀ ≈ 62 μm. For d < d*, the QVG correction to
the Casimir pressure exceeds the standard QFT value — this is a change
of physical regime, not a small perturbative correction. See `core/casimir.py`
for the full calculation.

---

## Citation

```bibtex
@book{jarry2026qvg,
  author    = {Jarry, Bertrand},
  title     = {Unification of 21st-Century Physics},
  subtitle  = {General Relativity, Quantum Mechanics, the Standard Model,
               and the Cosmological Constant — One Framework, Five Axioms,
               Zero Free Standard Model Parameters},
  year      = {2026},
  publisher = {Kindle Direct Publishing},
  note      = {arXiv: [to be added]}
}
```

---

## License

MIT License — see [LICENSE](LICENSE).

## Contact

Bertrand Jarry — berjarry71 on GitHub
