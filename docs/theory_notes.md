# QVG Theory Notes

## What is original in QVG vs Connes-Chamseddine NCG

### Connes-Chamseddine results (prior work, used as foundation)

The following results are established in the NCG literature and used
as the starting point of the QVG program:

| Result | Reference |
|--------|-----------|
| A_F = ℂ ⊕ ℍ ⊕ M₃(ℂ) from C1–C5 | Chamseddine-Connes (1996, 2008) |
| Gauge group SU(3)×SU(2)×U(1) from A_F | Chamseddine-Connes (1997) |
| sin²θ_W = 3/8 at Λ_GUT | Chamseddine-Connes (1997) |
| S_fond = Tr[f(D²/Λ²)] and Seeley-DeWitt expansion | Chamseddine-Connes (1997) |
| m_H ≈ 125 GeV prediction | Chamseddine-Connes (2010, 2012) |
| Neutrino mixing via see-saw in D_F | Chamseddine-Connes-Marcolli (2007) |

### QVG original contributions

| Result | Where |
|--------|-------|
| Self-consistent spectral equation F_ρ[T] | Chapter 5 |
| Fixed-point selection of SM parameters | Chapter 5–9 |
| Uniqueness via Morse theory on U(3)⁴ | Chapter 10 |
| Strong CP resolution via Z₂(A_F) | Chapter 11 |
| Spectral Planck law of the vacuum | Chapter 14 |
| Λ_cosmo from thermodynamic pressure | Chapter 15 |
| Spectral Casimir correction δP/P | Chapter 20 |

## External inputs of the program

The QVG program has **two external inputs**, not zero:

1. **Λ_GUT ≈ 2×10¹⁶ GeV** — the GUT unification scale.
   Not derived from C1–C5. Taken from phenomenology of gauge coupling
   unification. A derivation of Λ_GUT from the axioms remains open.

2. **The spectral cutoff function f** (and its moment f₂ ≈ 982).
   Different choices of f give different predictions for Λ_cosmo.
   The normalization condition on ρ*_ext fixes E₀ given f and Λ_GUT,
   but f itself is not uniquely determined by the axioms.

The claim "zero free parameters" in the book applies specifically to the
**19 Standard Model parameters** (fermion masses, mixing angles, gauge
couplings). These are derived, not input.

## Status of key results

### Uniqueness of the fixed point (Chapter 10)

**Status: Numerically established. Analytical proof incomplete.**

The uniqueness argument uses:
1. λ_min(Hess F) = 0.00218 > 0  ← verified numerically
2. χ(U(3)⁴) = 0  ← exact (topological)
3. ΔF_inst = 7×10⁻⁶ > 0  ← estimated numerically, not proven analytically
                              for ALL directions in U(3)⁴

The result should be considered a well-supported numerical conjecture
until the analytical proof of (3) is complete.

### Strong CP problem (Chapter 11)

**Status: Z₂ symmetry argument established. 4-loop calculation estimated.**

The argument that θ_QCD = 0 at Λ_GUT from Z₂(A_F) invariance of S_fond
is rigorous. The subsequent generation of θ_QCD ≈ 8×10⁻¹⁹ via 4-loop
corrections is estimated by power counting. The explicit 4-loop diagrams
have not been computed. An independent verification is required.

### Cosmological constant (Chapter 15)

**Status: 14% residual — source identified.**

The distinction energy/pressure is correct but the 14% residual comes
from the coupling K_ij between ρ*_ext and ρ*_int not yet calculated
precisely. The calculation is in progress.

## References

[1] A. Chamseddine, A. Connes. The spectral action principle. 
    Commun. Math. Phys. 186, 731 (1997).

[2] A. Chamseddine, A. Connes, M. Marcolli. Gravity and the Standard 
    Model with neutrino mixing. Adv. Theor. Math. Phys. 11, 991 (2007).

[3] A. Connes. On the spectral characterization of manifolds. 
    J. Noncommut. Geom. 7, 1 (2013).

[4] PDG. Review of Particle Physics (2024).
