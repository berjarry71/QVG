# QVG Theory Notes

Mathematical derivations and physical arguments underlying the
Quantum Vacuum Geometry programme. For the complete paper see
`QVG_paper.tex` in the root directory.

---

## 1. The finite algebra A_F

From five algebraic axioms C1–C5 on the finite spectral triple
(A_F, H_F, D_F, J_F, γ_F), the unique algebra of minimal dimension is:

```
A_F = ℂ ⊕ ℍ ⊕ M₃(ℂ)
```

This generates:
- Gauge group: SU(3) × SU(2) × U(1)
- Three generations of fermions: N_F = 96 modes
- Hypercharge trace: Tr_F(Y²) = 10 (exact, algebraic)
- Weinberg angle: sin²θ_W = 3/8 at Λ_GUT (exact, from trace formula)

**Reference:** Chamseddine, Connes, Marcolli (2007), Adv. Theor. Math. Phys. 11, 991.

---

## 2. The spectral action and Seeley–DeWitt expansion

The spectral action:
```
S_fond[D, Λ] = Tr[f(D²/Λ²)]
```

Seeley–DeWitt expansion for large Λ:
```
Tr[f(D²/Λ²)] = f₂ Λ⁴ a₀ + f₀ Λ² a₂ + a₄ + O(Λ⁻²)
```

where f₂ = ∫₀^∞ f(u) du, f₀ = f(0), and for D² = -(∇² + E):

```
a₀ = N/(16π²) Vol(M)
a₂ = 1/(16π²) ∫ Tr(E + R/6) √g d⁴x
a₄ = 1/(16π²) ∫ Tr[E²/2 - RE/6 + Ω_μν Ω^μν/12 + ...] √g d⁴x
```

**Reference:** Vassilevich (2003), Phys. Rept. 388, 279.

---

## 3. The exact decoupling K_ij = 0

**Proposition:** For the product triplet D = D_ext ⊗ 1 + γ⁵ ⊗ D_F:
```
D² = D²_ext ⊗ 1 + 1 ⊗ D²_F   (exactly)
```

**Proof:**
```
D² = D²_ext ⊗ 1 + {D_ext, γ⁵} ⊗ D_F + 1 ⊗ D²_F

{D_ext, γ⁵} = i {γᵘ, γ⁵} ∇_μ = 0   (Clifford algebra)
```

**Consequence:** The cosmological constant (from D²_ext) is independent
of the fermion masses (from D²_F) at all orders. No fine-tuning required.

---

## 4. The spectral free-energy functional

```
F_ρ[T] = Σᵢ ρᵢ λᵢ + E₀ Σᵢ ρᵢ ln ρᵢ
        = ⟨λ⟩ - E₀ S[ρ]
```

where λᵢ are the eigenvalues of T = D_F†D_F / Λ²_GUT.

**Fixed-point condition:**
```
ρ* = ρ*[T*(Y*(ρ*))]
```

**Fixed point:** ρ*ᵢ = 1/96 (uniform), S[ρ*] = ln(96) = maximum entropy.

**Why uniform?** All fermion masses m_f >> E₀ (lightest is electron at
0.511 MeV >> E₀ = 3.2 meV), so all Boltzmann weights are exponentially
suppressed and approximately equal. Correction: O(m_e²/Λ_GUT E₀) ~ 10⁻²¹.

---

## 5. The Koide relation K = 2/3

The fixed-point lepton Yukawa matrix Y*_e satisfies:

```
K = (m_e + m_μ + m_τ) / (√m_e + √m_μ + √m_τ)² = 2/3  ±  10⁻⁵
```

This is not imposed — it emerges from the entropy-maximum structure of Y*_e.
QVG gives K = 0.666661 vs 2/3 = 0.66667. Agreement: 10⁻⁵.

---

## 6. Newton's constant from spectral data

From the a₂ term of the Seeley–DeWitt expansion:

```
G = 3π / (f₂ × N_F × Λ²_GUT)
```

With N_F = 96, Λ_GUT = 2×10¹⁶ GeV, G = 6.674×10⁻¹¹ SI:
```
f₂ = 3π / (G × 96 × (2×10¹⁶ GeV)²) = 73,168  [natural units]
```

---

## 7. The cosmological constant

Active bosonic degrees of freedom at E₀ = 3.2 meV:
- Photon: g* = 2 (massless, 2 polarisations) ✓
- W±, Z⁰, H: inactive (masses >> E₀) ✗
- Gluons: confined at Λ_QCD >> E₀ ✗
- → **g* = 2** (derived from electroweak symmetry breaking of A_F)

Stefan–Boltzmann vacuum energy:
```
ρ_vac = g* π² E₀⁴ / (30 (ħc)³) = 1.437×10⁻¹⁰ J/m³
```

Cosmological constant:
```
Λ_cosmo = (8πG/c⁴) × ρ_vac = 0.990×10⁻⁵² m⁻²
```

Observed: 1.089×10⁻⁵² m⁻². **Agreement: 9%.** Residual from 2.3% uncertainty in E₀.

---

## 8. The CMB temperature (new, March 2026)

### 8.1 The a₄ electromagnetic coefficient

With the EM curvature Ω_μν = ie F_μν and hypercharge weighting:

```
a₄^EM = -Tr_F(Y²) × α / (12π) × ∫ F_μν F^μν √g d⁴x
```

**Key:** Use Tr_F(Y²) = 10, not N_F = 96. The photon couples to fermion i
with weight Yᵢ², so the effective count is Σᵢ nᵢ Yᵢ² = Tr_F(Y²) = 10.

### 8.2 The Caldeira–Leggett identification

Ohmic spectral density from a₄^EM:
```
J(ω) = 12π × c_EM × ω = Tr_F(Y²) × α × ω
```

Thermal coupling:
```
η = lim_{ω→0} J(ω)/ω = Tr_F(Y²) × α = 10/137 = 0.0730
```

### 8.3 The temperature formula

In the weak-coupling regime η << 1 (Caldeira–Leggett 1983):
```
T_CMB = T_vide × η + O(η²)
      = (E₀/k_B) × Tr_F(Y²) × α
      = 37.12 K × (10/137)
      = 2.709 K
```

Observed (Fixsen 2009): 2.72548 K. **Agreement: 0.61%.**

### 8.4 Why linear in η?

The linear law n=1 is the unique consistent power law (see seeley_dewitt.py).
Also forced by: (i) a₄^EM linear in α (first-order perturbation theory);
(ii) standard Caldeira–Leggett weak-coupling result.

---

## 9. The quantum gravity propagator

Path integral over Dirac operators:
```
Z_grav = ∫ 𝒟[D] exp(-Tr[f(D²/Λ²)] / ħ)
```

UV-finite because f suppresses |λ| > Λ². Saddle-point expansion gives
the Euclidean propagator:
```
G_E(r) = (Λ²_GUT / 16π²) × exp(-Λ²_GUT r²/4)
```

Finite at r=0: G_E(0) = Λ²/16π². Standard graviton propagator ~1/r² → ∞.

Quantum-corrected Einstein equations:
```
G_μν + Λ_cosmo g_μν
  + (α_g/Λ²) R_μν
  + (β_g/Λ²) R g_μν = (8πG/c⁴) T_μν
```

with exact Seeley–DeWitt coefficients:
```
α_g = -N/(48π²) = -0.00844
β_g = +N/(96π²) = +0.00422
```

---

## 10. Open frontiers

| Frontier | Statement | Method | Impact |
|----------|-----------|--------|--------|
| F1-bis | Z_ext(E₀) = N_F = 96 | Discrete mode sum on compact M | Closes all residuals |
| F2 | Unique minimum on U(3)⁴ | Equivariant Morse theory | Rigorous proof of 19 SM params |
| F5 | Derive G from geometry | Downstream of F1-bis | Zero external inputs |
| Kubo | Derive factor 12π in J(ω) | One-loop OPE calculation | Rigorous T_CMB derivation |

---

## 11. Corrected code (March 2026)

### RGE correction

The previous `rge.py` used α₁ = α₂ = α₃ at Λ_GUT as starting condition
(SU(5) assumption). This is **not** the Chamseddine–Connes condition.

**Correct method:** Run observed PDG 2024 values *upward* from m_Z:
```python
alpha_i(Lambda) = alpha_i(m_Z) / (1 - b_i * alpha_i(m_Z) / (2*pi) * ln(Lambda/m_Z))
```

Result: sin²θ_W(Λ_GUT) = 0.37569 vs 3/8 = 0.37500. Agreement: 0.18%.

### a₄ coefficient correction

Previous version used N_F × α/(12π) for the EM coupling.
Correct value: Tr_F(Y²) × α/(12π) (hypercharge weighting).
Factor: Tr_F(Y²)/N_F = 10/96 = 0.104 (not 1).

---

*Last updated: March 2026*
*Reference: Jarry, B. (2026). Quantum Vacuum Geometry. viXra preprint.*
*Code: github.com/berjarry71/QVG*
