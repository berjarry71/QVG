"""
tests/test_fixed_point.py — Unit tests for the QVG fixed-point algorithm.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core'))

import numpy as np
import pytest
from spectral_triplet import build_fermion_modes, compute_sin2_thetaW
from fixed_point import boltzmann_distribution, free_energy


def test_mode_count():
    """N = 96 fermionic modes."""
    modes = build_fermion_modes()
    assert len(modes) == 96


def test_connes_chamseddine_limit():
    """Uniform ρ → sin²θ_W = 3/8 (Connes-Chamseddine result)."""
    modes = build_fermion_modes()
    rho_uniform = np.ones(96) / 96
    s = compute_sin2_thetaW(rho_uniform, modes)
    assert abs(s - 3/8) < 0.01, f"Expected 3/8=0.375, got {s:.4f}"


def test_boltzmann_normalization():
    """ρ* must be a probability distribution."""
    T = np.diag(np.random.exponential(1.0, 96))
    rho, _ = boltzmann_distribution(T)
    assert abs(rho.sum() - 1.0) < 1e-10
    assert np.all(rho >= 0)


def test_free_energy_minimum():
    """F_ρ[T] must be finite and real."""
    T = np.diag(np.abs(np.random.randn(96)))
    rho, _ = boltzmann_distribution(T)
    F = free_energy(rho, T)
    assert np.isfinite(F)
    assert np.isreal(F)


def test_rho_sum_to_one():
    """Distribution must sum to 1 after update."""
    T = np.diag(np.sort(np.abs(np.random.randn(96))))
    rho, _ = boltzmann_distribution(T)
    assert abs(rho.sum() - 1.0) < 1e-12


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
