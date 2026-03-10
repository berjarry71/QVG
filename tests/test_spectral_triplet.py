"""Tests for spectral_triplet.py"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'core'))
import numpy as np
import pytest
from spectral_triplet import (build_fermion_modes, compute_sin2_thetaW,
                               verify_connes_identity)

def test_N_96():
    assert len(build_fermion_modes()) == 96

def test_connes_chamseddine_identity():
    """Tr(Y²)/Tr(T₃²) = 5/3 — structural property of A_F."""
    modes = build_fermion_modes()
    ratio, expected, ok = verify_connes_identity(modes)
    assert ok, f"Connes identity failed: {ratio:.6f} ≠ {expected:.6f}"

def test_sin2_uniform():
    """Uniform ρ → sin²θ_W = 3/8 (Connes-Chamseddine 1997)."""
    modes = build_fermion_modes()
    rho = np.ones(96) / 96
    s = compute_sin2_thetaW(rho, modes)
    assert abs(s - 3/8) < 1e-6, f"Expected 0.375, got {s:.6f}"

def test_symmetry_particle_antiparticle():
    """Particles and antiparticles contribute equally to Tr(Y²)."""
    modes = build_fermion_modes()
    Y2_part  = sum(m['Y']**2 for m in modes if m['ptype']=='particle')
    Y2_anti  = sum(m['Y']**2 for m in modes if m['ptype']=='antiparticle')
    assert abs(Y2_part - Y2_anti) < 1e-10

def test_rho_normalization():
    """ρ must sum to 1."""
    rho = np.ones(96) / 96
    assert abs(rho.sum() - 1.0) < 1e-12

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
