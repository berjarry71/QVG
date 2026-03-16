"""
mcmc_lattice.py
===============
Prototype for Lattice Non-Commutative Geometry via MCMC.

This script simulates the path integral Z = int dD exp(-S[D]) over random
Dirac matrix configurations, avoiding the perturbative block of QCD transition
by brute-force stochastic sampling of the geometric space.

Reference: Jarry, B. (2026). Quantum Vacuum Geometry.
"""

import numpy as np
import time

# --- Parameters ---
N_dim = 15          # Matrix rank (truncated fuzzy space approximation)
g_coupling = 0.5    # Coupling constant for Tr(D^4)
sweeps = 10000      # Number of Metropolis-Hastings sweeps
delta = 0.15        # Proposal step size

def spectral_action(D):
    """
    Computes S(D) = Tr(D^2) + (g/N) * Tr(D^4).
    This is the Seeley-DeWitt polynomial expansion of the full spectral action.
    """
    D2 = D @ D
    D4 = D2 @ D2
    return np.real(np.trace(D2) + (g_coupling / N_dim) * np.trace(D4))

def run_mcmc():
    print(f"Initializing random Dirac operator (Rank {N_dim})...")
    # Initial random Hermitian matrix
    A = np.random.randn(N_dim, N_dim) + 1j * np.random.randn(N_dim, N_dim)
    D = (A + A.conj().T) / 2.0
    
    current_S = spectral_action(D)
    accepted = 0
    
    start_time = time.time()
    print(f"Starting MCMC simulation for {sweeps} sweeps...")
    
    for i in range(sweeps):
        # Generate random Hermitian fluctuation
        H = np.random.randn(N_dim, N_dim) + 1j * np.random.randn(N_dim, N_dim)
        dH = delta * (H + H.conj().T) / 2.0
        D_new = D + dH
        
        new_S = spectral_action(D_new)
        dS = new_S - current_S
        
        # Metropolis accept/reject
        if dS < 0 or np.exp(-dS) > np.random.rand():
            D = D_new
            current_S = new_S
            accepted += 1

    end_time = time.time()
    
    print(f"\nSimulation finished in {end_time - start_time:.2f}s")
    print(f"Acceptance Rate: {accepted / sweeps * 100:.1f}%")
    print(f"Final Action S(D): {current_S:.4f}")
    
    # Analyze the spectrum
    evals = np.linalg.eigvalsh(D)
    print("\nPhysical Spectrum (First 5 eigenvalues):")
    print(np.round(evals[:5], 3))

if __name__ == "__main__":
    run_mcmc()
