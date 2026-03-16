import numpy as np
from scipy.integrate import solve_ivp
import json

# PDG 2024 Initial Conditions
m_Z = 91.1876
alpha_s_mZ = 0.1179
alpha_mZ = 1.0 / 127.95
sin2thW_mZ = 0.23122

# Couplings at m_Z
alpha_1_mZ = (5.0 / 3.0) * (alpha_mZ / (1.0 - sin2thW_mZ)) # GUT normalized
alpha_2_mZ = alpha_mZ / sin2thW_mZ
alpha_3_mZ = alpha_s_mZ

inv_a = np.array([1.0/alpha_1_mZ, 1.0/alpha_2_mZ, 1.0/alpha_3_mZ])

# Coefficients
b = np.array([41.0/10.0, -19.0/6.0, -7.0])
B = np.array([
    [199.0/50.0, 27.0/10.0, 44.0/5.0],
    [9.0/10.0,   35.0/6.0,  12.0],
    [11.0/10.0,  9.0/2.0,  -26.0]
])

def rge_2loop_inv(t, inv_a_vec):
    """
    d(1/alpha_i) / dt = - b_i / (2pi) - sum_j B_ij / (8pi^2 * alpha_j)
    """
    da_inv = np.zeros(3)
    for i in range(3):
        term1 = -b[i] / (2 * np.pi)
        # alpha_j = 1.0 / inv_a_vec[j]
        term2 = -sum(B[i, j] / inv_a_vec[j] for j in range(3)) / (8 * np.pi**2)
        da_inv[i] = term1 + term2
    return da_inv

t_start = 0.0
t_end = np.log(2e16 / m_Z)
t_eval = np.linspace(t_start, t_end, 500)

sol = solve_ivp(rge_2loop_inv, (t_start, t_end), inv_a, t_eval=t_eval, method='RK45', rtol=1e-11, atol=1e-13)

alpha1 = 1.0 / sol.y[0]
alpha2 = 1.0 / sol.y[1]
alpha3 = 1.0 / sol.y[2]

sin2_thW = (3.0/5.0 * alpha1) / (alpha2 + 3.0/5.0 * alpha1)

print(f"sin2thW(Lambda_GUT) at 2-loops: {sin2_thW[-1]:.5f}")

# Also calculate what 1-loop analytic would have done (and if it hits the pole)
# Analytic 1-loop: alpha_i(t) = alpha_i(0) / (1 - (b_i/(2pi)) * alpha_i(0) * t)
# To hit a pole, denominator <= 0 -> t = 2pi / (b_i * alpha_i(0))
t_pole = 2 * np.pi / (b[0] * alpha_1_mZ)
print(f"Pole t: {t_pole:.2f}, GUT t: {t_end:.2f}")

data = {
    "t": t_eval.tolist(),
    "alpha_1_inv": sol.y[0].tolist(),
    "alpha_2_inv": sol.y[1].tolist(),
    "alpha_3_inv": sol.y[2].tolist(),
    "sin2_thW": sin2_thW.tolist()
}

with open('/home/user/Desktop/data/tmp/rge_results.json', 'w') as f:
    json.dump(data, f)
