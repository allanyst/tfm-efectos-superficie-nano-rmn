import numpy as np
from scipy import special
import math

def zeros_J0(p_max):
    return special.jn_zeros(0, p_max)  # zeros of J0 (array length p_max)

def gauss_nodes_weights(n):
    # returns nodes x_i in [-1,1] and weights w_i
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w

def scale_nodes_weights(x, w, a, b):
    # scale from [-1,1] to [a,b]
    nodes = 0.5*(b-a)*x + 0.5*(b+a)
    weights = 0.5*(b-a)*w
    return nodes, weights

def precompute_modes_gauss(R, L, d, D, m_max, p_max,
                           n_rho=40, n_z=40,
                           spectral_cutoff=None):
    """
    R, L, d in nm. D in nm^2/µs.
    n_rho, n_z: nodes for Gauss-Legendre (recommended 24..64)

    spectral_cutoff:
        float or None.
        If float: only modes with lambda^2 < spectral_cutoff are computed.
        If None: all modes up to (m_max, p_max) are computed.
    """

    # -------------------------------------------------
    # zeros and betas
    # -------------------------------------------------
    zeros = zeros_J0(p_max)        # zeros (dimensionless)
    betas = zeros / R              # beta in 1/nm

    J1_at_zeros = special.j1(zeros)
    denom = J1_at_zeros**2

    # -------------------------------------------------
    # Precompute spectral values (cheap)
    # -------------------------------------------------
    m_vals = np.arange(1, m_max+1)

    eta_e = (2.0 * m_vals * np.pi / L)
    eta_o = ((2.0*m_vals - 1.0) * np.pi / L)

    beta_matrix = np.tile(betas, (m_max, 1))
    eta_e_matrix = np.tile(eta_e[:, None], (1, p_max))
    eta_o_matrix = np.tile(eta_o[:, None], (1, p_max))

    lambda2_e = beta_matrix**2 + eta_e_matrix**2
    lambda2_o = beta_matrix**2 + eta_o_matrix**2

    # spectral masks
    if spectral_cutoff is None:
        mode_mask_e = np.ones_like(lambda2_e, dtype=bool)
        mode_mask_o = np.ones_like(lambda2_o, dtype=bool)
    else:
        mode_mask_e = (lambda2_e < spectral_cutoff)
        mode_mask_o = (lambda2_o < spectral_cutoff)

    # -------------------------------------------------
    # Gauss nodes/weights
    # -------------------------------------------------
    x_r, w_r = gauss_nodes_weights(n_rho)
    x_z, w_z = gauss_nodes_weights(n_z)

    rho_nodes, rho_weights = scale_nodes_weights(x_r, w_r, 0.0, R)
    z_nodes, z_weights = scale_nodes_weights(x_z, w_z, d, d+L)

    RHO = rho_nodes[:, None]
    Z = z_nodes[None, :]

    W_rz = W_rz_vec(RHO, Z)
    Wquad = (rho_weights[:, None]) * (z_weights[None, :])

    # -------------------------------------------------
    # Prepare output arrays
    # -------------------------------------------------
    I_sin = np.zeros((m_max, p_max), dtype=np.float64)
    I_cos = np.zeros((m_max, p_max), dtype=np.float64)

    # -------------------------------------------------
    # Precompute J0
    # -------------------------------------------------
    J0_vals = np.zeros((p_max, n_rho), dtype=np.float64)
    for p in range(p_max):
        J0_vals[p, :] = special.j0(betas[p] * rho_nodes)

    # -------------------------------------------------
    # Precompute axial functions
    # -------------------------------------------------
    sin_axial = np.zeros((m_max, n_z), dtype=np.float64)
    cos_axial = np.zeros((m_max, n_z), dtype=np.float64)

    for m in range(1, m_max+1):
        eta_e_m = 2.0 * m * math.pi / L
        eta_o_m = (2.0*m - 1.0) * math.pi / L
        sin_axial[m-1, :] = np.sin(eta_e_m * (z_nodes - d - L/2.0))
        cos_axial[m-1, :] = np.cos(eta_o_m * (z_nodes - d - L/2.0))

    # -------------------------------------------------
    # Compute integrals with spectral masking
    # -------------------------------------------------
    for p in range(p_max):
        J0_r = J0_vals[p, :]
        M = W_rz * (J0_r[:, None]) * Wquad

        for m in range(1, m_max+1):
            m_idx = m - 1

            # skip expensive calculation if mode not relevant
            if (not mode_mask_e[m_idx, p]) and (not mode_mask_o[m_idx, p]):
                continue

            val_sin = np.sum(M * sin_axial[m_idx, None, :])
            val_cos = np.sum(M * cos_axial[m_idx, None, :])

            I_sin[m_idx, p] = val_sin
            I_cos[m_idx, p] = val_cos

    # -------------------------------------------------
    # Pack results
    # -------------------------------------------------
    modes = {
        'I_sin': I_sin,
        'I_cos': I_cos,
        'denom': denom,
        'lambda2_e': lambda2_e,
        'lambda2_o': lambda2_o,
        'R': R,
        'L': L,
        'd': d,
        'D': D,
        'betas': betas,
        'zeros': zeros,
        'J1_at_zeros': J1_at_zeros,
        'mode_mask_e': mode_mask_e,
        'mode_mask_o': mode_mask_o
    }

    return modes

def W_rz_vec(rho, z):
    # rho, z can be numpy arrays (same shape broadcasting)
    denom = rho*rho + z*z
    # avoid divide by zero (set denom->inf) if necessary
    denom = np.where(denom == 0.0, np.inf, denom)
    return rho * (3.0 * z*z * denom**(-2.5) - denom**(-1.5))

def precompute_surface_integrals(R, L, d, modes, n_rho=80, n_z=80):
    """
    Calcula:
      - B_bottom[p] = ∫_0^R W(rho, d)      J0(beta_p rho) d rho
      - B_top[p]    = ∫_0^R W(rho, d+L)    J0(beta_p rho) d rho
      - Z_sin_wall[m] = ∫_d^{d+L} W(R, z) sin(eta_e (z - d - L/2)) dz
      - Z_cos_wall[m] = ∫_d^{d+L} W(R, z) cos(eta_o (z - d - L/2)) dz
    usando Gauss-Legendre.
    """
    betas = modes["betas"]     # shape (p_max,)
    m_max, p_max = modes["I_sin"].shape

    # Gauss en rho
    x_r, w_r = np.polynomial.legendre.leggauss(n_rho)
    rho_nodes = 0.5*(R - 0.0)*x_r + 0.5*(R + 0.0)
    rho_weights = 0.5*(R - 0.0)*w_r

    # Gauss en z
    x_z, w_z = np.polynomial.legendre.leggauss(n_z)
    z_nodes  = 0.5*((d+L) - d)*x_z + 0.5*((d+L) + d)
    z_weights = 0.5*((d+L) - d)*w_z

    W_bottom = W_rz_vec(rho_nodes, d)       # shape (n_rho,)
    W_top    = W_rz_vec(rho_nodes, d+L)     # shape (n_rho,)

    # W en la pared
    W_wall_z = W_rz_vec(np.full_like(z_nodes, R), z_nodes)  # shape (n_z,)

    # Precomputar J0(beta_p * rho) para todas las p
    J0_r = np.zeros((p_max, n_rho))
    for p in range(p_max):
        J0_r[p, :] = special.j0(betas[p] * rho_nodes)

    # B_bottom, B_top
    B_bottom = np.zeros(p_max)
    B_top    = np.zeros(p_max)
    for p in range(p_max):
        B_bottom[p] = np.sum(W_bottom * J0_r[p, :] * rho_weights)
        B_top[p]    = np.sum(W_top    * J0_r[p, :] * rho_weights)

    # Z_sin_wall, Z_cos_wall (solo dependen de m)
    Z_sin_wall = np.zeros(m_max)
    Z_cos_wall = np.zeros(m_max)
    for m in range(1, m_max+1):
        eta_e = 2.0 * m * math.pi / L
        eta_o = (2.0*m - 1.0) * math.pi / L
        phase = z_nodes - d - 0.5*L
        sin_ax = np.sin(eta_e * phase)
        cos_ax = np.cos(eta_o * phase)

        Z_sin_wall[m-1] = np.sum(W_wall_z * sin_ax * z_weights)
        Z_cos_wall[m-1] = np.sum(W_wall_z * cos_ax * z_weights)

    return {
        "B_bottom": B_bottom,
        "B_top": B_top,
        "Z_sin_wall": Z_sin_wall,
        "Z_cos_wall": Z_cos_wall,
    }

# # ---------- Pesos modales F_{m,p} de superficie ----------

def precompute_surface_modes(R, L, d, modes, surface_ints):
    """
    Construye F^{even}_{m,p} y F^{odd}_{m,p} tal como en la expresión analítica.
    """
    I_sin = modes["I_sin"]      # even
    I_cos = modes["I_cos"]      # odd
    betas = modes["betas"]
    lambda2_e = modes["lambda2_e"]
    lambda2_o = modes["lambda2_o"]

    # zeros = special.jn_zeros(0, I_sin.shape[1])
    # J1_at_zeros = special.j1(zeros)

    J1_at_zeros = modes["J1_at_zeros"]

    B_bottom = surface_ints["B_bottom"]
    B_top    = surface_ints["B_top"]
    Z_sin_wall = surface_ints["Z_sin_wall"]
    Z_cos_wall = surface_ints["Z_cos_wall"]

    m_max, p_max = I_sin.shape
    F_even = np.zeros((m_max, p_max))
    F_odd  = np.zeros((m_max, p_max))

    for m in range(1, m_max+1):
        mi = m - 1

        for p in range(p_max):
            beta = betas[p]
            lam2_e = lambda2_e[mi, p]
            lam2_o = lambda2_o[mi, p]
            J1 = J1_at_zeros[p]

            Is = I_sin[mi, p]   # even
            Ic = I_cos[mi, p]   # odd

            # Coeficientes axiales
            C_sin = ((-1)**m) * (2.0*m*math.pi) / (L * lam2_e)
            C_cos = ((-1)**m) * ((2.0*m - 1.0)*math.pi) / (L * lam2_o)

            B_delta = B_bottom[p] - B_top[p]

            # Caps
            F_caps_even = (B_delta / (J1*J1)) * (C_sin * Is)
            F_caps_odd  = (B_delta / (J1*J1)) * (C_cos * Ic)

            # Wall
            Zs = Z_sin_wall[mi]
            Zc = Z_cos_wall[mi]

            F_wall_even = (beta / J1) * (Is * Zs / lam2_e)
            F_wall_odd  = (beta / J1) * (Ic * Zc / lam2_o)

            F_even[mi, p] = F_caps_even + F_wall_even
            F_odd [mi, p] = F_caps_odd  + F_wall_odd

    return {
        "F_even": F_even,
        "F_odd":  F_odd,
        "lambda2_even": lambda2_e,
        "lambda2_odd":  lambda2_o,
    }

# ---------- G_surface(t) y G_surface(∞) ----------

from numba import njit, prange

@njit(parallel=True, fastmath=True)

def sum_surface_modes_numeric(F_even, lambda2_even, F_odd, lambda2_odd, D, times):
    m_max, p_max = F_even.shape
    n_times = times.shape[0]
    G = np.zeros(n_times, dtype=np.float64)
    for ti in prange(n_times):
        t = times[ti]
        acc = 0.0
        for mi in range(m_max):
            for pi in range(p_max):
                x_e = D * lambda2_even[mi, pi] * t
                x_o = D * lambda2_odd[mi, pi] * t
                # numéricamente estable: 1 - exp(-x) = -expm1(-x) (evita cancelación para x<<1)
                # math.expm1 puede no estar disponible en todas las versiones de numba; si no, usar 1-math.exp(-x)
                acc += F_even[mi, pi] * (-math.expm1(-x_e))
                acc += F_odd[mi, pi]  * (-math.expm1(-x_o))
        G[ti] = acc
    return G
    
def G_surface_vectorized(times_us, modes, surf_modes):
        times_us = np.asarray(times_us, dtype=np.float64)
        # F = surf_modes["F"]
        # lambda2_o = surf_modes["lambda2_o"]
        F_e = surf_modes["F_even"]
        lam2_e = surf_modes["lambda2_even"]
        F_o = surf_modes["F_odd"]
        lam2_o = surf_modes["lambda2_odd"]
        D = modes["D"]
        R = modes["R"]
        L = modes["L"]
        V = math.pi * R*R * L

        # Prefactor -8π^2 / V según tu ecuación de G^{surface}(t)
        # prefactor = -8.0 * math.pi**2 / V
        # return prefactor * sums
        sums = sum_surface_modes_numeric(F_e, lam2_e, F_o, lam2_o, D, times_us)
        prefactor = -8.0 * math.pi**2 / V
        return prefactor * sums


def G_surface_infty_from_modes(modes, surf_modes):
    """
    Medida armónica de superficie (límite t → ∞):

      G_∞ = - (8π^2 / V) * sum_{m,p} [ F^{even}_{m,p} + F^{odd}_{m,p} ]
    """
    F_even = surf_modes["F_even"]
    F_odd  = surf_modes["F_odd"]

    R = modes["R"]
    L = modes["L"]
    V = math.pi * R * R * L

    return -8.0 * math.pi**2 * (np.sum(F_even) + np.sum(F_odd)) / V