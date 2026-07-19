#
# Fast interpolation of LOTOS-EUROS data onto LES grid
# Uses the Numba library with JIT compilation to speed the calculations up.
# Arseni Doyennel (IHS EUR/VUA), Feb. 2025
#

import numpy as np
from numba import jit #for compiling the python 
from numba import prange   #for parallel execution
    
# --- Interpolation function ---
@jit(nopython=True, nogil=True)
def interpolate_vert_mass_conserving(conc_data, altitude_data, zt):

    """
    Interpolates concentrations from LOTOS-EUROS altitude levels to DALES vertical grid
    using a mass-conserving approach with smoother transitions and extrapolation.

    Args:
        conc_data (numpy.ndarray): concentrations, shape (z, lat, lon)
        altitude_data (numpy.ndarray): Corresponding altitude above sea level, shape (z, lat, lon)
        zt (numpy.ndarray): Target DALES vertical coordinates (height above ground)

    Returns:
        numpy.ndarray: Interpolated values, shape (lat, lon, zt)
    """
    nz, ny, nx = conc_data.shape
    nz_new = len(zt)  # Number of DALES target levels (center levels)

    interpolated_val = np.zeros((ny, nx, nz_new))

    # --- Build target edges from target centers (zt) ---
    z_edges_tgt = np.zeros(nz_new + 1)
    # interior edges = midpoints between centers
    for ii in range(1, nz_new):
        z_edges_tgt[ii] = 0.5 * (zt[ii - 1] + zt[ii])
    # bottom and top extrapolated edges
    z_edges_tgt[0] = zt[0] - (z_edges_tgt[1] - zt[0])
    z_edges_tgt[-1] = zt[-1] + (zt[-1] - z_edges_tgt[-2])

    # --- Loop over horizontal grid points ---
    for lat in range(ny):
        for lon in range(nx):
            altitude_profile = altitude_data[:, lat, lon]
            conc_profile = conc_data[:, lat, lon]

            # Estimate surface elevation (terrain height)
            if nz > 1:
                z_surface = altitude_profile[0] - (altitude_profile[1] - altitude_profile[0])
            else:
                z_surface = altitude_profile[0]

            # Shift altitudes so surface starts at 0 m
            altitude_profile = altitude_profile - z_surface

            # Enforce monotonic increase (safety against small numeric noise)
            for kk in range(1, nz):
                if altitude_profile[kk] < altitude_profile[kk - 1]:
                    altitude_profile[kk] = altitude_profile[kk - 1]

            # --- Build source layer edges from source centers ---
            z_edges_src = np.zeros(nz + 1)
            for ii in range(1, nz):
                z_edges_src[ii] = 0.5 * (altitude_profile[ii - 1] + altitude_profile[ii])
            z_edges_src[0] = altitude_profile[0] - (z_edges_src[1] - altitude_profile[0])
            z_edges_src[-1] = altitude_profile[-1] + (altitude_profile[-1] - z_edges_src[-2])

            # --- Interpolate mass conservatively ---
            for k in range(nz_new):
                z_low = z_edges_tgt[k]
                z_high = z_edges_tgt[k + 1]
                mass_sum = 0.0
                overlap_sum = 0.0

                for j in range(nz):
                    z1 = z_edges_src[j]
                    z2 = z_edges_src[j + 1]
                    # overlap region
                    overlap_low = z_low if z_low > z1 else z1
                    overlap_high = z_high if z_high < z2 else z2
                    overlap = overlap_high - overlap_low
                    if overlap > 0.0:
                        mass_sum += conc_profile[j] * overlap
                        overlap_sum += overlap

                if overlap_sum > 0.0:
                    interpolated_val[lat, lon, k] = mass_sum / overlap_sum
                else:
                    # Extrapolate with nearest layer (physically consistent)
                    if z_high <= altitude_profile[0]:
                        interpolated_val[lat, lon, k] = conc_profile[0]
                    elif z_low >= altitude_profile[-1]:
                        interpolated_val[lat, lon, k] = conc_profile[-1]
                    else:
                        interpolated_val[lat, lon, k] = 0.0

    return interpolated_val
    
@jit(nopython=True, nogil=True)
def interpolate_vert_mass_conserving_CAMS(conc_data, altitude_data, zt):
    """
    Interpolates CAMS concentrations (given at model levels) to LES vertical grid
    using a mass-conserving approach.

    Args:
        conc_data (numpy.ndarray): CAMS concentrations, shape (z, lat, lon)
        altitude_data (numpy.ndarray): CAMS altitude of model levels (m above ellipsoid), shape (z, lat, lon)
        zt (numpy.ndarray): LES target vertical coordinates (m above ground or ellipsoid)

    Returns:
        numpy.ndarray: Interpolated concentrations, shape (lat, lon, zt)
    """
    nz, ny, nx = conc_data.shape
    nz_new = len(zt)
    interpolated_val = np.zeros((ny, nx, nz_new))

    # --- Target edges ---
    z_edges_tgt = np.zeros(nz_new + 1)
    for ii in range(1, nz_new):
        z_edges_tgt[ii] = 0.5 * (zt[ii - 1] + zt[ii])
    z_edges_tgt[0] = zt[0] - (z_edges_tgt[1] - zt[0])
    z_edges_tgt[-1] = zt[-1] + (zt[-1] - z_edges_tgt[-2])

    for lat in range(ny):
        for lon in range(nx):
            altitude_profile = altitude_data[:, lat, lon]
            conc_profile = conc_data[:, lat, lon]

            # Skip invalid profiles
            if np.isnan(altitude_profile).any() or np.isnan(conc_profile).any():
                continue

            # Ensure monotonic altitude
            for kk in range(1, nz):
                if altitude_profile[kk] < altitude_profile[kk - 1]:
                    altitude_profile[kk] = altitude_profile[kk - 1]

            # --- Source edges ---
            z_edges_src = np.zeros(nz + 1)
            for ii in range(1, nz):
                z_edges_src[ii] = 0.5 * (altitude_profile[ii - 1] + altitude_profile[ii])
            z_edges_src[0] = altitude_profile[0] - (z_edges_src[1] - altitude_profile[0])
            z_edges_src[-1] = altitude_profile[-1] + (altitude_profile[-1] - z_edges_src[-2])

            # --- Conservative interpolation ---
            for k in range(nz_new):
                z_low, z_high = z_edges_tgt[k], z_edges_tgt[k + 1]
                mass_sum, overlap_sum = 0.0, 0.0

                for j in range(nz):
                    z1, z2 = z_edges_src[j], z_edges_src[j + 1]
                    overlap_low = max(z_low, z1)
                    overlap_high = min(z_high, z2)
                    overlap = overlap_high - overlap_low
                    if overlap > 0.0:
                        mass_sum += conc_profile[j] * overlap
                        overlap_sum += overlap

                if overlap_sum > 0.0:
                    interpolated_val[lat, lon, k] = mass_sum / overlap_sum
                else:
                    if z_high <= altitude_profile[0]:
                        interpolated_val[lat, lon, k] = conc_profile[0]
                    elif z_low >= altitude_profile[-1]:
                        interpolated_val[lat, lon, k] = conc_profile[-1]
                    else:
                        interpolated_val[lat, lon, k] = 0.0

    return interpolated_val

@jit(nopython=True, nogil=True)
def calc_horz_interpolation_factors_unaligned(i0, j0, fi, fj, x, y, x_LS, y_LS, x_sw, y_sw):
    epsilon = 1e-10
    for i in range(x.shape[0]):
        for j in range(y.shape[1]):
            # Adjust coordinates relative to SW corner (turned off as we interpolate in geo lat and lon)
            x_rel = x[i, j] #- x_sw
            y_rel = y[i, j] #- y_sw

            # Find nearest indices in LS grid
            i0[i, j] = np.where(x_LS - x_rel <= 0)[0][-1]
            j0[i, j] = np.where(y_LS - y_rel <= 0)[0][-1]

            # Compute interpolation factors
            x_diff = x_LS[i0[i, j] + 1] - x_LS[i0[i, j]] if i0[i, j] + 1 < len(x_LS) else 1.0
            y_diff = y_LS[j0[i, j] + 1] - y_LS[j0[i, j]] if j0[i, j] + 1 < len(y_LS) else 1.0

            fi[i, j] = 1. - ((x_rel - x_LS[i0[i, j]]) / (x_diff + epsilon))
            fj[i, j] = 1. - ((y_rel - y_LS[j0[i, j]]) / (y_diff + epsilon))


@jit(nopython=True, nogil=True)
def interpolate_kernel_2d_unaligned(field_LES, field_LS, i0, j0, ifac, jfac):
    itot, jtot, ktot = field_LES.shape
    for i in range(itot):
        for j in range(jtot):
            for k in range(ktot):
                il, jl = i0[i, j], j0[i, j]
                il_next = il + 1 if il + 1 < field_LS.shape[0] else il
                jl_next = jl + 1 if jl + 1 < field_LS.shape[1] else jl
                
                field_LES[i, j, k] = (
                    jfac[i, j] * (
                        ifac[i, j] * field_LS[il, jl, k] +
                        (1 - ifac[i, j]) * field_LS[il_next, jl, k]
                    ) +
                    (1 - jfac[i, j]) * (
                        ifac[i, j] * field_LS[il, jl_next, k] +
                        (1 - ifac[i, j]) * field_LS[il_next, jl_next, k]
                    )
                )
                
                
def ensure_monotonic(arr, name, flip=False):
    """Ensure that the coordinate array is monotonically increasing."""
    if np.all(np.diff(arr) > 0):
        return arr
    else:
        print(f"Warning: {name} is not increasing, flipping the order.")
        arr = arr[::-1]
    
    # Optionally flip based on the flip flag
    if flip:
        arr = arr[::-1]
    
    return arr

class GridInterpolatorLE:
    def __init__(self, x_LS, y_LS, z_LS, x, y, z, x_sw, y_sw):
        self.x_LS = ensure_monotonic(x_LS, "x_LS")
        self.y_LS = ensure_monotonic(y_LS, "y_LS")
        
        self.i0 = np.zeros_like(x, dtype=int)
        self.ifac = np.zeros_like(x, dtype=float)
        self.j0 = np.zeros_like(y, dtype=int)
        self.jfac = np.zeros_like(y, dtype=float)
        
        calc_horz_interpolation_factors_unaligned(
            self.i0, self.j0, self.ifac, self.jfac, x, y, self.x_LS, self.y_LS, x_sw, y_sw
        )

    def interpolate_2d(self, field_LS):
        field_LES = np.zeros((self.i0.shape[0], self.j0.shape[1], field_LS.shape[2]), dtype=field_LS.dtype)
        interpolate_kernel_2d_unaligned(field_LES, field_LS, self.i0, self.j0, self.ifac, self.jfac)
        return field_LES

