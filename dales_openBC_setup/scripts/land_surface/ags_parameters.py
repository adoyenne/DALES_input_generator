import math
import warnings

def compute_dmax(gm25, a, b):
    """
    Compute Dmax using:
    gm25 in m/s → converted to mm/s internally
    """
    Dmax_gkg = math.exp((a - math.log(gm25 * 1000.0)) / b) #g/kg
    
    Dmax_kgkg = Dmax_gkg / 1000.0
    return Dmax_kgkg


# ============================================================
# Vegetation lookup table
# ============================================================

ags_params   = ['gm25','Ammax25','f0','alpha0','co2_comp298','T1gm','T1Am'] #,'Dmax'

VEG_LOOKUP = {

    # ===== Forests =====
    "fce": {               #the temperate forest (Loobos values) evergreen needleleaf forest (C3)
        "gm25": 3.e-3,
        "Dmax": 0.124,
        "Ammax25": 2.2,
        "f0": 0.89,
        "a": 0.0,
        "b": 0.0,
        "T1gm": 278.0,
        "T1Am": 281.0,
        "co2_comp298": 82.2,
        "alpha0": 0.017,
    },
    
    "fbd": {
        "gm25": 2.0e-3,
        "Dmax": 0.109,        # (broadleaf slightly more sensitive than fce)
        "Ammax25": 1.83,    
        "f0": 0.80,          
        "a": 0.0,
        "b": 0.0,
        "T1gm": 278.0,
        "T1Am": 281.0,
        "co2_comp298": 42.0,  # (C3 broadleaf)
        "alpha0": 0.0142,    
    },
    
    # ===== Low vegetation (C3 default) =====
    "ara": {
        "gm25": 1.3e-3,
        "Ammax25": 2.20,
        "f0": 0.85,
        "alpha0": 0.0142,
        "a": 2.381,
        "b": 0.6103,
        "T1gm": 278.0,
        "T1Am": 281.0,
        "co2_comp298": 42.0,
    },

    "crp": {
        "gm25": 1.4e-3,
        "Ammax25": 1.83,
        "f0": 0.92,
        "alpha0": 0.0142,
        "a": 2.381,
        "b": 0.6103,
        "T1gm": 278.0,
        "T1Am": 281.0,
        "co2_comp298": 42.0,
    },

    "sem": {
        "gm25": 1.0e-3,
        "Ammax25": 1.83,
        "f0": 0.80,
        "alpha0": 0.0142,
        "a": 2.381,
        "b": 0.6103,
        "T1gm": 278.0,
        "T1Am": 281.0,
        "co2_comp298": 42.0,
    },

    # ===== Grass (special handling: C3 vs C4) =====
    "grs": {
        "C3": {
            "gm25": 1.5e-3,
            "Ammax25": 2.0,
            "f0": 0.75,
            "alpha0": 0.0142,
            "a": 2.381,
            "b": 0.6103,
            "T1gm": 278.0,
            "T1Am": 281.0,
            "co2_comp298": 42.0,
        },
        "C4": {
            "gm25": 2.3e-3,
            "Ammax25": 1.83,
            "f0": 0.70,
            "alpha0": 0.0117,
            "a": 5.323,
            "b": 0.8923,
            "T1gm": 286.15,
            "T1Am": 286.15,
            "co2_comp298": 2.6,
        }
    }
}


# ============================================================
# Public API
# ============================================================
    
def get_veg_params(lu, planttype=None, only_ags=False):

    """
    Retrieve vegetation parameters.

    Parameters
    ----------
    lu : str
        Land-use short name ('fce', 'grs', etc.)
    planttype : int, optional
        If lu == 'grs':
            planttype == 4 → C4
            else → C3

    Returns
    -------
    dict
        Dictionary with vegetation parameters
    """

    lu = lu.lower()

    # --- handle non-vegetation safely ---
    if lu not in VEG_LOOKUP:
        warnings.warn(f"Unknown land-use type (treated as non-vegetation): {lu}")
        return {k: 0.0 for k in ags_params}   # SAFE fallback

    # --- Grass special case ---
    if lu == "grs":
        key = "C4" if planttype == 4 else "C3"
        params = VEG_LOOKUP["grs"][key].copy()
    else:
        params = VEG_LOOKUP[lu].copy()

    # --- Compute Dmax if needed ---
    #if "Dmax" not in params or params["Dmax"] in [None, "auto"]:
        #params["Dmax"] = compute_dmax(
        #    params["gm25"],
        #    params["a"],
        #    params["b"]
        #)

    # --- AGS filter ---
    if only_ags:
        return {k: params.get(k, 0.0) for k in ags_params}

    return params