# Simple petrophysical log analysis
# Reads LAS (or CSV) and computes:
# - density porosity
# - Vsh from GR
# - effective porosity (shale-corrected)
# - Archie Sw from Rt

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Optional: import lasio if you have LAS files
try:
    import lasio
    HAS_LASIO = True
except Exception:
    HAS_LASIO = False

# PARAMETERS (tune for your data)
RHO_MATRIX = 2.65   # g/cc (typical for clean sandstone)
RHO_FLUID = 1.0     # g/cc (water)
GR_MIN = 10.0       # API (clean sand baseline)
GR_MAX = 120.0      # API (shale baseline)
ARCHIE_A = 1.0
ARCHIE_M = 2.0
ARCHIE_N = 2.0
RW = 0.1            # ohm-m (formation water resistivity)

def read_logs_las_or_csv(path):
    """Return DataFrame indexed by depth with columns: GR, RHOB, Rt (optional)."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.las' and HAS_LASIO:
        las = lasio.read(path)
        df = las.df()  # index is depth by default
        # normalize typical mnemonics
        mapping = {}
        for col in df.columns:
            ll = col.lower()
            if 'gr' in ll or 'gamma' in ll:
                mapping[col] = 'GR'
            elif 'rho' in ll or 'rhob' in ll:
                mapping[col] = 'RHOB'
            elif 'res' in ll or 'rt' in ll or 'resist' in ll:
                mapping[col] = 'Rt'
        df = df.rename(columns=mapping)
        return df
    else:
        # expect CSV with a depth column and standard mnemonics
        df = pd.read_csv(path)
        if 'DEPTH' in df.columns:
            df = df.set_index('DEPTH')
        elif 'Depth' in df.columns:
            df = df.set_index('Depth')
        # ensure columns upper-cased
        df.columns = [c.strip() for c in df.columns]
        return df

def compute_density_porosity(rhob, rho_matrix=RHO_MATRIX, rho_fluid=RHO_FLUID):
    # phi = (rho_matrix - rho_bulk) / (rho_matrix - rho_fluid)
    phi = (rho_matrix - rhob) / (rho_matrix - rho_fluid)
    return np.clip(phi, 0, 1)

def compute_vsh(gr, gr_min=GR_MIN, gr_max=GR_MAX):
    vsh = (gr - gr_min) / (gr_max - gr_min)
    return np.clip(vsh, 0, 1)

def archie_sw(rt, phi_e, a=ARCHIE_A, m=ARCHIE_M, n=ARCHIE_N, rw=RW):
    # Sw = ((a * Rw) / (phi_e**m * Rt)) ** (1/n)
    # Avoid divide-by-zero and negative phi_e
    phi_e_safe = np.maximum(phi_e, 1e-6)
    rt_safe = np.maximum(rt, 1e-6)
    sw = ((a * rw) / (phi_e_safe**m * rt_safe)) ** (1.0 / n)
    return np.clip(sw, 0, 1)

def main(path):
    df = read_logs_las_or_csv(path)
    # Check required logs
    if 'RHOB' not in df.columns:
        raise ValueError("RHOB (bulk density) log not found in file.")
    if 'GR' not in df.columns:
        print("Warning: GR log not found. Vsh will be set to 0.")
    # compute
    df['PHI_DENS'] = compute_density_porosity(df['RHOB'])
    if 'GR' in df.columns:
        df['VSH'] = compute_vsh(df['GR'])
    else:
        df['VSH'] = 0.0
    df['PHI_E'] = df['PHI_DENS'] * (1 - df['VSH'])  # simple shale correction
    if 'Rt' in df.columns:
        df['Sw'] = archie_sw(df['Rt'], df['PHI_E'])
    else:
        df['Sw'] = np.nan
    # simple summary
    print("Summary statistics:")
    print(df[['PHI_DENS', 'VSH', 'PHI_E', 'Sw']].describe())

    # plot logs
    fig, axes = plt.subplots(ncols=4, nrows=1, figsize=(10, 8), sharey=True)
    depth = df.index
    axes[0].plot(df.get('GR', np.zeros_like(depth)), depth, 'g')
    axes[0].set_xlabel('GR')
    axes[1].plot(df['RHOB'], depth, 'k')
    axes[1].set_xlabel('RHOB')
    axes[2].plot(df['PHI_DENS'], depth, 'b', label='Phi_dens')
    axes[2].plot(df['PHI_E'], depth, 'c', label='Phi_eff')
    axes[2].set_xlabel('Porosity')
    axes[2].legend()
    if 'Rt' in df.columns:
        axes[3].semilogx(df['Rt'], depth, 'r')
        axes[3].set_xlabel('Rt (ohm-m)')
        ax4b = axes[3].twiny()
        ax4b.plot(df['Sw'], depth, 'm')
        ax4b.set_xlabel('Sw')
    else:
        axes[3].text(0.5, 0.5, 'No Rt log', transform=axes[3].transAxes, ha='center')
        axes[3].set_xlabel('Rt / Sw')

    axes[0].invert_yaxis()
    axes[0].set_ylabel('Depth')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: python petrophysics_example.py <path_to_las_or_csv>")
        sys.exit(1)
    main(sys.argv[1])
