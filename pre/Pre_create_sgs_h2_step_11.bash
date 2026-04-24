#!/bin/bash +x
################################################################################
# PASSO 0 - MOdulos Cray a carregar
################################################################################
#
# Adicione em seu ~/.bashrc ou script de ambiente:
#
#   module purge
#   module load PrgEnv-gnu
module unload PrgEnv-cray/8.6.0
module load   PrgEnv-gnu/8.6.0
#   module load cray-netcdf        # NetCDF C + Fortran
module load cray-netcdf-hdf5parallel/4.9.0.15
#   module load cray-hdf5          # HDF5
module load cray-hdf5-parallel/1.14.3.3
#   module load cray-mpich         # MPI
module load cray-mpich/8.1.31
#   module load cray-pio           # Parallel I/O (necessario para MOM6)
#module load pio/2.6.5
#   module load esmf/8.4.0         # ESMF (ajuste versao disponivel no seu Cray)
module load cray-pals/1.6.1
module load cray-python/3.11.7
module load autoconf/2.72
module load libfabric/1.22.0
module load nco/5.3.6


module list
export FC=ftn
export CC=cc
export PATH=/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/automake/automake-1.16.5/bin:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf/bin:\
/p/projetos/monan_atm/paulo.kubota/coupler/lib_couplers/fms/FRE-NCtools/build/bin:$PATH
export LD_LIBRARY_PATH=/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/libssl/lib:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/automake/automake-1.16.5/lib:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf/lib:\
/p/projetos/monan_atm/paulo.kubota/coupler/lib_couplers/fms/FRE-NCtools/build/lib:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/lib:$LD_LIBRARY_PATH
path_fre=/p/projetos/monan_atm/paulo.kubota/coupler/lib_couplers/fms/FRE-NCtools/build/bin
rm ./INPUT/sgs_h2.nc
cat << 'EOF'>create_sgs_h2.py
"""
create_sgs_h2.py
-----------------
Cria ./INPUT/sgs_h2.nc (rugosidade topografica sub-grade h^2) para o MOM6.
Resolucao generica.

Uso rapido:
    python create_sgs_h2.py --preset one --synthetic
    python create_sgs_h2.py --preset quarter --raw-file gebco_2023.nc

Dependencias:
    pip install numpy netCDF4 scipy
"""

import argparse
import sys
import numpy as np
from pathlib import Path
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator

try:
    import netCDF4 as nc
except ImportError:
    sys.exit("Instale: pip install netCDF4 scipy")

from Mom6_grid_config import (GridConfig, add_grid_args, add_vertical_args,
                               print_presets, KAPPA_ITIDES, KAPPA_H2_FACTOR)


H2_MAX = 4000.0**2   # [m^2] - amplitude maxima 4 km
H2_MIN = 1.0         # [m^2]


def read_highres_topo(filepath):
    ds  = nc.Dataset(filepath)
    lon = next((ds.variables[v] for v in ["lon","longitude","x"] if v in ds.variables), None)
    lat = next((ds.variables[v] for v in ["lat","latitude","y"] if v in ds.variables), None)
    dep = next((ds.variables[v] for v in ["elevation","z","depth","topo"] if v in ds.variables), None)
    if not all([lon, lat, dep]):
        ds.close()
        raise ValueError(f"Variaveis nao encontradas: {list(ds.variables.keys())}")
    lon  = np.array(lon[:]); lat = np.array(lat[:])
    elev = np.array(dep[:], dtype=np.float64)
    ds.close()
    print(f"  Alta resolucao: {elev.shape}")
    return lon, lat, elev


def compute_h2_block(lon_hr, lat_hr, elev_hr, cfg, lons_m, lats_m):
    """
    Calcula h^2 por bloco: para cada celula do modelo,
    h^2 = variancia da topografia oceanica sub-grade.
    Metodo vetorizado em latitude.
    """
    h2 = np.zeros((cfg.njglobal, cfg.niglobal))
    dlat_h = cfg.dlat / 2.0
    dlon_h = cfg.dlon / 2.0

    if lon_hr[0] > lon_hr[-1]:
        lon_hr = lon_hr[::-1]; elev_hr = elev_hr[:, ::-1]
    if lat_hr[0] > lat_hr[-1]:
        lat_hr = lat_hr[::-1]; elev_hr = elev_hr[::-1, :]

    for j, lat_c in enumerate(lats_m):
        jj = np.where((lat_hr >= lat_c - dlat_h) &
                      (lat_hr <  lat_c + dlat_h))[0]
        if len(jj) == 0:
            continue
        for i, lon_c in enumerate(lons_m):
            ii = np.where((lon_hr >= lon_c - dlon_h) &
                          (lon_hr <  lon_c + dlon_h))[0]
            if len(ii) == 0:
                continue
            sub = elev_hr[np.ix_(jj, ii)]
            pts = sub[sub < 0]
            if len(pts) >= 4:
                h2[j, i] = np.var(pts)

        if j % 50 == 0:
            print(f"  j={j}/{cfg.njglobal}  ({100*j/cfg.njglobal:.0f}%)")

    return h2


def synthetic_h2(cfg, lons_m, lats_m):
    print("  Gerando h^2 sintetico...")
    lo, la = np.meshgrid(lons_m, lats_m)

    h2 = np.random.exponential(500.0, (cfg.njglobal, cfg.niglobal))**2

    # Dorsais mesoceanicas - alta rugosidade
    for lon_r, amp in [(-25, 8000), (-110, 3000), (70, 5000)]:
        h2 += amp**2 * np.exp(-((lo - lon_r)**2) / (3.0**2))

    # Terra = 0
    land = (
        (la >  60) | (la < -70) |
        ((lo > -10) & (lo <  50) & (la > 35)) |
        ((lo > -30) & (lo <  60) & (la > -10) & (la < 35)) |
        ((lo >  60) & (lo < 150) & (la >  5)) |
        ((lo > -120) & (lo < -60) & (la > 10)) |
        ((lo > -80) & (lo < -35) & (la > -55) & (la < 10))
    )
    h2[land] = 0.0

    # Plataformas - baixa rugosidade
    shelf = (la > 55) & ~land
    h2[shelf] *= 0.1

    return gaussian_filter(h2, sigma=1.0)


def apply_limits(h2):
    ocean = h2 > 0
    h2 = np.where(ocean & (h2 < H2_MIN), H2_MIN, h2)
    h2 = np.where(ocean & (h2 > H2_MAX), H2_MAX, h2)
    return h2


def write_sgs_h2(cfg, h2, lons_m, lats_m, outpath):
    ds = nc.Dataset(outpath, "w", format="NETCDF4_CLASSIC")
    ds.createDimension("nx", cfg.niglobal)
    ds.createDimension("ny", cfg.njglobal)

    vl = ds.createVariable("nx", "f8", ("nx",))
    vl.units = "degrees_east"; vl[:] = lons_m
    vl = ds.createVariable("ny", "f8", ("ny",))
    vl.units = "degrees_north"; vl[:] = lats_m

    vh2 = ds.createVariable("h2", "f4", ("ny","nx"), fill_value=0.0)
    vh2.units            = "m2"
    vh2.long_name        = "Sub-grid topographic height variance"
    vh2.KAPPA_ITIDES     = KAPPA_ITIDES
    vh2.KAPPA_H2_FACTOR  = KAPPA_H2_FACTOR
    vh2[:]               = h2.astype(np.float32)

    ds.title    = (f"MOM6 sgs_h2 {cfg.niglobal}X{cfg.njglobal} "
                   f"({cfg.dlon:.4g} deg)")
    ds.NIGLOBAL = cfg.niglobal
    ds.NJGLOBAL = cfg.njglobal
    ds.history  = "Created by create_sgs_h2.py"
    ds.close()

    ocean = h2 > 0
    print(f"[OK] {outpath}")
    print(f"     h^2: [{h2[ocean].min():.1f}, {h2[ocean].max():.2e}] m^2")
    print(f"     h RMS: [{np.sqrt(h2[ocean].min()):.1f}, "
          f"{np.sqrt(h2[ocean].max()):.0f}] m")


def main():
    parser = argparse.ArgumentParser(
        description="Cria ./INPUT/sgs_h2.nc para o MOM6 - resolucao generica",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--list-presets", action="store_true")
    parser.add_argument("--raw-file",  help="Batimetria alta resolucao (GEBCO 30-arcsec)")
    parser.add_argument("--synthetic", action="store_true")
    parser.add_argument("--output-dir", default="INPUT")

    add_grid_args(parser)
    args = parser.parse_args()

    if args.list_presets:
        print_presets(); sys.exit(0)
    if not args.raw_file and not args.synthetic:
        parser.error("Informe --raw-file ou --synthetic")

    cfg = GridConfig(args)
    print("\nCriando ./INPUT/sgs_h2.nc")
    print(cfg.summary())

    lons_m, lats_m = cfg.model_grid()

    if args.synthetic:
        h2 = synthetic_h2(cfg, lons_m, lats_m)
    else:
        print(f"\nLendo {args.raw_file}...")
        lon_hr, lat_hr, elev_hr = read_highres_topo(args.raw_file)
        print("Calculando variancia sub-grade (pode demorar)...")
        h2 = compute_h2_block(lon_hr, lat_hr, elev_hr, cfg, lons_m, lats_m)

    print("Aplicando limites...")
    h2 = apply_limits(h2)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_sgs_h2(cfg, h2, lons_m, lats_m, outdir / "sgs_h2.nc")


if __name__ == "__main__":
    main()
EOF
python3 create_sgs_h2.py --preset one1  --raw-file GEBCO/GEBCO_2025_sub_ice.nc
rm create_sgs_h2.py
