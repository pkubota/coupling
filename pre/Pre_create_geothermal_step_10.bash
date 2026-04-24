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
rm ./INPUT/geothermal_heating_cm2g.nc
cat << 'EOF'>create_geothermal.py
"""
create_geothermal.py
---------------------
Cria ./INPUT/geothermal_heating_cm2g.nc (fluxo de calor geotermico) para o MOM6.
Resolucao generica.

Uso rapido:
    python create_geothermal.py --preset one --synthetic --topog-file INPUT/topog.nc
    python create_geothermal.py --preset quarter --raw-file global_hf.nc \\
                                --topog-file INPUT/topog.nc

Dependencias:
    pip install numpy netCDF4 scipy
"""

import argparse
import sys
import numpy as np
from pathlib import Path
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

try:
    import netCDF4 as nc
except ImportError:
    sys.exit("Instale: pip install netCDF4 scipy")

from Mom6_grid_config import (GridConfig, add_grid_args, add_vertical_args,
                               print_presets, GEOTHERMAL_SCALE)


MEAN_HF = 0.086   # W m^-2
MAX_HF  = 10.0    # W m^-2
MIN_HF  = 0.03    # W m^-2


def read_topog_mask(topog_file, cfg):
    """Le mascara oceanica do ./INPUT/topog.nc, verificando consistencia de grade."""
    ds    = nc.Dataset(topog_file)
    depth = np.array(ds.variables["depth"][:])

    ny_f, nx_f = depth.shape
    ds.close()

    if ny_f != cfg.njglobal or nx_f != cfg.niglobal:
        raise ValueError(
            f"./INPUT/topog.nc tem shape ({ny_f}, {nx_f}), "
            f"mas a grade configurada e ({cfg.njglobal}, {cfg.niglobal}).\n"
            f"Use o mesmo --preset / --niglobal / --njglobal que usou em create_topog.py."
        )
    ocean = depth > 0.0
    print(f"  Mascara: {ocean.sum()} pontos oceanicos ({100*ocean.mean():.1f}%)")
    return ocean


def read_raw_hf(filepath):
    ds  = nc.Dataset(filepath)
    lon = next((ds.variables[v] for v in ["lon","longitude","x"] if v in ds.variables), None)
    lat = next((ds.variables[v] for v in ["lat","latitude","y"] if v in ds.variables), None)
    hf  = next((ds.variables[v] for v in ["heat_flow","hf","q","geo_heat","z"]
                if v in ds.variables), None)
    if not all([lon, lat, hf]):
        ds.close()
        raise ValueError(f"Variaveis nao encontradas: {list(ds.variables.keys())}")
    lon = np.array(lon[:]); lat = np.array(lat[:])
    hf  = np.array(hf[:],  dtype=np.float64)
    ds.close()
    if hf.mean() > 10.0:   # detecta mW m^-2
        print("  Convertendo mW m^-2 - W m^-2")
        hf /= 1000.0
    print(f"  Dado bruto: {hf.shape}  [{hf.min():.4f},{hf.max():.3f}] W m^-2")
    return lon, lat, hf


def interpolate_hf(lon_r, lat_r, hf_r, lons_m, lats_m, nj, ni):
    if lon_r[0] > lon_r[-1]:
        lon_r = lon_r[::-1]; hf_r = hf_r[:, ::-1]
    if lat_r[0] > lat_r[-1]:
        lat_r = lat_r[::-1]; hf_r = hf_r[::-1, :]
    interp = RegularGridInterpolator(
        (lat_r, lon_r), hf_r,
        method="linear", bounds_error=False, fill_value=MEAN_HF
    )
    lo2d, la2d = np.meshgrid(lons_m, lats_m)
    return interp(np.column_stack([la2d.ravel(), lo2d.ravel()])).reshape(nj, ni)


def synthetic_hf(cfg, lons_m, lats_m):
    """Fluxo geotermico sintetico baseado em caracteristicas tectonicas."""
    print("  Gerando campo geotermico sintetico...")
    lo, la = np.meshgrid(lons_m, lats_m)
    hf     = np.ones((cfg.njglobal, cfg.niglobal)) * MEAN_HF

    # Dorsais mesoceanicas
    for lon_r, amp in [(-20, 0.40), (-110, 0.35), (70, 0.30), (60, 0.25)]:
        hf += amp * np.exp(-((lo - lon_r)**2) / (8.0**2))

    # Pontos quentes
    for lon_h, lat_h, amp in [(-156, 20, 0.80), (-18, 65, 0.60), (70, -50, 0.40)]:
        hf += amp * np.exp(-((lo-lon_h)**2+(la-lat_h)**2) / (2.0**2))

    # Oceano antigo (Pacifico NW) - mais frio
    old = ((lo > 140.0) | (lo < -180.0)) & (la > 10) & (la < 45)
    hf[old] *= 0.6

    hf = gaussian_filter(hf, sigma=max(1.0, 2.0 / cfg.dlon))
    return np.maximum(hf, MIN_HF)


def apply_mask_and_limits(hf, ocean_mask):
    hf[~ocean_mask] = 0.0
    ocean = ocean_mask & (hf > 0)
    hf = np.where(ocean & (hf < MIN_HF), MIN_HF, hf)
    hf = np.where(ocean & (hf > MAX_HF), MAX_HF, hf)
    return hf


def write_geothermal(cfg, hf, lons_m, lats_m, outpath):
    ds = nc.Dataset(outpath, "w", format="NETCDF4_CLASSIC")
    ds.createDimension("nx", cfg.niglobal)
    ds.createDimension("ny", cfg.njglobal)

    vl = ds.createVariable("nx", "f8", ("nx",)); vl.units = "degrees_east";  vl[:] = lons_m
    vl = ds.createVariable("ny", "f8", ("ny",)); vl.units = "degrees_north"; vl[:] = lats_m

    vhf = ds.createVariable("geo_heat", "f4", ("ny","nx"), fill_value=0.0)
    vhf.units              = "W m-2"
    vhf.long_name          = "Geothermal heat flux"
    vhf.standard_name      = "upward_geothermal_heat_flux_at_sea_floor"
    vhf.positive           = "upward"
    vhf.GEOTHERMAL_SCALE   = cfg.geothermal_scale
    vhf[:]                 = hf.astype(np.float32)

    ds.title    = (f"MOM6 geothermal {cfg.niglobal}x{cfg.njglobal} "
                   f"({cfg.dlon:.4g} deg)")
    ds.NIGLOBAL = cfg.niglobal
    ds.NJGLOBAL = cfg.njglobal
    ds.history  = "Created by create_geothermal.py"
    ds.close()

    ocean = hf > 0
    print(f"[OK] {outpath}")
    print(f"     Fluxo: [{hf[ocean].min():.4f}, {hf[ocean].max():.3f}] W m^-2")
    print(f"     Media oceanica: {hf[ocean].mean():.4f} W m^-2  (tipico: ~{MEAN_HF})")
    if abs(hf[ocean].mean() - MEAN_HF) / MEAN_HF > 0.5:
        print("     AVISO: Media muito diferente do valor tipico - verifique unidades!")


def main():
    parser = argparse.ArgumentParser(
        description="Cria ./INPUT/geothermal_heating_cm2g.nc para o MOM6 - resolucao generica",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--list-presets", action="store_true")
    parser.add_argument("--raw-file",   help="Dado bruto de fluxo geotermico (NetCDF)")
    parser.add_argument("--synthetic",  action="store_true")
    parser.add_argument("--topog-file", required=True,
                        help="./INPUT/topog.nc para mascara oceanica (obrigatorio)")
    parser.add_argument("--output-dir", default="INPUT")

    add_grid_args(parser)
    args = parser.parse_args()

    if args.list_presets:
        print_presets(); sys.exit(0)
    if not args.raw_file and not args.synthetic:
        parser.error("Informe --raw-file ou --synthetic")

    cfg = GridConfig(args)
    print("\nCriando ./INPUT/geothermal_heating_cm2g.nc")
    print(cfg.summary())

    lons_m, lats_m = cfg.model_grid()

    print(f"\nLendo mascara de {args.topog_file}...")
    ocean_mask = read_topog_mask(args.topog_file, cfg)

    if args.synthetic:
        hf = synthetic_hf(cfg, lons_m, lats_m)
    else:
        print(f"Lendo {args.raw_file}...")
        lon_r, lat_r, hf_r = read_raw_hf(args.raw_file)
        print("Interpolando...")
        hf = interpolate_hf(lon_r, lat_r, hf_r, lons_m, lats_m,
                            cfg.njglobal, cfg.niglobal)

    print("Aplicando mascara e limites...")
    hf = apply_mask_and_limits(hf, ocean_mask)

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_geothermal(cfg, hf, lons_m, lats_m,
                     outdir / "geothermal_heating_cm2g.nc")

    print(f"""
Fontes recomendadas:
  Davies (2013):   https://doi.org/10.1002/ggge.20271
  Lucazeau (2019): https://doi.org/10.1029/2019GC008464

GEOTHERMAL_SCALE = {GEOTHERMAL_SCALE} (no MOM_input)
  Arquivo em W m^-2 - SCALE=1e-3 -
  Arquivo em mW m^-2 - SCALE=1e-6
""")


if __name__ == "__main__":
    main()
EOF
python3 create_geothermal.py --preset one1 --synthetic --topog-file INPUT/topog.nc
rm create_geothermal.py
