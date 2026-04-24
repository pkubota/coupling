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
rm INPUT/mom6_cmems_mod_glo_phy_20260101.nc
################################################################################
# PASSO 1 - Visao Geral e Dependencias
################################################################################
cat << 'EOF'>create_ic_cmems.py
"""
create_ic_cmems.py
------------------
Cria o arquivo de condicoes iniciais do MOM6/SIS2 a partir dos
arquivos diarios do CMEMS/GLORYS12 (sistema NEMO grade ORCA025).

Interpolacao para a malha do modelo:
  - Horizontal: usa ocean_hgrid.nc (supergrid FMS) para extrair
    lon/lat dos pontos T (centros de celula) e interpola com
    scipy RegularGridInterpolator.
  - Vertical: usa coordenada sigma2/dz do arquivo hycom (ex:
    hycom1_75_800m.nc) para calcular as profundidades-alvo das
    camadas HYCOM1.  Se o arquivo nao for fornecido, usa grade
    z uniforme.
  - Fill/NaN: apos a interpolacao horizontal, pontos oceanicos
    que ficaram com NaN ou fill_value (terra no CMEMS mas oceano
    no modelo) sao preenchidos com o vizinho oceanico valido mais
    proximo usando cKDTree em coordenadas cartesianas 3D.

Estrutura esperada dos arquivos de entrada (CMEMS anfc 0.083deg):
  cmems_mod_glo_phy_anfc_*_multi-vars_*_DATE.nc   - ssh, sic, sit
  cmems_mod_glo_phy-thetao_*_thetao_*_DATE.nc     - temperatura
  cmems_mod_glo_phy-so_*_so_*_DATE.nc             - salinidade
  cmems_mod_glo_phy-cur_*_uo-vo_*_DATE.nc         - velocidades u, v

Estrutura do arquivo de saida (padrao MOM6):
  Dimensoes: time(1), depth(NK), ny(NY_MOD), nx(NX_MOD),
             lonq(NX_MOD+1), lath(NY_MOD), latq(NY_MOD+1)
  Variaveis: ptemp, salt, u, v, ssh, sic, sit

Uso:
    # Interpolando para a malha do modelo (recomendado):
    python create_ic_cmems.py \\
        --cmems-dir /p/projetos/ioper/data/external/copernicus/2026/01/01/00 \\
        --date 2026-01-01 \\
        --hgrid INPUT/ocean_hgrid.nc \\
        --vcoord INPUT/hycom1_75_800m.nc \\
        --topo   INPUT/topog.nc \\
        --output-dir INPUT

    # Modo sintetico para teste:
    python create_ic_cmems.py --synthetic --date 2026-01-01 \\
        --hgrid INPUT/ocean_hgrid.nc \\
        --vcoord INPUT/hycom1_75_800m.nc \\
        --topo   INPUT/topog.nc \\
        --output-dir INPUT

Dependencias:
    pip install numpy netCDF4 scipy

Referencia de grade:
    ORCA025 (NEMO): 1/12graus ~ 0.083graus global
    Grade C de Arakawa: T/S/SSH em (ny,nx), u em (ny,nx+1), v em (ny+1,nx)
"""

import argparse
import sys
import os
import re
import numpy as np
from pathlib import Path
from datetime import datetime

try:
    import netCDF4 as nc
    from scipy.interpolate import RegularGridInterpolator
    from scipy.spatial import cKDTree
except ImportError:
    sys.exit("Instale: pip install netCDF4 scipy")


FILL_VAL = np.float32(9.96921e+36)


# =============================================================
# Leitura da malha do modelo (ocean_hgrid.nc)
# =============================================================

def read_model_grid(hgrid_path):
    """
    Le ocean_hgrid.nc (supergrid FMS, nx=2*NIGLOBAL+1) e retorna
    as coordenadas dos pontos T (centros de celula).

    Convencao supergrid FMS:
      x(ny, nx)  com ny = 2*NJGLOBAL+1, nx = 2*NIGLOBAL+1
      Pontos T estao nos indices impares: [1,3,5,...] em ambas as
      direcoes ? subgrade [1::2, 1::2]

    Retorna
    -------
    lon_t : (NY_MOD, NX_MOD)  longitudes dos pontos T
    lat_t : (NY_MOD, NX_MOD)  latitudes dos pontos T
    lon_t_1d : (NX_MOD,)  longitude 1-D (linha central)
    lat_t_1d : (NY_MOD,)  latitude 1-D (coluna central)
    NX_MOD, NY_MOD : int
    """
    ds = nc.Dataset(str(hgrid_path))
    x_super = ds["x"][:]   # (ny_super, nx_super)
    y_super = ds["y"][:]
    ds.close()

    # Pontos T: indices impares do supergrid
    lon_t = x_super[1::2, 1::2]   # (NY_MOD, NX_MOD)
    lat_t = y_super[1::2, 1::2]

    NY_MOD, NX_MOD = lon_t.shape

    # 1-D para uso como eixos de interpolacao
    lon_t_1d = lon_t[NY_MOD // 2, :]   # linha do meio (latitude media)
    lat_t_1d = lat_t[:, NX_MOD // 2]   # coluna do meio

    print(f"  Grade do modelo (pontos T): {NY_MOD} x {NX_MOD}")
    print(f"  lon T: {lon_t.min():.2f} a {lon_t.max():.2f}")
    print(f"  lat T: {lat_t.min():.2f} a {lat_t.max():.2f}")

    return lon_t, lat_t, lon_t_1d, lat_t_1d, NX_MOD, NY_MOD


def read_ocean_mask(topo_path, ny, nx):
    """
    Le topog.nc e retorna mascara booleana (True = oceano).
    Usa depth > 0 como criterio.
    """
    if topo_path is None or not Path(topo_path).exists():
        print("  [AVISO] topog.nc nao fornecido - assumindo tudo oceano")
        return np.ones((ny, nx), dtype=bool)

    ds = nc.Dataset(str(topo_path))
    depth = np.array(ds["depth"][:])
    ds.close()

    if depth.shape != (ny, nx):
        print(f"  [AVISO] topog shape {depth.shape} != ({ny},{nx})"
              " - assumindo tudo oceano")
        return np.ones((ny, nx), dtype=bool)

    mask = depth > 0
    print(f"  Mascara oceano: {mask.sum()} pontos oceanicos de {ny*nx}")
    return mask


def read_model_vcoord(vcoord_path, nk):
    """
    Le o arquivo de coordenada vertical HYCOM (hycom1_75_800m.nc).
    Le coordenada vertical HYCOM. Usa 'dz' para calcular
    profundidades dos centros das camadas. Fallback: grade z uniforme.

    Espera dimensao 'Layer' (NK) ou 'Interface' (NK+1).
    Usa 'dz' (espessura minima) para calcular profundidades-alvo
    dos centros de camada.

    Se o arquivo nao existir ou nao tiver 'dz', usa grade z uniforme.

    Retorna
    -------
    depth_out : (nk,) profundidades dos centros das camadas [m]
    """
    if vcoord_path is None or not Path(vcoord_path).exists():
        print("  [AVISO] vcoord nao fornecido - usando grade z uniforme")
        return build_depth_uniform(nk)

    ds = nc.Dataset(str(vcoord_path))

    if "dz" in ds.variables:
        dz = ds["dz"][:]
        # remove Interface extra se existir
        if len(dz) == nk + 1:
            dz = dz[:nk]
        elif len(dz) != nk:
            print(f"  [AVISO] dz tem {len(dz)} elementos, esperado {nk} "
                  "- usando grade z uniforme")
            ds.close()
            return build_depth_uniform(nk)
        # centros = interfaces acumuladas - dz/2
        zi = np.cumsum(dz)
        depth_out = zi - dz / 2.0
        print(f"  Coordenada vertical (HYCOM dz): "
              f"{depth_out[0]:.2f} a {depth_out[-1]:.2f} m")
    elif "sigma2" in ds.variables:
        # fallback: distribui uniformemente ate max_depth=6500m
        print("  [AVISO] 'dz' ausente, usando grade z uniforme")
        ds.close()
        return build_depth_uniform(nk)
    else:
        print("  [AVISO] vcoord sem 'dz' - usando grade z uniforme")
        ds.close()
        return build_depth_uniform(nk)

    ds.close()
    return depth_out


def build_depth_uniform(nk, max_depth=6500.0):
    """Grade z uniforme com centros das camadas."""
    dz = max_depth / nk
    return np.arange(0.5 * dz, max_depth, dz)


# =============================================================
# Preenchimento de pontos fill/NaN via cKDTree
# =============================================================

def ll2xyz(lon, lat):
    """Lon/lat (graus) -> coordenadas cartesianas 3D na esfera unitaria."""
    lon_r = np.deg2rad(lon)
    lat_r = np.deg2rad(lat)
    return np.column_stack([
        np.cos(lat_r).ravel() * np.cos(lon_r).ravel(),
        np.cos(lat_r).ravel() * np.sin(lon_r).ravel(),
        np.sin(lat_r).ravel()])


def fill_missing_kdtree(data, lon_t, lat_t, ocean_mask,
                         fill_threshold=None, default_val=0.0):
    """
    Preenche pontos oceanicos com NaN ou fill_value usando o vizinho
    oceanico valido mais proximo (cKDTree em coordenadas 3D).

    Parametros
    ----------
    data          : (nz, ny, nx) ou (ny, nx)
    lon_t         : (ny, nx)
    lat_t         : (ny, nx)
    ocean_mask    : (ny, nx) True = oceano no modelo
    fill_threshold: limiar acima do qual considera fill (default FILL_VAL*0.9)
    default_val   : valor usado se nenhum ponto valido existir no nivel

    Retorna
    -------
    data_out : mesmo shape, sem NaN/fill em pontos oceanicos
    """
    if fill_threshold is None:
        fill_threshold = float(FILL_VAL) * 0.9

    is_2d = (data.ndim == 2)
    if is_2d:
        data = data[np.newaxis]

    nz, ny, nx = data.shape
    data_out = data.copy().astype(np.float32)
    n_fixed_total = 0

    for k in range(nz):
        lev = data_out[k]
        ok  = ocean_mask & (lev < fill_threshold) & ~np.isnan(lev)
        bad = ocean_mask & ((lev >= fill_threshold) | np.isnan(lev))

        if bad.sum() == 0:
            continue

        if ok.sum() == 0:
            # sem pontos validos: herda do nivel anterior ou usa default
            src = data_out[k-1] if k > 0 else None
            if src is not None:
                ok_prev = ocean_mask & (src < fill_threshold) & ~np.isnan(src)
                data_out[k][bad] = src[bad] if ok_prev.sum() > 0 else default_val
            else:
                data_out[k][bad] = default_val
            n_fixed_total += int(bad.sum())
            continue

        xyz_ok  = ll2xyz(lon_t[ok],  lat_t[ok])
        xyz_bad = ll2xyz(lon_t[bad], lat_t[bad])
        tree    = cKDTree(xyz_ok)
        _, idx  = tree.query(xyz_bad, k=1, workers=-1)
        data_out[k][bad] = lev[ok][idx]
        n_fixed_total += int(bad.sum())

    # pontos de terra ficam com fill_value
    for k in range(nz):
        data_out[k][~ocean_mask] = float(FILL_VAL)

    if n_fixed_total > 0:
        print(f"    cKDTree: {n_fixed_total} pontos corrigidos em {nz} niveis")

    return data_out[0] if is_2d else data_out


# =============================================================
# Localizacao e leitura dos arquivos CMEMS
# =============================================================

def find_cmems_files(cmems_dir, date_str):
    d = Path(cmems_dir)
    date_nodash = date_str.replace("-", "")

    patterns = {
        "multi":  [f"*multi-vars*{date_str}*", f"*multi*{date_nodash}*"],
        "thetao": [f"*thetao*{date_str}*",     f"*thetao*{date_nodash}*"],
        "so":     [f"*-so_*{date_str}*",       f"*-so_*{date_nodash}*",
                   f"*_so_*{date_str}*"],
        "cur":    [f"*uo-vo*{date_str}*",      f"*cur*{date_str}*",
                   f"*uo-vo*{date_nodash}*"],
    }

    found = {}
    for key, pats in patterns.items():
        for pat in pats:
            matches = list(d.glob(pat))
            if matches:
                found[key] = matches[0]
                break
        if key not in found:
            all_nc = list(d.glob("*.nc"))
            print(f"  [AVISO] Arquivo '{key}' nao encontrado em {d}")
            for f in all_nc:
                print(f"            {f.name}")

    return found


def open_var(filepath, *var_names):
    ds = nc.Dataset(str(filepath))
    for name in var_names:
        if name in ds.variables:
            return ds, ds.variables[name], name
    available = list(ds.variables.keys())
    ds.close()
    raise ValueError(
        f"Nenhuma de {var_names} encontrada em {filepath}.\n"
        f"Variaveis disponiveis: {available}"
    )


class CmemsReader:
    def __init__(self, files):
        self.files = files
        self.lon_t = None
        self.lat_t = None
        self.depth = None

    def _squeeze(self, arr):
        arr = np.ma.filled(arr, np.nan)
        while arr.ndim > 3 and arr.shape[0] == 1:
            arr = arr[0]
        return arr

    def read_coords(self, ds):
        lon_names   = ["longitude", "lon", "nx", "x", "nav_lon"]
        lat_names   = ["latitude",  "lat", "ny", "y", "nav_lat",
                       "lath"]
        depth_names = ["depth", "lev", "olevel", "deptht", "z"]

        lon = next((ds.variables[v][:] for v in lon_names
                    if v in ds.variables), None)
        lat = next((ds.variables[v][:] for v in lat_names
                    if v in ds.variables), None)
        dep = next((ds.variables[v][:] for v in depth_names
                    if v in ds.variables), None)

        if lon is None or lat is None:
            raise ValueError(
                f"Coordenadas lon/lat nao encontradas. "
                f"Vars: {list(ds.variables.keys())}")

        return np.array(lon), np.array(lat), \
               (np.array(dep) if dep is not None else None)

    def read_thetao(self):
        fp = self.files.get("thetao") or self.files.get("multi")
        if fp is None:
            raise FileNotFoundError("Arquivo de temperatura nao encontrado")
        ds, var, name = open_var(fp, "thetao", "ptemp", "votemper", "temp")
        lon, lat, dep = self.read_coords(ds)
        data = self._squeeze(var[:])
        ds.close()
        if self.lon_t is None:
            self.lon_t = lon if lon.ndim == 1 else lon[0, :]
            self.lat_t = lat if lat.ndim == 1 else lat[:, 0]
            self.depth = dep
        print(f"  ptemp ({name}): shape={data.shape}  "
              f"[{np.nanmin(data):.2f}, {np.nanmax(data):.2f}] grausC")
        return data

    def read_so(self):
        fp = self.files.get("so") or self.files.get("multi")
        if fp is None:
            raise FileNotFoundError("Arquivo de salinidade nao encontrado")
        ds, var, name = open_var(fp, "so", "salt", "vosaline", "salinity")
        data = self._squeeze(var[:])
        ds.close()
        print(f"  salt  ({name}): shape={data.shape}  "
              f"[{np.nanmin(data):.2f}, {np.nanmax(data):.2f}] psu")
        return data

    def read_ssh(self):
        fp = self.files.get("multi")
        if fp is None:
            raise FileNotFoundError("Arquivo multi-vars nao encontrado")
        ds, var, name = open_var(fp, "ssh", "zos", "sossheig",
                                 "sea_surface_height")
        data = self._squeeze(var[:])
        ds.close()
        print(f"  ssh   ({name}): shape={data.shape}  "
              f"[{np.nanmin(data):.3f}, {np.nanmax(data):.3f}] m")
        return data

    def read_sic_sit(self):
        fp = self.files.get("multi")
        sic = sit = None
        if fp is not None:
            ds = nc.Dataset(str(fp))
            for name in ["siconc", "sic", "ice_concentration", "CN"]:
                if name in ds.variables:
                    sic = self._squeeze(ds.variables[name][:])
                    break
            for name in ["sithick", "sit", "ice_thickness", "HI"]:
                if name in ds.variables:
                    sit = self._squeeze(ds.variables[name][:])
                    break
            ds.close()
        return sic, sit

    def read_uv(self):
        fp = self.files.get("cur") or self.files.get("multi")
        if fp is None:
            raise FileNotFoundError("Arquivo de correntes nao encontrado")
        ds = nc.Dataset(str(fp))
        u_names = ["uo", "u", "vozocrtx", "uocrtx", "ut"]
        v_names = ["vo", "v", "vomecrty", "vocrtx",  "vt"]
        u_var = next((ds.variables[n] for n in u_names
                      if n in ds.variables), None)
        v_var = next((ds.variables[n] for n in v_names
                      if n in ds.variables), None)
        if u_var is None or v_var is None:
            ds.close()
            raise ValueError(
                f"u/v nao encontrados. Vars: {list(ds.variables.keys())}")
        u_data = self._squeeze(u_var[:])
        v_data = self._squeeze(v_var[:])
        ds.close()
        print(f"  u: shape={u_data.shape}  "
              f"[{np.nanmin(u_data):.3f}, {np.nanmax(u_data):.3f}] m/s")
        print(f"  v: shape={v_data.shape}  "
              f"[{np.nanmin(v_data):.3f}, {np.nanmax(v_data):.3f}] m/s")
        return u_data, v_data


# =============================================================
# Interpolacao horizontal para a malha do modelo
# =============================================================

def normalize_lon(lon, ref_lon_min):
    """
    Normaliza longitudes do CMEMS para o mesmo intervalo do modelo.
    ref_lon_min: longitude minima do modelo (ex: -280)
    """
    lon = np.array(lon, dtype=float)
    # converte para o intervalo [ref_lon_min, ref_lon_min+360)
    lon_out = ((lon - ref_lon_min) % 360.0) + ref_lon_min
    return lon_out


def interp_horiz_to_model(data_cmems, lon_cmems_1d, lat_cmems_1d,
                           lon_model_2d, lat_model_2d):
    """
    Interpola dado CMEMS (grade regular 1-D) para a malha 2-D do modelo.
    Pontos fora do dominio CMEMS ficam NaN (tratados por fill_missing_kdtree).

    Parametros
    ----------
    data_cmems   : (nz, ny_src, nx_src) ou (ny_src, nx_src)
    lon_cmems_1d : (nx_src,)  longitudes CMEMS ordenadas
    lat_cmems_1d : (ny_src,)  latitudes CMEMS ordenadas
    lon_model_2d : (NY_MOD, NX_MOD)  longitudes dos pontos T do modelo
    lat_model_2d : (NY_MOD, NX_MOD)  latitudes dos pontos T do modelo

    Retorna
    -------
    data_out : (nz, NY_MOD, NX_MOD) ou (NY_MOD, NX_MOD)
    """
    NY_MOD, NX_MOD = lon_model_2d.shape
    pts = np.column_stack([lat_model_2d.ravel(),
                           lon_model_2d.ravel()])

    # normaliza lon CMEMS para o mesmo intervalo do modelo
    lon_min_model = lon_model_2d.min()
    lon_src = normalize_lon(lon_cmems_1d, lon_min_model)

    # garante ordenacao crescente apos normalizacao
    order = np.argsort(lon_src)
    lon_src = lon_src[order]

    if data_cmems.ndim == 2:
        data_s = data_cmems[:, order]
        interp = RegularGridInterpolator(
            (lat_cmems_1d, lon_src), data_s,
            method="linear", bounds_error=False,
            fill_value=np.nan)
        return interp(pts).reshape(NY_MOD, NX_MOD)

    nz = data_cmems.shape[0]
    data_out = np.full((nz, NY_MOD, NX_MOD), np.nan, dtype=np.float32)
    for k in range(nz):
        data_s = data_cmems[k][:, order]
        interp = RegularGridInterpolator(
            (lat_cmems_1d, lon_src), data_s,
            method="linear", bounds_error=False,
            fill_value=np.nan)
        data_out[k] = interp(pts).reshape(NY_MOD, NX_MOD)
        if (k + 1) % 10 == 0 or k == nz - 1:
            print(f"    camada {k+1}/{nz} ok", end="\r")
    print()
    return data_out


# =============================================================
# Interpolacao vertical
# =============================================================

def normalize_depth_cmems(depth_raw):
    """
    Detecta e converte unidade de profundidade do CMEMS para metros.
    Heuristica: se o valor maximo for > 10000, assume centimetros.
    """
    depth = np.array(depth_raw, dtype=float)
    if depth.max() > 10000.0:
        print(f"  [INFO] depth max={depth.max():.0f} ? convertendo cm?m")
        depth /= 100.0
    return depth


def interp_vertical(data_in, depth_in, depth_out):
    """
    Interpola verticalmente de depth_in para depth_out.
    data_in  : (nz_in, ny, nx)
    depth_in : (nz_in,) em metros, positivo para baixo
    depth_out: (nz_out,) em metros, positivo para baixo
    Retorna  : (nz_out, ny, nx)
    """
    nz_in, ny, nx = data_in.shape
    nz_out = len(depth_out)
    data_out = np.full((nz_out, ny, nx), np.nan, dtype=np.float32)

    for j in range(ny):
        for i in range(nx):
            col = data_in[:, j, i]
            valid = ~np.isnan(col)
            if valid.sum() < 2:
                continue
            data_out[:, j, i] = np.interp(
                depth_out,
                depth_in[valid], col[valid],
                left=col[valid][0],
                right=col[valid][-1])

    return data_out

# =============================================================
# Construcao da grade C a partir dos dados em grade T
# CORRIGIDO: media mascarada para evitar contaminacao por fill_value
# =============================================================

def t_to_c_grid(data_t, axis):
    """
    Desloca campo da grade T para grade U (axis='x') ou V (axis='y').
    Usa media mascarada: se um vizinho for fill/NaN, usa o outro.
    Se ambos forem fill/NaN, o resultado fica NaN (tratado depois).

    data_t : (nz, ny, nx)  -- valores fill_value ou NaN em terra
    Retorna: (nz, ny, nx+1) para axis='x'
             (nz, ny+1, nx) para axis='y'
    """
    FILL = float(FILL_VAL)
    THRESH = FILL * 0.9

    nz, ny, nx = data_t.shape

    # Trabalha com float64 internamente para evitar overflow na media
    d = data_t.astype(np.float64)

    # Mascara de terra: True onde e fill ou NaN
    land = (d >= THRESH) | np.isnan(d)
    # Substitui fill/NaN por NaN para operar com np.nanmean
    d_masked = np.where(land, np.nan, d)

    if axis == 'x':
        # Pontos U estao entre T[i] e T[i+1]
        out = np.full((nz, ny, nx + 1), np.nan, dtype=np.float64)
        # Bordas: replica o ponto T mais proximo
        out[:, :, 0]    = d_masked[:, :, 0]
        out[:, :, -1]   = d_masked[:, :, -1]
        # Interior: media dos dois vizinhos T (ignorando NaN)
        left  = d_masked[:, :, :-1]
        right = d_masked[:, :, 1:]
        # Onde ambos sao validos: media; onde so um e valido: usa o valido
        both_valid = ~np.isnan(left) & ~np.isnan(right)
        only_left  = ~np.isnan(left) &  np.isnan(right)
        only_right =  np.isnan(left) & ~np.isnan(right)
        mid = np.full_like(left, np.nan)
        mid[both_valid]  = 0.5 * (left[both_valid] + right[both_valid])
        mid[only_left]   = left[only_left]
        mid[only_right]  = right[only_right]
        out[:, :, 1:-1] = mid

    else:  # 'y'
        # Pontos V estao entre T[j] e T[j+1]
        out = np.full((nz, ny + 1, nx), np.nan, dtype=np.float64)
        out[:, 0, :]    = d_masked[:, 0, :]
        out[:, -1, :]   = d_masked[:, -1, :]
        bot = d_masked[:, :-1, :]
        top = d_masked[:, 1:, :]
        both_valid = ~np.isnan(bot) & ~np.isnan(top)
        only_bot   = ~np.isnan(bot) &  np.isnan(top)
        only_top   =  np.isnan(bot) & ~np.isnan(top)
        mid = np.full_like(bot, np.nan)
        mid[both_valid] = 0.5 * (bot[both_valid] + top[both_valid])
        mid[only_bot]   = bot[only_bot]
        mid[only_top]   = top[only_top]
        out[:, 1:-1, :] = mid

    return out.astype(np.float32)
    
# =============================================================
# Construcao da grade C a partir dos dados em grade T
# =============================================================

def t_to_c_grid_old(data_t, axis, nx_out, ny_out):
    """
    Desloca campo da grade T para grade U (axis='x') ou V (axis='y').
    data_t : (nz, ny, nx)
    Retorna: (nz, ny, nx+1) para axis='x'
             (nz, ny+1, nx) para axis='y'
    """
    nz, ny, nx = data_t.shape
    if axis == 'x':
        out = np.zeros((nz, ny, nx + 1), dtype=np.float32)
        out[:, :, 1:-1] = 0.5 * (data_t[:, :, :-1] + data_t[:, :, 1:])
        out[:, :, 0]    = data_t[:, :, 0]
        out[:, :, -1]   = data_t[:, :, -1]
    else:  # 'y'
        out = np.zeros((nz, ny + 1, nx), dtype=np.float32)
        out[:, 1:-1, :] = 0.5 * (data_t[:, :-1, :] + data_t[:, 1:, :])
        out[:, 0, :]    = data_t[:, 0, :]
        out[:, -1, :]   = data_t[:, -1, :]
    return out


# =============================================================
# Dados sinteticos
# =============================================================

def make_synthetic(lon_t_2d, lat_t_2d, depth_out):
    """Gera campos sinteticos realistas para teste."""
    print("  Gerando campos sinteticos...")
    ny, nx = lon_t_2d.shape
    nk = len(depth_out)

    ptemp = np.zeros((nk, ny, nx), dtype=np.float32)
    salt  = np.zeros((nk, ny, nx), dtype=np.float32)
    for k, z in enumerate(depth_out):
        sst = 28.0 * np.exp(-(lat_t_2d / 30.0)**2)
        ptemp[k] = (sst * np.exp(-z / 500.0) + 2.0 * np.exp(-z / 2000.0)
                    ).astype(np.float32)
        sss = 35.0 + 1.5 * np.cos(np.deg2rad(lat_t_2d))
        salt[k] = (sss + 0.3 * np.exp(-z / 1000.0)).astype(np.float32)

    ssh = (0.5 * np.exp(-(lat_t_2d / 20.0)**2)
           * np.sin(np.deg2rad(lon_t_2d / 2))).astype(np.float32)

    u_t = np.zeros((nk, ny, nx), dtype=np.float32)
    v_t = np.zeros((nk, ny, nx), dtype=np.float32)
    for k, z in enumerate(depth_out):
        ek = np.exp(-z / 300.0)
        u_t[k] = (0.3 * ek * np.sin(np.deg2rad(lat_t_2d * 2))).astype(
            np.float32)
        v_t[k] = (0.1 * ek * np.cos(np.deg2rad(lon_t_2d))).astype(
            np.float32)

    u = t_to_c_grid(u_t, 'x', nx + 1, ny)
    v = t_to_c_grid(v_t, 'y', nx, ny + 1)

    sic = np.where(lat_t_2d > 70, 0.8, 0.0).astype(np.float32)
    sit = np.where(lat_t_2d > 70, 1.5, 0.0).astype(np.float32)

    return ptemp, salt, u, v, ssh, sic, sit


# =============================================================
# Escrita do arquivo de saida
# =============================================================

def write_ic(outpath, date_str, lon_t_1d, lat_t_1d, depth,
             ptemp, salt, u, v, ssh, sic, sit):
    """
    Escreve o arquivo de CI no formato padrao MOM6/SIS2.
    lon_t_1d, lat_t_1d : coordenadas 1-D dos pontos T do modelo
    """
    date_str_local = date_str
    nk, ny, nx = ptemp.shape
    nx_u = u.shape[2]   # nx+1
    ny_v = v.shape[1]   # ny+1

    dlon = (lon_t_1d[-1] - lon_t_1d[0]) / (nx - 1)
    dlat = (lat_t_1d[-1] - lat_t_1d[0]) / (ny - 1)
    lonq = np.append(lon_t_1d - dlon / 2, lon_t_1d[-1] + dlon / 2)
    latq = np.append(lat_t_1d - dlat / 2, lat_t_1d[-1] + dlat / 2)
    print(date_str_local)
    date_str_local = "2016-01-01 01:30:00"
    print(date_str_local)

    # Data -> horas desde 1950-01-01
    # Aceita: "YYYY-MM-DD", "YYYY-MM-DD HH:MM", "YYYY-MM-DDTHH:MM", "YYYY-MM-DD HH:MM:SS"
    t0   = datetime(1950, 1, 1)
    for fmt in ("%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S",
                "%Y-%m-%dT%H:%M",    "%Y-%m-%d %H:%M",
                "%Y-%m-%dT%H",       "%Y-%m-%d"):
        try:
            tval = datetime.strptime(date_str_local, "%Y-%m-%d %H:%M:%S")
            break
        except ValueError:
            continue
    else:
        raise ValueError(
            f"Formato de data nao reconhecido: '{date_str}'\n"
            f"Use: YYYY-MM-DD | YYYY-MM-DD HH:MM | YYYY-MM-DDTHH:MM"
        )

    hours_since = (tval - t0).total_seconds() / 3600.0

    ds = nc.Dataset(str(outpath), "w", format="NETCDF4_CLASSIC")

    ds.createDimension("time",  1)
    ds.createDimension("nx",    nx)
    ds.createDimension("ny",    ny)
    ds.createDimension("lonq",  nx_u)
    ds.createDimension("lath",  ny)
    ds.createDimension("latq",  ny_v)
    ds.createDimension("depth", nk)

    def coord(name, dim, data, **attrs):
        v = ds.createVariable(name, "f8", (dim,))
        for k, val in attrs.items():
            setattr(v, k, val)
        v[:] = data

    vt = ds.createVariable("time", "f4", ("time",))
    vt.standard_name = "time"
    vt.units         = "hours since 1950-01-01"
    vt.calendar      = "gregorian"
    vt.axis          = "T"
    vt[:]            = [hours_since]

    coord("nx",   "nx",   lon_t_1d,
          standard_name="longitude",
          long_name="Nominal Longitude of T-cell center",
          units="degrees_east", axis="X")
    coord("ny",   "ny",   lat_t_1d,
          standard_name="latitude",
          long_name="Nominal Latitude of T-cell center",
          units="degrees_north", axis="Y")
    coord("lonq", "lonq", lonq,
          standard_name="longitude",
          long_name="Longitude of U-cell center",
          units="degrees_east", axis="X")
    coord("lath", "lath", lat_t_1d,
          standard_name="latitude",
          long_name="Latitude of V-cell center (same as T)",
          units="degrees_north", axis="Y")
    coord("latq", "latq", latq,
          standard_name="latitude",
          long_name="Latitude of corner points",
          units="degrees_north", axis="Y")

    vd = ds.createVariable("depth", "f8", ("depth",))
    vd.units    = "m"
    vd.positive = "down"
    vd.axis     = "Z"
    vd[:]       = depth

    def ocean_var(name, dims, data, **attrs):
        v = ds.createVariable(name, "f4", dims, fill_value=FILL_VAL)
        v.missing_value = FILL_VAL
        for k, val in attrs.items():
            setattr(v, k, val)
        d = np.where(np.isnan(data), float(FILL_VAL), data).astype(np.float32)
        if d.ndim == 2:
            v[:] = d[np.newaxis]
        else:
            v[:] = d[np.newaxis]

    ocean_var("ptemp", ("time","depth","ny","nx"), ptemp,
              long_name="Potential temperature", units="degrees_C")
    ocean_var("salt",  ("time","depth","ny","nx"), salt,
              long_name="Salinity", units="psu")
    ocean_var("u",     ("time","depth","lath","lonq"), u,
              long_name="Zonal velocity", units="m s-1")
    ocean_var("v",     ("time","depth","latq","nx"), v,
              long_name="Meridional velocity", units="m s-1")
    ocean_var("ssh",   ("time","ny","nx"),
              ssh[np.newaxis] if ssh.ndim == 2 else ssh,
              standard_name="sea_surface_height_above_geoid",
              long_name="Sea surface height", units="m")
    ocean_var("sic",   ("time","ny","nx"),
              sic[np.newaxis] if sic.ndim == 2 else sic,
              standard_name="sea_ice_area_fraction",
              long_name="Ice concentration", units="1")
    ocean_var("sit",   ("time","ny","nx"),
              sit[np.newaxis] if sit.ndim == 2 else sit,
              standard_name="sea_ice_thickness",
              long_name="Sea ice thickness", units="m")

    ds.title   = f"MOM6/SIS2 initial conditions from CMEMS - {date_str}"
    ds.history = "Created by create_ic_cmems.py"
    ds.source  = "CMEMS GLORYS12 / Global Ocean Physics Analysis"
    ds.date    = date_str
    ds.NK = nk; ds.NX = nx; ds.NY = ny
    ds.close()

    size_mb = outpath.stat().st_size / 1e6
    print(f"[OK] {outpath}  ({size_mb:.1f} MB)")
    print(f"     shape: time=1, depth={nk}, ny={ny}, nx={nx}")
    print(f"     u: (time,depth,lath={ny},lonq={nx_u})")
    print(f"     v: (time,depth,latq={ny_v},nx={nx})")


# =============================================================
# Validacao
# =============================================================

def validate_output(outpath):
    print(f"\nValidando {outpath.name}...")
    ds = nc.Dataset(str(outpath))
    required_vars = ["ptemp","salt","u","v","ssh","sic","sit"]
    required_dims = ["time","nx","ny","depth","lonq","lath","latq"]
    ok = True
    for v in required_vars:
        if v in ds.variables:
            d = ds.variables[v][:]
            nnan = int(np.sum(np.isnan(d)))
            print(f"  [OK] {v}: {ds.variables[v].shape}  NaN={nnan}")
        else:
            print(f"  [ERRO] '{v}' ausente"); ok = False
    for d in required_dims:
        if d in ds.dimensions:
            print(f"  [OK] dim '{d}' = {len(ds.dimensions[d])}")
        else:
            print(f"  [ERRO] dim '{d}' ausente"); ok = False
    nx  = len(ds.dimensions["nx"])
    nxu = len(ds.dimensions["lonq"])
    ny  = len(ds.dimensions["ny"])
    nyv = len(ds.dimensions["latq"])
    print(f"  Grade C: lonq={nxu} ({'OK' if nxu==nx+1 else 'ERRO'}), "
          f"latq={nyv} ({'OK' if nyv==ny+1 else 'ERRO'})")
    ds.close()
    return ok


# =============================================================
# Main
# =============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Cria CI do MOM6/SIS2 interpolada para a malha do modelo",
        formatter_class=argparse.RawTextHelpFormatter
    )

    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--cmems-dir",
                     help="Diretorio com os arquivos CMEMS da data")
    src.add_argument("--synthetic", action="store_true",
                     help="Gerar campo sintetico para teste")

    parser.add_argument("--file-multi",  help="Arquivo multi-vars (ssh,sic,sit)")
    parser.add_argument("--file-thetao", help="Arquivo de temperatura")
    parser.add_argument("--file-so",     help="Arquivo de salinidade")
    parser.add_argument("--file-cur",    help="Arquivo de correntes (u,v)")

    parser.add_argument("--date", default="2026-01-01",
                        help="Data da CI (YYYY-MM-DD)")
    parser.add_argument("--hgrid", default="INPUT/ocean_hgrid.nc",
                        help="ocean_hgrid.nc (supergrid FMS) "
                             "[default: INPUT/ocean_hgrid.nc]")
    parser.add_argument("--vcoord", default="INPUT/hycom1_75_800m.nc",
                        help="Arquivo de coordenada vertical HYCOM "
                             "[default: INPUT/hycom1_75_800m.nc]")
    parser.add_argument("--topo",   default="INPUT/topog.nc",
                        help="topog.nc para mascara oceano/terra"
                             "[default: INPUT/topog.nc]")
    parser.add_argument("--nk",   type=int, default=75,
                        help="Numero de camadas verticais [default: 75]")
    parser.add_argument("--output-dir",  default="INPUT")
    parser.add_argument("--output-name", default=None)

    args = parser.parse_args()

    date_nodash = args.date.replace("-", "")
    outname = (args.output_name or
               f"mom6_cmems_mod_glo_phy_{date_nodash}.nc")
    outdir  = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / outname

    print(f"\nCriando condicoes iniciais MOM6/SIS2")
    print(f"  Data    : {args.date}")
    print(f"  hgrid   : {args.hgrid}")
    print(f"  vcoord  : {args.vcoord}")
    print(f"  topo    : {args.topo}")
    print(f"  NK      : {args.nk}")
    print(f"  Saida   : {outpath}")

    # ---- Leitura da malha do modelo -------------------------
    print("\nLendo malha do modelo...")
    lon_t_2d, lat_t_2d, lon_t_1d, lat_t_1d, NX_MOD, NY_MOD = \
        read_model_grid(args.hgrid)

    # ---- Leitura da mascara oceano/terra ---------------------
    print("\nLendo mascara oceano/terra...")
    ocean_mask = read_ocean_mask(args.topo, NY_MOD, NX_MOD)

    # ---- Leitura da coordenada vertical ---------------------
    print("\nLendo coordenada vertical...")
    depth_out = read_model_vcoord(args.vcoord, args.nk)

    # ---- Modo sintetico -------------------------------------
    if args.synthetic:
        print("\nGerando campos sinteticos...")
        ptemp, salt, u, v, ssh, sic, sit = \
            make_synthetic(lon_t_2d, lat_t_2d, depth_out)
        write_ic(outpath, args.date, lon_t_1d, lat_t_1d, depth_out,
                 ptemp, salt, u, v, ssh, sic, sit)
        validate_output(outpath)
        return

    # ---- Localizacao dos arquivos CMEMS ---------------------
    print("\nLocalizando arquivos CMEMS...")
    if args.cmems_dir:
        files = find_cmems_files(args.cmems_dir, args.date)
    else:
        files = {}
        if args.file_multi:  files["multi"]  = Path(args.file_multi)
        if args.file_thetao: files["thetao"] = Path(args.file_thetao)
        if args.file_so:     files["so"]     = Path(args.file_so)
        if args.file_cur:    files["cur"]    = Path(args.file_cur)
    for k, v in files.items():
        print(f"  {k:8s}: {v.name}")

    reader = CmemsReader(files)

    # ---- Leitura dos dados CMEMS ----------------------------
    print("\nLendo variaveis CMEMS...")
    ptemp_raw = reader.read_thetao()
    salt_raw  = reader.read_so()
    ssh_raw   = reader.read_ssh()
    sic_raw, sit_raw = reader.read_sic_sit()
    u_raw, v_raw     = reader.read_uv()

    lon_cmems = reader.lon_t
    lat_cmems = reader.lat_t
    depth_in  = normalize_depth_cmems(
        reader.depth if reader.depth is not None
        else np.linspace(0, 6500.0, ptemp_raw.shape[0]))

    # ---- Interpolacao horizontal para a malha do modelo -----
    print("\nInterpolando horizontalmente para a malha do modelo...")
    print("  ptemp...")
    ptemp_h = interp_horiz_to_model(ptemp_raw, lon_cmems, lat_cmems,
                                     lon_t_2d, lat_t_2d)
    print("  salt...")
    salt_h  = interp_horiz_to_model(salt_raw,  lon_cmems, lat_cmems,
                                     lon_t_2d, lat_t_2d)
    print("  ssh...")
    ssh_h   = interp_horiz_to_model(ssh_raw,   lon_cmems, lat_cmems,
                                     lon_t_2d, lat_t_2d)

    sic_raw = sic_raw if sic_raw is not None else \
              np.zeros(ssh_raw.shape, dtype=np.float32)
    sit_raw = sit_raw if sit_raw is not None else \
              np.zeros(ssh_raw.shape, dtype=np.float32)
    print("  sic/sit...")
    sic_h   = interp_horiz_to_model(sic_raw, lon_cmems, lat_cmems,
                                     lon_t_2d, lat_t_2d)
    sit_h   = interp_horiz_to_model(sit_raw, lon_cmems, lat_cmems,
                                     lon_t_2d, lat_t_2d)
    print("  u...")
    # u pode ter nx+1 na grade C ? usa mesma lon_cmems (grade T approx)
    u_lon = lon_cmems
    if u_raw.shape[-1] != len(lon_cmems):
        # grade C deslocada: usa lon_cmems sem o ultimo ponto
        u_lon = lon_cmems[:u_raw.shape[-1]]
    u_h = interp_horiz_to_model(u_raw, u_lon, lat_cmems,
                                  lon_t_2d, lat_t_2d)
    print("  v...")
    v_lat = lat_cmems
    if v_raw.shape[-2] != len(lat_cmems):
        v_lat = lat_cmems[:v_raw.shape[-2]]
    v_h = interp_horiz_to_model(v_raw, lon_cmems, v_lat,
                                  lon_t_2d, lat_t_2d)

    # --- Preenchimento fill/NaN via cKDTree (1a passagem) ---
    print("\nPreenchendo fill/NaN com vizinho mais proximo (cKDTree)...")
    print("  ptemp..."); ptemp_h = fill_missing_kdtree(
        ptemp_h, lon_t_2d, lat_t_2d, ocean_mask, default_val=2.0)
    print("  salt...");  salt_h  = fill_missing_kdtree(
        salt_h,  lon_t_2d, lat_t_2d, ocean_mask, default_val=35.0)
    print("  ssh...");   ssh_h   = fill_missing_kdtree(
        ssh_h,   lon_t_2d, lat_t_2d, ocean_mask, default_val=0.0)
    print("  sic...");   sic_h   = fill_missing_kdtree(
        sic_h,   lon_t_2d, lat_t_2d, ocean_mask, default_val=0.0)
    print("  sit...");   sit_h   = fill_missing_kdtree(
        sit_h,   lon_t_2d, lat_t_2d, ocean_mask, default_val=0.0)
    print("  u...");     u_h     = fill_missing_kdtree(
        u_h,     lon_t_2d, lat_t_2d, ocean_mask, default_val=0.0)
    print("  v...");     v_h     = fill_missing_kdtree(
        v_h,     lon_t_2d, lat_t_2d, ocean_mask, default_val=0.0)

    # ---- Interpolacao vertical para a coord do modelo -------
    print("\nInterpolando verticalmente para a coordenada do modelo...")
    print("  ptemp...")
    ptemp = interp_vertical(ptemp_h, depth_in, depth_out)
    print("  salt...")
    salt  = interp_vertical(salt_h,  depth_in, depth_out)
    print("  u...")
    u_3d  = interp_vertical(u_h,     depth_in, depth_out)
    print("  v...")
    v_3d  = interp_vertical(v_h,     depth_in, depth_out)

    # --- 2a passagem fill apos interpolacao vertical ---
    # (niveis profundos em colunas rasas podem ficar NaN)
    print("\nSegunda passagem fill pos-interpolacao vertical...")
    print("  ptemp..."); ptemp   = fill_missing_kdtree(
        ptemp, lon_t_2d, lat_t_2d, ocean_mask, default_val=2.0)
    print("  salt...");  salt    = fill_missing_kdtree(
        salt,  lon_t_2d, lat_t_2d, ocean_mask, default_val=35.0)
    print("  u...");     u_3d     = fill_missing_kdtree(
        u_3d,     lon_t_2d, lat_t_2d, ocean_mask, default_val=0.0)
    print("  v...");     v_3d     = fill_missing_kdtree(
        v_3d,     lon_t_2d, lat_t_2d, ocean_mask, default_val=0.0)
#
#PK    # converte para grade C
#PK    u = t_to_c_grid(u_3d, 'x', NX_MOD + 1, NY_MOD)
#PK    v = t_to_c_grid(v_3d, 'y', NX_MOD, NY_MOD + 1)
#
# ---- Converte para grade C ANTES do fill final de bordas ----------------
    # t_to_c_grid ja usa media mascarada; pontos de borda ficam NaN
    print("\nConvertendo u/v para grade C (media mascarada)...")
    u = t_to_c_grid(u_3d, 'x')   # (nk, ny, nx+1)
    v = t_to_c_grid(v_3d, 'y')   # (nk, ny+1, nx)

    # ---- Mascara oceano para as grades U e V --------------------------------
    # Grade U: oceano se pelo menos um dos dois T vizinhos e oceano
    ocean_mask_u = np.zeros((NY_MOD, NX_MOD + 1), dtype=bool)
    ocean_mask_u[:, 0]    = ocean_mask[:, 0]
    ocean_mask_u[:, -1]   = ocean_mask[:, -1]
    ocean_mask_u[:, 1:-1] = ocean_mask[:, :-1] | ocean_mask[:, 1:]

    # Grade V: oceano se pelo menos um dos dois T vizinhos e oceano
    ocean_mask_v = np.zeros((NY_MOD + 1, NX_MOD), dtype=bool)
    ocean_mask_v[0, :]    = ocean_mask[0, :]
    ocean_mask_v[-1, :]   = ocean_mask[-1, :]
    ocean_mask_v[1:-1, :] = ocean_mask[:-1, :] | ocean_mask[1:, :]

    # lon/lat aproximadas para as grades U e V (para o cKDTree)
    lon_u = np.zeros((NY_MOD, NX_MOD + 1))
    lon_u[:, 0]    = lon_t_2d[:, 0]
    lon_u[:, -1]   = lon_t_2d[:, -1]
    lon_u[:, 1:-1] = 0.5 * (lon_t_2d[:, :-1] + lon_t_2d[:, 1:])
    lat_u = np.zeros((NY_MOD, NX_MOD + 1))
    lat_u[:, 0]    = lat_t_2d[:, 0]
    lat_u[:, -1]   = lat_t_2d[:, -1]
    lat_u[:, 1:-1] = 0.5 * (lat_t_2d[:, :-1] + lat_t_2d[:, 1:])

    lon_v = np.zeros((NY_MOD + 1, NX_MOD))
    lon_v[0, :]    = lon_t_2d[0, :]
    lon_v[-1, :]   = lon_t_2d[-1, :]
    lon_v[1:-1, :] = 0.5 * (lon_t_2d[:-1, :] + lon_t_2d[1:, :])
    lat_v = np.zeros((NY_MOD + 1, NX_MOD))
    lat_v[0, :]    = lat_t_2d[0, :]
    lat_v[-1, :]   = lat_t_2d[-1, :]
    lat_v[1:-1, :] = 0.5 * (lat_t_2d[:-1, :] + lat_t_2d[1:, :])

    # ---- Fill final nas grades U e V ----------------------------------------
    print("  u (grade C)...")
    u = fill_missing_kdtree(u, lon_u, lat_u, ocean_mask_u, default_val=0.0)
    print("  v (grade C)...")
    v = fill_missing_kdtree(v, lon_v, lat_v, ocean_mask_v, default_val=0.0)
    # ---- Escrita --------------------------------------------
    print(f"\nEscrevendo {outpath}...")
    write_ic(outpath, args.date, lon_t_1d, lat_t_1d, depth_out,
             ptemp, salt, u, v, ssh_h, sic_h, sit_h)

    validate_output(outpath)

    print(f"""
Trecho do MOM_input para esta CI:
----------------------------------
TEMP_SALT_Z_INIT_FILE = "{outname}"
Z_INIT_FILE_PTEMP_VAR = "ptemp"
Z_INIT_FILE_SALT_VAR  = "salt"
TS_FILE               = "{outname}"
SURFACE_HEIGHT_IC_FILE = "{outname}"
SURFACE_HEIGHT_IC_VAR  = "ssh"
----------------------------------
""")


if __name__ == "__main__":
    main()


EOF
python3 create_ic_cmems.py \
        --cmems-dir /p/projetos/ioper/data/external/copernicus/2026/01/01/00 \
        --date 2026-01-01 \
        --hgrid INPUT/ocean_hgrid.nc \
        --vcoord INPUT/hycom1_75_800m.nc \
        --topo   INPUT/topog.nc \
        --output-dir INPUT
