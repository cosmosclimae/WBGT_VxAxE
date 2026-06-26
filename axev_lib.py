#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
axev_lib.py — functions for the A x E x V Snakemake workflow.
Adds: multi-variable hazard (WBGT-level sweep) + windowed population aggregation
(memory-flat, no full 1 km raster in RAM).
"""
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray  # noqa: F401
import rasterio
from rasterio.windows import Window
from rasterio import features
from affine import Affine
import geopandas as gpd


# ----------------------------------------------- hazard -> A (sum of listed vars)
def load_hazard(path, varnames, n_years=30, lon_0_360=True):
    """A = annual-mean days = sum over listed monthly vars, summed over 360 months, /n_years.
    `varnames` is a list: 1 element for a cumulative ge-var, several to sum binned classes."""
    ds = xr.open_dataset(path)
    da = sum(ds[v] for v in varnames)                    # elementwise sum of the listed vars
    dims = list(da.dims)
    lon = next(c for c in dims if c.lower() in ("lon", "longitude", "x"))
    lat = next(c for c in dims if c.lower() in ("lat", "latitude", "y"))
    tdim = next(c for c in dims if c not in (lon, lat))
    A = da.sum(dim=tdim, skipna=True) / n_years
    A = A.rename({lon: "x", lat: "y"})
    if lon_0_360:
        A = A.assign_coords(x=(((A["x"] + 180) % 360) - 180)).sortby("x")
    A = A.sortby("y").rio.write_crs("EPSG:4326")
    A.name = "A_days"
    ds.close()
    return A


def save_grid(A, path):
    A.to_netcdf(path)


def load_grid(path):
    ds = xr.open_dataset(path)
    name = "A_days" if "A_days" in ds.data_vars else \
           [v for v in ds.data_vars if v not in ("spatial_ref", "crs")][0]
    return ds[name].rio.write_crs("EPSG:4326")


# ------------------------------------- population 1 km -> 0.25 grid, by blocks
def load_pop_on_grid(pop_path, grid, block=2048):
    """Sum 1 km population COUNTS onto the 0.25 deg `grid` cells, streaming by blocks.
    RAM stays flat (~ one block), independent of raster size."""
    gx = grid["x"].values; gy = grid["y"].values
    dx = float(gx[1] - gx[0]); dy = float(gy[1] - gy[0])
    x_edges = np.concatenate([gx - dx / 2.0, [gx[-1] + dx / 2.0]])
    y_asc = gy if gy[0] < gy[-1] else gy[::-1]
    ddy = abs(dy)
    y_edges = np.concatenate([y_asc - ddy / 2.0, [y_asc[-1] + ddy / 2.0]])
    nx = len(gx)
    acc = np.zeros((len(y_asc), nx), dtype="float64")

    with rasterio.open(pop_path) as src:
        nodata = src.nodata; T = src.transform; H, W = src.height, src.width
        for row0 in range(0, H, block):
            h = min(block, H - row0)
            for col0 in range(0, W, block):
                w = min(block, W - col0)
                arr = src.read(1, window=Window(col0, row0, w, h)).astype("float64")
                if nodata is not None:
                    arr = np.where(arr == nodata, 0.0, arr)
                arr = np.where(np.isfinite(arr), arr, 0.0)
                if not arr.any():
                    continue
                xs = T.c + T.a * (col0 + np.arange(w) + 0.5)
                ys = T.f + T.e * (row0 + np.arange(h) + 0.5)
                ix = np.digitize(xs, x_edges) - 1
                iy = np.digitize(ys, y_edges) - 1
                vx = (ix >= 0) & (ix < nx)
                vy = (iy >= 0) & (iy < len(y_asc))
                if not vx.any() or not vy.any():
                    continue
                sub = arr[np.ix_(vy, vx)]
                flat = (iy[vy][:, None] * nx + ix[vx][None, :]).ravel()
                np.add.at(acc.ravel(), flat, sub.ravel())

    if gy[0] >= gy[-1]:
        acc = acc[::-1, :]
    P = xr.DataArray(acc, dims=("y", "x"), coords={"y": gy, "x": gx}, name="P")
    return P.rio.write_crs("EPSG:4326")


# ------------------------------------------------------- country id raster
def build_country_id(grid, shp, iso_field):
    g = gpd.read_file(shp)[[iso_field, "geometry"]].rename(columns={iso_field: "iso3"})
    g = g[g["iso3"].notna() & (g["iso3"] != "-99")].reset_index(drop=True)
    iso2id = {iso: i + 1 for i, iso in enumerate(sorted(g["iso3"].unique()))}
    id2iso = {v: k for k, v in iso2id.items()}
    x = grid["x"].values; y = grid["y"].values
    dx = float(x[1] - x[0]); dy = float(y[1] - y[0])
    flip = dy > 0
    if flip:
        y = y[::-1]
    transform = Affine.translation(x[0] - dx / 2, y[0] + abs(dy) / 2) * Affine.scale(dx, -abs(dy))
    shapes = ((geom, iso2id[iso]) for geom, iso in zip(g.geometry, g["iso3"]))
    cid = features.rasterize(shapes, out_shape=(len(y), len(x)), transform=transform,
                             fill=0, dtype="int32", all_touched=False)
    if flip:
        cid = cid[::-1, :]
    return cid, id2iso


# --------------------------------------------------------------- GVI table
def load_gvi(csv, iso="iso_code", year="year", ssp="ssp", val="pgvi"):
    df = pd.read_csv(csv)
    df = df[pd.to_numeric(df[year], errors="coerce").notna()].copy()
    df[year] = df[year].astype(int); df[ssp] = df[ssp].astype(int)
    df[val] = pd.to_numeric(df[val], errors="coerce")
    df = df[df[val].notna() & (df[val] != -9)]
    return df.rename(columns={iso: "iso3", year: "year", ssp: "ssp_gvi", val: "pgvi"})


# ----------------------------------------------- zonal stats for one combo
def zonal(A, P, cid, id2iso, taus):
    b = P * A                                   # person-days (the one multiplication)
    valid = np.isfinite(b) & np.isfinite(P) & (cid > 0)
    c, bb, pp, aa = cid[valid], b[valid], P[valid], A[valid]
    order = np.argsort(c)
    c, bb, pp, aa = c[order], bb[order], pp[order], aa[order]
    uniq, idx = np.unique(c, return_index=True)
    df = pd.DataFrame({"cid": uniq,
                       "charge_c": np.add.reduceat(bb, idx),
                       "pop_c": np.add.reduceat(pp, idx)})
    df["iso3"] = df["cid"].map(id2iso)
    df["a_bar_c"] = np.where(df["pop_c"] > 0, df["charge_c"] / df["pop_c"], np.nan)
    sens = []
    for tau in taus:
        m = aa >= tau
        df[f"N_{tau}"] = np.add.reduceat(np.where(m, pp, 0.0), idx)
        sens.append({"tau": int(tau), "global_N": float(np.where(m, pp, 0.0).sum())})
    return df, sens


# --------------------------------------- concentration index + curve
def concentration(dfc):
    """Pop-weighted, ranked by pgvi ascending.
    CI>0 <=> exposure burden concentrated on HIGH-pgvi (vulnerable) populations."""
    d = dfc.dropna(subset=["pgvi", "a_bar_c", "pop_c"])
    d = d[(d["pop_c"] > 0)].sort_values("pgvi").reset_index(drop=True)
    pop = d["pop_c"].to_numpy(float); a = d["a_bar_c"].to_numpy(float)
    w = pop / pop.sum()
    r = np.cumsum(w) - 0.5 * w
    mu = np.sum(w * a); rbar = np.sum(w * r)
    CI = 2.0 * np.sum(w * (a - mu) * (r - rbar)) / mu if mu > 0 else np.nan
    charge = pop * a
    curve = pd.DataFrame({"cum_pop": np.concatenate([[0.0], np.cumsum(w)]),
                          "cum_charge": np.concatenate([[0.0], np.cumsum(charge) / charge.sum()])})
    return CI, curve


def tercile(s):
    try:
        return pd.qcut(s, 3, labels=["low", "mid", "high"])
    except ValueError:
        return pd.Series(["mid"] * len(s), index=s.index)
