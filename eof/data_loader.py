"""Carga las variables atmosfericas (ERA5, 1940-2024) usadas por la demo."""

from pathlib import Path

import xarray as xr

DATA_DIR = Path(__file__).resolve().parent / "data"

VARIABLES = {
    "Temperatura (2 m)": dict(
        file=DATA_DIR / "t2m1940_2024.nc",
        var="t2m",
        unit="°C",
        offset=-273.15,
        scale=1.0,
        cmap="RdBu_r",
    ),
    "Precipitación": dict(
        file=DATA_DIR / "pp1940_2024.nc",
        var="tp",
        unit="mm/día",
        offset=0.0,
        scale=1000.0,
        cmap="BrBG",
    ),
}


def load_variable(key):
    """Retorna (valores (time,lat,lon), lat, lon, time, meses, unidad, cmap)."""
    cfg = VARIABLES[key]
    ds = xr.open_dataset(cfg["file"])
    values = ds[cfg["var"]].values * cfg["scale"] + cfg["offset"]
    lat = ds["latitude"].values
    lon = ds["longitude"].values
    time = ds["valid_time"].values
    months = ds["valid_time"].dt.month.values
    return values, lat, lon, time, months, cfg["unit"], cfg["cmap"]
