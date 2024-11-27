""" utilities """

import logging
import os
from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray

fmt = "%Y%m%d%H"

na_values = {
    "ID": -1,
    "ZMIN(m)": -9999.90,
    "ZLAT(N)": -9999.90,
    "ZLON(E)": -9999.90,
    "DY(km)": -9999.90,
    "DX(km)": -9999.90,
    "DIST(km)": -9999.90,
    "MAXDUR(h)": -9999.90,
    "PTY-OVR": 9999.90,
    "FERRY(km)": -99999.90,
    "FERRX(km)": -99999.90,
    "FERR(km)": -99999.90,
    "VLat(N)": -9999.90,
    "VLon(E)": -9999.90,
    "VSo": -9999.90,
    "VRo": -9999.90,
    "VZmin": -9999.90,
}

# Mean radius of the earth (in km)
EARTH_RADIUS = 6371.009


@lru_cache(maxsize=15)
def get_obsds(time, **kwargs):
    """
    get GFS 0-h forecast

    convert longitude range to -180, +180
    sort latitude and longitude
    """
    logging.info(f"get_obsds {time} {kwargs}")

    # Translate var to cfVarName
    if "var" in kwargs:
        var = kwargs.pop("var")
        if var == "z":
            kwargs.update(cfVarName="gh")
        else:
            kwargs.update(cfVarName=var)

    obs_file = (
        "/glade/campaign/collections/rda/data/ds084.1/"
        f"{time.strftime('%Y')}/{time.strftime('%Y%m%d')}/"
        f"gfs.0p25.{time.strftime('%Y%m%d%H')}.f000.grib2"
    )
    obsds = xarray.open_dataset(
        obs_file,
        engine="cfgrib",
        backend_kwargs={"indexpath": f"{os.getenv('TMPDIR')}/cfgrib_index_{hash(obs_file)}"},
        filter_by_keys={"typeOfLevel": "isobaricInhPa", **kwargs},
    )

    # don't need `step`. For f000 it is always zero. `time` is same as `valid_time`.
    #obsds = obsds.drop_vars(["step", "time"])

    obsds = obsds.rename(longitude="lon", latitude="lat", gh="z")

    # Convert longitudes from 0-360 to -180-180
    obsds = obsds.assign_coords(lon=((obsds["lon"] + 180) % 360) - 180)
    # Sort lat and lon to maintain order and let slices work
    obsds = obsds.sortby(["lon", "lat"])
    return obsds


def haversine(point1, point2):
    """
    Calculate the great circle distance between two points
    on the Earth, specified as (lon, lat), where lon and lat
    are in degrees.

    Returns: distance between points in km
    """
    # convert decimal degrees to radians
    point1 = np.radians(point1)
    point2 = np.radians(point2)
    lon1 = point1.iloc[:, 0]
    lat1 = point1.iloc[:, 1]
    lon2 = point2.iloc[:, 0]
    lat2 = point2.iloc[:, 1]

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS * np.arcsin(np.sqrt(a))


def handleTimestamp(timestamp):
    # Case 1: If it's a pandas Timestamp or datetime-like object
    if isinstance(timestamp, pd.Timestamp):
        return timestamp

    # Case 2: If it's a numpy datetime64
    elif isinstance(timestamp, np.datetime64):
        # Convert to pandas Timestamp for strftime support
        return pd.to_datetime(timestamp)

    # Case 3: If it's an xarray DataArray
    elif isinstance(timestamp, xarray.DataArray):
        # Assuming the DataArray contains datetime64 values, convert to pandas Timestamps
        return pd.to_datetime(timestamp.values)

    else:
        raise TypeError("Unsupported type for timestamp argument")


def getfcst(itime, valid_time, workdir: Path, isensemble=False):
    itime = handleTimestamp(itime)
    valid_time = handleTimestamp(valid_time)
    # workdirs is 1-element list for deterministic
    workdirs = [workdir / itime.strftime(fmt)]
    if isensemble:  # all the members
        workdirs = workdir.glob(f"E{itime.strftime(fmt)}.p??.*")

    fhr = (valid_time - itime) / pd.Timedelta(hours=1)

    fcst = []
    for workdir in workdirs:
        ifile = workdir / f"diag_TroughsCutoffs.{itime.strftime(fmt)}.f{fhr:03.0f}.track"
        print(".", end="")
        fcst.append(
            pd.read_csv(
                ifile,
                header=0,
                sep=r"\s+",
                na_values=na_values,
            )
        )
    return pd.concat(fcst)


def getobs(valid_time):
    valid_time = handleTimestamp(valid_time)
    obs_path = valid_time.strftime(
        "/glade/u/home/klupo/work_new/postdoc/kasugaEA21/version9/HGT_500mb/"
        f"gfs.0p25.%Y%m%d%H.f000.track"
    )
    obs = pd.read_csv(
        obs_path,
        header=0,
        sep=r"\s+",
        na_values=na_values,
    )
    return obs


def stack(ds: xarray.Dataset):
    # group same variable from different levels into new variable with vertical dimension
    potential_plevs = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50, 10]
    for statevar in ["omg", "q", "rh", "t", "u", "v", "z"]:
        for tend in [
            "",
            "cnvgwd",
            "deepcnv",
            "lw",
            "mp",
            "nonphys",
            "ogwd",
            "pbl",
            "rdamp",
            "shalcnv",
            "sw",
        ]:
            potentialvars = [f"d{statevar}3dt{plev}_{tend}" for plev in potential_plevs]
            newvar = f"d{statevar}3dt_{tend}"
            if tend == "":
                potentialvars = [f"{statevar}{plev}" for plev in potential_plevs]
                newvar = statevar
            logging.info(newvar)
            this_vars = [var for var in potentialvars if var in ds]
            if not this_vars:
                continue
            stack = []
            for var in this_vars:
                logging.info(f"{var} {statevar} {tend}")
                plev = int(
                    var.removeprefix(f"d{statevar}3dt")
                    .removeprefix(statevar)
                    .removesuffix(f"_{tend}")
                )
                assert plev, f"unexpected plev {plev}"
                tmpvar = ds[var].expand_dims(pfull=[plev])
                tmpvar.attrs["description"] = tmpvar.attrs["description"].lstrip(f"{plev}-mb ")
                stack.append(tmpvar)
                ds = ds.drop_vars(var)
            ds[newvar] = xarray.concat(stack, dim="pfull")
    return ds


def tissot(ax: plt.axes, df: pd.DataFrame, **kwargs) -> None:
    """
    Draw circle with Ro(km) radius around (LON(E), LAT(N)).

    Parameters
    ----------
    ax : plt.axes
        axes to draw on
    df : pd.DataFrame
        DataFrame with columns ["Ro(km)", "LON(E)", "LAT(N)"]

    Returns
    -------
    None
    """

    if df.empty:
        return
    for i, row in df.iterrows():
        if "color" in row:
            kwargs.update({"color": row["color"]})
        t = ax.tissot(
            rad_km=row["Ro(km)"],
            lons=row["LON(E)"],
            lats=row["LAT(N)"],
            **kwargs,
        )
    return t
