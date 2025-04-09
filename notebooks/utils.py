""" utilities """

import logging
import os
import warnings
from collections.abc import Iterable
from functools import lru_cache
from pathlib import Path
from typing import Any, List, Optional, Tuple

import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray
from pint import Quantity
from scipy.interpolate import griddata

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

TMPDIR = Path(os.getenv("TMPDIR"))


loc2color = {"any": "C0", "CONUS": "C1", "Lupo2023": "C2"}


def get_location(df: pd.DataFrame) -> pd.Series:
    """
    Return descriptive location for each row in df
    given lat/lon, time, and ID
    """

    lon = df["LON(E)"].copy()
    ineg = lon < 0
    lon[ineg] = lon[ineg] + 360
    location = pd.Series("any", index=df.index)
    location[df["LAT(N)"].between(23, 50) & lon.between(360 - 127, 360 - 68)] = "CONUS"
    # TODO: map ID to valid_time, not ITIME. These are only good for f096
    ee = (df["ITIME"] == 2019102206) & (df["ID"] == 32379)
    ee |= (df["ITIME"] == 2019112306) & (df["ID"] == 33027)
    ee |= (df["ITIME"] == 2019121900) & (df["ID"] == 33453)
    ee |= (df["ITIME"] == 2020020806) & (df["ID"] == 34433)
    ee |= (df["ITIME"] == 2020022900) & (df["ID"] == 34916)
    ee |= (df["ITIME"] == 2020040512) & (df["ID"] == 35563)
    ee |= (df["ITIME"] == 2020040812) & (df["ID"] == 35563)
    ee |= (df["ITIME"] == 2020051512) & (df["ID"] == 36483)
    ee |= (df["ITIME"] == 2020090600) & (df["ID"] == 38671)
    ee |= (df["ITIME"] == 2020112406) & (df["ID"] == 40236)
    ee |= (df["ITIME"] == 2021031318) & (df["ID"] == 42385)
    ee |= (df["ITIME"] == 2021042706) & (df["ID"] == 43351)
    ee |= (df["ITIME"] == 2021123006) & (df["ID"] == 48354)
    ee |= (df["ITIME"] == 2022041006) & (df["ID"] == 50591)
    ee |= (df["ITIME"] == 2022050112) & (df["ID"] == 51124)
    ee |= (df["ITIME"] == 2022061018) & (df["ID"] == 52004)
    location[ee] = "Lupo2023"
    return location


@lru_cache(maxsize=15)
def get_obsds(time, **kwargs):
    """
    get GFS 0-h forecast

    convert longitude range to -180, +180
    sort latitude and longitude
    """
    logging.info(f"get_obsds {time} {kwargs}")

    cfVarName = dict(q="r", z="gh")

    # Translate var to cfVarName
    shortName = None
    if "shortName" in kwargs:
        shortName = kwargs.pop("shortName")
        if shortName in cfVarName:
            kwargs.update(cfVarName=cfVarName[shortName])
        else:
            kwargs.update(cfVarName=shortName)

    obs_file = (
        "/glade/campaign/collections/rda/data/d084001/"
        f"{time.strftime('%Y')}/{time.strftime('%Y%m%d')}/"
        f"gfs.0p25.{time.strftime('%Y%m%d%H')}.f000.grib2"
    )
    if shortName is None:
        open_function = xarray.open_dataset
    else:
        open_function = xarray.open_dataarray

    # If a dictionary value is a Quantity, convert to magnitude. i.e. {"level": 500*units.hPa}
    filter_by_keys = {
        key: value.m if isinstance(value, Quantity) else value for key, value in kwargs.items()
    }

    obs = open_function(
        obs_file,
        engine="cfgrib",
        backend_kwargs={"indexpath": f"{os.getenv('TMPDIR')}/cfgrib_index_{hash(obs_file)}"},
        filter_by_keys={"typeOfLevel": "isobaricInhPa", **filter_by_keys},
    )

    obs = obs.rename(longitude="lon", latitude="lat", isobaricInhPa="pfull")
    # Turn cfVarName back to original name `shortName`.
    if isinstance(obs, xarray.Dataset):
        cfVarName_reverse = {v: k for k, v in cfVarName.items()}
        obs = obs.rename(cfVarName_reverse)
    elif isinstance(obs, xarray.DataArray):
        if obs.name in cfVarName.values():
            obs.name = {value: key for key, value in cfVarName.items()}[obs.name]

    # Convert longitudes from 0-360 to -180-180
    obs = obs.assign_coords(lon=((obs["lon"] + 180) % 360) - 180)
    # Sort lat and lon to maintain order and let slices work
    obs = obs.sortby(["lon", "lat"])
    return obs


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


def get_contains_zmin(df: pd.DataFrame) -> pd.Series:
    """Return "contains_zmin" column"""
    lloc = df[["LON(E)", "LAT(N)"]]
    zloc = df[["ZLON(E)", "ZLAT(N)"]]
    contains_zmin = pd.Series(haversine(lloc, zloc) < df["Ro(km)"], index=df.index)
    return contains_zmin


def label_id(ax: plt.Axes, df: pd.DataFrame, **kwargs: Any) -> List[plt.Artist]:
    """
    Plot Zmin locations.
    Draw line from obs to matched fcst (if one exists).

    Args:
        ax (plt.Axes): The matplotlib axes object to plot on.
        df (pd.DataFrame): DataFrame containing location data, including 'LON(E)', 'LAT(N)',
                           'ZLON(E)', 'ZLAT(N)', 'VLon(E)', 'VLat(N)', 'color', and 'ID' columns.
        **kwargs: Additional keyword arguments to pass to the text annotations and plot functions.

    Returns:
        List[plt.Artist]: A list of matplotlib artists (plot objects) created.
    """

    markers: List[str] = [">", "<", "P", "D", "^", "v", "*", "s", "X", "d"]
    ts: List[plt.Artist] = []

    for _, row in df.iterrows():
        x = row["LON(E)"]
        y = row["LAT(N)"]
        color = row["color"]

        if False:
            # Plot Zmin location
            (zmin_plot,) = ax.plot(
                row["ZLON(E)"],
                row["ZLAT(N)"],
                "o",
                transform=ccrs.PlateCarree(),
                color=color,
                markerfacecolor="none",
                label="Zmin",
            )
            ts.append(zmin_plot)

        if "VLon(E)" in row:
            vx = row["VLon(E)"]
            vy = row["VLat(N)"]
        else:
            # analysis
            vx = x
            vy = y

        # Plot line and marker for ID
        logging.info(f"marker and line {vx,x} {vy, y}")
        marker: str = markers[int(row["ID"] % len(markers))]
        (id_plot,) = ax.plot(
            [vx, x],
            [vy, y],
            marker=marker,
            markersize=10,
            color="k",
            markeredgecolor="none",
            alpha=0.5,
            label=int(row["ID"]),
            transform=ccrs.Geodetic(),
        )
        ts.append(id_plot)

    return ts


warnings.filterwarnings("ignore", category=UserWarning, module="cartopy")


def animate(
    valid_time,
    ax,
    itime,
    workdir: Path,
    forecast_length: int,
    ids: Iterable = None,
    isensemble=False,
):
    print(valid_time, end="")
    obs = getobs(valid_time, ids=ids)

    dc_in = getattr(ax, "_annotations", [])
    # Clear only the data plots, not the static features
    for coll in dc_in:
        logging.info(f"remove {coll}")
        coll.remove()

    ax._annotations = []

    alpha = 0.8
    fhr = (valid_time - itime) / pd.Timedelta(hours=1)
    workdirs = [workdir / f"{itime.strftime(fmt)}.F{forecast_length:03d}.C768"]
    if isensemble:  # all the members
        workdirs = workdir.glob(f"E{itime.strftime(fmt)}.p??.*")
        alpha = 0.2

    text_kw = dict(
        fontsize="xx-small",
        horizontalalignment="center",
        verticalalignment="center",
        transform=cartopy.crs.PlateCarree(),
        clip_on=True,
    )

    analysis_blob = tissot(ax, obs, alpha=0.4)
    ax._annotations.extend(analysis_blob)

    for workdir in workdirs:
        print(".", end="")
        df = getfcst(itime, valid_time, workdir.parent, isensemble=isensemble, ids=ids)

        ax.set_title(f"{valid_time} f{fhr:03.0f}")
        fcst_blob = tissot(ax, df, alpha=alpha)
        ax._annotations.extend(fcst_blob)

        visible = df["location"].isin(["CONUS", "Lupo2023"])
        visible &= ~df.ID.isnull()

        if "VLon(E)" in df:
            visible &= ~df["VLon(E)"].isnull()

        text_ids = label_id(
            ax,
            df[visible],
            **text_kw,
        )
        ax._annotations.extend(text_ids)

    # don't show duplicates
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    # legend for Zmin and ID
    ax.legend(by_label.values(), by_label.keys())

    return (ax,)


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


def getfcst(itime, valid_time, workdir: Path, isensemble=False, ids: Iterable = None):
    itime = handleTimestamp(itime)
    valid_time = handleTimestamp(valid_time)
    assert itime <= valid_time, f"valid_time {valid_time} before initialization time {itime}"
    # workdirs is 1-element list for deterministic
    workdirs = [workdir / itime.strftime(fmt)]
    if isensemble:  # all the members
        workdirs = workdir.glob(f"E{itime.strftime(fmt)}.p??.*")

    fhr = (valid_time - itime) / pd.Timedelta(hours=1)

    fcst = []
    fname = (
        f"gfs.0p25.{itime.strftime(fmt)}.f000.track"
        if itime == valid_time
        else f"diag_TroughsCutoffs.{itime.strftime(fmt)}.f{fhr:03.0f}.track"
    )
    for workdir in workdirs:
        ifile = workdir / fname
        # print(".", end="")
        fcst.append(
            pd.read_csv(
                ifile,
                header=0,
                sep=r"\s+",
                na_values=na_values,
            )
        )
    fcst = pd.concat(fcst)
    fcst["location"] = get_location(fcst)
    fcst["color"] = fcst["location"].map(loc2color)

    fcst["contains_zmin"] = get_contains_zmin(fcst)
    fcst["valid_time"] = valid_time

    if ids:
        if not fcst.ID.isin(ids).any():
            logging.warning(f"no ID {ids} in fcst {valid_time}. Return empty DataFrame")
        fcst = fcst[fcst.ID.isin(ids)]

    return fcst


def getobs(valid_time, ids: Iterable = None):
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
    if ids:
        if not obs.ID.isin(ids).any():
            logging.warning(f"no ID {ids} in fcst {valid_time}. Return empty DataFrame")
        obs = obs[obs.ID.isin(ids)]
    obs["color"] = "k"
    obs["contains_zmin"] = get_contains_zmin(obs)

    return obs


def rotate_field_xarray(field: xarray.DataArray, lat_a, lon_a, lat_b, lon_b):
    """
    Rotate a field defined on a set of lat/lon points and interpolate it back to the original grid.

    Parameters:
    - field: xarray.DataArray, field values with 2D coordinates `XLAT` (latitude) and `XLONG` (longitude),
             and optionally a `valid_time` dimension.
    - lat_a, lon_a: Latitude and longitude of the original reference point A.
    - lat_b, lon_b: Latitude and longitude of the new position B.

    Returns:
    - interpolated_field: xarray.DataArray with the same dimensions and coordinates as `field`.
    """
    if (lat_a, lon_a) == (lat_b, lon_b):
        print("Start and end locations are equal. Returning the original array.")
        return field

    # Vectorized lat/lon to Cartesian conversion
    def latlon_to_cartesian(lat, lon):
        lat_rad = np.radians(lat)
        lon_rad = np.radians(lon)
        x = np.cos(lat_rad) * np.cos(lon_rad)
        y = np.cos(lat_rad) * np.sin(lon_rad)
        z = np.sin(lat_rad)
        return np.stack([x, y, z], axis=-1)

    # Vectorized Cartesian to lat/lon conversion
    def cartesian_to_latlon(cartesian):
        x, y, z = cartesian.T
        lat = np.degrees(np.arcsin(z))
        lon = np.degrees(np.arctan2(y, x))
        return np.stack([lat, lon], axis=-1)

    # Compute rotation matrix
    def compute_rotation_matrix(v1, v2):
        v1 = v1 / np.linalg.norm(v1)
        v2 = v2 / np.linalg.norm(v2)
        axis = np.cross(v1, v2)
        axis_norm = np.linalg.norm(axis)
        if axis_norm == 0:  # v1 and v2 are collinear
            return np.eye(3)
        axis = axis / axis_norm
        angle = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
        ux, uy, uz = axis
        c = np.cos(angle)
        s = np.sin(angle)
        C = 1 - c
        return np.array(
            [
                [c + ux**2 * C, ux * uy * C - uz * s, ux * uz * C + uy * s],
                [uy * ux * C + uz * s, c + uy**2 * C, uy * uz * C - ux * s],
                [uz * ux * C - uy * s, uz * uy * C + ux * s, c + uz**2 * C],
            ]
        )

    # Extract latitude and longitude arrays
    lat_grid = field["XLAT"].values
    lon_grid = field["XLONG"].values
    lat_flat = lat_grid.ravel()
    lon_flat = lon_grid.ravel()

    # Calculate rotation
    A = latlon_to_cartesian(lat_a, lon_a)
    B = latlon_to_cartesian(lat_b, lon_b)
    R = compute_rotation_matrix(A, B)

    # Rotate all points in one operation
    S_cartesian = latlon_to_cartesian(lat_flat, lon_flat)
    S_prime_cartesian = S_cartesian @ R.T
    S_prime = cartesian_to_latlon(S_prime_cartesian)
    lat_prime, lon_prime = S_prime[:, 0], S_prime[:, 1]

    # Interpolation function
    def interpolate_field(field_slice):
        return griddata(
            points=(lat_prime, lon_prime),
            values=field_slice.values.ravel(),
            xi=(lat_flat, lon_flat),
            method="linear",
        ).reshape(lat_grid.shape)

    # Handle interpolation for time-dependent or static fields
    if "valid_time" in field.dims:
        interpolated_field = xarray.apply_ufunc(
            interpolate_field,
            field,
            input_core_dims=[["lat", "lon"]],
            vectorize=True,
            output_core_dims=[["lat", "lon"]],
            dask="parallelized",
            output_dtypes=[field.dtype],
        )
    else:
        interpolated_data = interpolate_field(field)
        interpolated_field = xarray.DataArray(
            interpolated_data,
            dims=["lat", "lon"],
            coords={"lat": field["lat"], "lon": field["lon"]},
            attrs=field.attrs,
        )

    # Update metadata
    interpolated_field.attrs.update(dict(lat_a=lat_a, lon_a=lon_a, lat_b=lat_b, lon_b=lon_b))
    return interpolated_field


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


def tissot(ax: plt.axes, df: pd.DataFrame, **kwargs):
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

    ts = []
    if df.empty:
        return ts
    for i, row in df.iterrows():
        if "color" in row:
            kwargs.update({"color": row["color"]})
        t = ax.tissot(
            rad_km=row["Ro(km)"],
            lons=row["LON(E)"],
            lats=row["LAT(N)"],
            **kwargs,
        )
        ts.append(t)
    return ts
