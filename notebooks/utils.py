""" utilities """

from metpy.calc import mixing_ratio_from_relative_humidity
from metpy.units import units
import logging
import os
import warnings
from collections.abc import Iterable, Sized, Iterator
from functools import lru_cache
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union, Dict, TYPE_CHECKING

import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray
from pint import Quantity
from scipy.interpolate import griddata

# Avoid circular import issues for type hinting
if TYPE_CHECKING:
    from typing import Self

# Configure logging to use logging.info instead of print
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

fmt = "%Y%m%d%H"

na_values: Dict[str, Union[int, float]] = {
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
EARTH_RADIUS: float = 6371.009

# Define temporary directory path
# Fallback to /tmp if TMPDIR env var is not set, which is a common default for temporary files.
TMPDIR: Path = Path(os.getenv("TMPDIR", "/tmp"))

loc2color: Dict[str, str] = {"any": "C0", "CONUS": "C1", "Lupo2023": "C2"}

# Specific IDs and timestamps for cutoff events (likely related to research cases or known events)
cutoffs: Dict[int, int] = {
    2019102606: 32379,
    2019112706: 33027,
    2019122300: 33453,
    2020021206: 34433,
    2020022312: 34695,
    2020030400: 34916,
    2020040912: 35563,
    2020041212: 35563,
    2020051912: 36483,
    # 2020090118: 38541,  # west, not east at f048, f072; okay f096
    2020091000: 38671,
    2020102818: 39553,  # ignore 39600 TC
    2020112806: 40236,
    2021031718: 42385,
    # 2021041200: 42948,  # northwest, not east at f072
    # 2021041412: 42904(f072), 43078(f096) # matching is messy
    2021050106: 43351,
    2021101918: 46789,
    2021102906: 46950,
    2022010306: 48354,  # unmatches at f066
    2022021612: 49289,
    2022041406: 50591,
    2022050512: 51124,
    2022061418: 52004,
}


def get_location(df: pd.DataFrame) -> pd.Series:
    """
    Return descriptive geographical location for each row in a DataFrame
    given latitude, longitude, initialization time (ITIME), and ID.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing geographical and meteorological feature data,
        expected to have 'LAT(N)', 'LON(E)', 'ITIME', and 'ID' columns.

    Returns
    -------
    pd.Series
        A Series where each element is a string indicating the location:
        "any", "CONUS", or "Lupo2023".
    """
    lon = df["LON(E)"].copy()
    # Normalize longitude to 0-360 for CONUS check
    ineg = lon < 0
    lon[ineg] = lon[ineg] + 360

    location = pd.Series("any", index=df.index)

    # Identify features within the CONUS bounding box
    location[df["LAT(N)"].between(23, 50) & lon.between(360 - 127, 360 - 68)] = "CONUS"

    # Identify specific "Lupo2023" research cases
    # Note: These ITIME/ID combinations are specific to certain forecast lengths (f096)
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
def get_obsds(time: pd.Timestamp, **kwargs: Any) -> Union[xarray.DataArray, xarray.Dataset]:
    """
    Get GFS 0-h forecast (analysis) data from a GRIB2 file.

    Parameters
    ----------
    time : pd.Timestamp
        The valid time for which to retrieve the GFS data.
    **kwargs : Any
        Additional keyword arguments to filter the GRIB messages.
        `shortName` (str): The short name of the variable to load (e.g., 'q', 't', 'z').
        Other `filter_by_keys` arguments can be passed (e.g., `typeOfLevel`).

    Returns
    -------
    Union[xarray.DataArray, xarray.Dataset]
        An xarray DataArray if `shortName` is specified, otherwise an xarray Dataset.
        Longitudes are converted to -180 to +180 range and sorted.
        If `shortName='q'`, mixing ratio is computed.

    Raises
    ------
    FileNotFoundError
        If the GRIB file for the specified time does not exist.
    ValueError
        If required variables for mixing ratio computation are missing.
    """
    logging.info(f"get_obsds {time} {kwargs}")

    # Mapping of MetPy shortNames to cfgrib's cfVarName for filtering
    cfVarName: Dict[str, str] = dict(q="r", z="gh")  # 'q' uses 'r' (RH) as proxy for filtering

    shortName: Optional[str] = kwargs.pop("shortName", None)

    # Resolve CF variable name for non-'q' fields
    filter_by_keys: Dict[str, Any] = {"typeOfLevel": "isobaricInhPa"}
    for key, value in kwargs.items():
        filter_by_keys[key] = value.m if isinstance(value, Quantity) else value

    if shortName is not None and shortName != "q":
        if shortName in cfVarName:
            filter_by_keys["cfVarName"] = cfVarName[shortName]
        else:
            filter_by_keys["cfVarName"] = shortName

    # Determine the file path
    obs_file: Path = (
        Path("/glade/campaign/collections/rda/data/d084001/")
        / f"{time.strftime('%Y')}"
        / f"{time.strftime('%Y%m%d')}"
        / f"gfs.0p25.{time.strftime('%Y%m%d%H')}.f000.grib2"
    )

    if not obs_file.exists():
        raise FileNotFoundError(f"GRIB file not found: {obs_file}")

    # If shortName is 'q', compute mixing ratio
    if shortName == "q":
        # Load RH (r), Temperature (t), and Pressure
        ds = xarray.open_dataset(
            obs_file,
            engine="cfgrib",
            backend_kwargs={
                "indexpath": f"{TMPDIR}/cfgrib_index_{hash(obs_file)}"},
            filter_by_keys=filter_by_keys, # filter_by_keys is already updated
        ).rename(longitude="lon", latitude="lat", isobaricInhPa="pfull")

        # Ensure required variables exist
        required_vars = {"r", "t", "pfull"}
        if not required_vars.issubset(ds.variables):
            raise ValueError(
                f"Required variables missing from dataset for mixing ratio: {required_vars - set(ds.variables)}"
            )

        rh = ds["r"].metpy.quantify()
        T = ds["t"].metpy.quantify()
        p = ds["pfull"].metpy.quantify()

        # Align pressure with dimensions
        p_broadcast = p.expand_dims(
            {"lat": rh.lat, "lon": rh.lon}).transpose(*rh.dims)

        # Compute mixing ratio
        q = mixing_ratio_from_relative_humidity(p_broadcast, T, rh)
        q.name = "q"
        q = q.metpy.dequantify() # Remove units for consistency with other xarray operations

        # Reattach coords from the original DataArray
        q = q.assign_coords(lat=ds.lat, lon=ds.lon, pfull=ds.pfull)
        q = q.sortby(["lon", "lat"])
        # Normalize longitude to -180 to 180
        q = q.assign_coords(lon=((q["lon"] + 180) % 360) - 180)

        return q

    # If another variable is requested, or full dataset
    obs = xarray.open_dataset( # Always open as dataset first to ensure correct variable names
        obs_file,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": f"{TMPDIR}/cfgrib_index_{hash(obs_file)}"},
        filter_by_keys=filter_by_keys,
    )

    # Rename dimensions
    obs = obs.rename(longitude="lon", latitude="lat", isobaricInhPa="pfull")

    # Rename cfVarName back to shortName for consistency if it was mapped
    if isinstance(obs, xarray.Dataset):
        cfVarName_reverse = {v: k for k, v in cfVarName.items()}
        for cf_name, short_name in cfVarName_reverse.items():
            if cf_name in obs:
                obs = obs.rename({cf_name: short_name})
        if shortName is not None: # If a specific shortName was requested, return that DataArray
            if shortName in obs:
                obs = obs[shortName]
            else:
                logging.warning(f"Requested variable '{shortName}' not found in dataset.")
                # For now, return original dataset if specified shortName is missing.
    elif isinstance(obs, xarray.DataArray):
        if obs.name in cfVarName.values():
            obs.name = {v: k for k, v in cfVarName.items()}[obs.name]

    # Normalize longitude to -180 to 180 and sort
    obs = obs.assign_coords(lon=((obs["lon"] + 180) % 360) - 180)
    obs = obs.sortby(["lon", "lat"])
    return obs


def haversine(point1: pd.DataFrame, point2: pd.DataFrame) -> pd.Series:
    """
    Calculate the great circle distance between two sets of points
    on the Earth, specified as (lon, lat), where lon and lat
    are in degrees.

    Parameters
    ----------
    point1 : pd.DataFrame
        DataFrame with columns for longitude and latitude for the first set of points.
        Expected columns are the first (longitude) and second (latitude).
    point2 : pd.DataFrame
        DataFrame with columns for longitude and latitude for the second set of points.
        Expected columns are the first (longitude) and second (latitude).

    Returns
    -------
    pd.Series
        Distance between corresponding points in km.
    """
    # convert decimal degrees to radians
    point1_rad = np.radians(point1)
    point2_rad = np.radians(point2)

    lon1 = point1_rad.iloc[:, 0]
    lat1 = point1_rad.iloc[:, 1]
    lon2 = point2_rad.iloc[:, 0]
    lat2 = point2_rad.iloc[:, 1]

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * \
        np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS * np.arcsin(np.sqrt(a))


def get_contains_zmin(df: pd.DataFrame) -> pd.Series:
    """
    Return "contains_zmin" column indicating if a feature's radius
    of influence ('Ro(km)') contains its Zmin (minimum geopotential height) location.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing feature data, expected to have
        'LON(E)', 'LAT(N)', 'ZLON(E)', 'ZLAT(N)', and 'Ro(km)' columns.

    Returns
    -------
    pd.Series
        A boolean Series, True if the feature's radius contains Zmin location.
    """
    lloc = df[["LON(E)", "LAT(N)"]]
    zloc = df[["ZLON(E)", "ZLAT(N)"]]
    contains_zmin = pd.Series(haversine(lloc, zloc) < df["Ro(km)"], index=df.index)
    return contains_zmin


class Feature:
    """
    Represents a single geographical or meteorological feature,
    holding its attributes in a Pandas Series.
    Designed to be initialized with either a pd.Series or a single-row pd.DataFrame.
    """
    _data: pd.Series

    def __init__(self, data: Union[pd.Series, pd.DataFrame]):
        """
        Initializes a Feature object.

        Parameters
        ----------
        data : pd.Series or pd.DataFrame
            The attributes of the feature.
            If a DataFrame, it must contain exactly one row.

        Raises
        ------
        ValueError
            If a DataFrame with more than one row is provided.
        TypeError
            If data is not a pandas Series or DataFrame.
        """
        if isinstance(data, pd.Series):
            self._data = data
        elif isinstance(data, pd.DataFrame):
            if len(data) == 1:
                self._data = data.iloc[0]  # Take the first (and only) row as a Series
            else:
                raise ValueError("DataFrame provided for Feature initialization must contain exactly one row.")
        else:
            raise TypeError("Feature must be initialized with a pandas Series or a single-row DataFrame.")

    def __getattr__(self, name: str) -> Any:
        """
        Allows direct attribute access to the underlying Pandas Series data.
        e.g., feature.ID, feature.Type.

        Parameters
        ----------
        name : str
            The name of the attribute to retrieve.

        Returns
        -------
        Any
            The value of the attribute from the underlying Pandas Series.

        Raises
        ------
        AttributeError
            If the attribute is not found in the Feature object or its underlying data.
        """
        # This is called only if the attribute is not found in the usual places
        if name in self._data.index:
            return self._data[name]
        # Fallback to default behavior for other missing attributes
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def __getitem__(self, key: Union[str, Any]) -> Any:
        """
        Allows subscript notation (e.g., feature["FERRX(km)"]) to access
        the underlying Pandas Series data.

        Parameters
        ----------
        key : str or Any
            The key (column name or index label) to retrieve from the underlying Pandas Series.

        Returns
        -------
        Any
            The value associated with the given key from the underlying Pandas Series.

        Raises
        ------
        KeyError
            If the key is not found in the underlying data.
        """
        try:
            return self._data[key]
        except KeyError:
            raise KeyError(f"'{key}' not found in the feature data.")

    def __repr__(self) -> str:
        """
        Provides a developer-friendly string representation of the Feature object.
        """
        return f"Feature(\n{self._data.__repr__()}\n)"

    def get_latlon(self) -> Tuple[float, float]:
        """
        Retrieves the latitude and longitude of the feature.

        Returns
        -------
        tuple: A tuple containing (latitude, longitude).

        Raises
        ------
        KeyError: If 'LAT(N)' or 'LON(E)' columns are not found in the feature data.
        """
        if "LAT(N)" in self._data.index and "LON(E)" in self._data.index:
            return float(self._data["LAT(N)"]), float(self._data["LON(E)"])
        else:
            raise KeyError("Attributes 'LAT(N)' and/or 'LON(E)' not found in feature data.")

    def to_series(self) -> pd.Series:
        """
        Returns the underlying Pandas Series data.
        """
        return self._data

    def to_dataframe(self) -> pd.DataFrame:
        """
        Returns the underlying data as a single-row DataFrame.
        """
        return pd.DataFrame([self._data])

    def subset_lonlat(self, extent: Tuple[float, float, float, float]) -> Optional['Feature']:
        """
        Checks if the feature's location falls within the specified longitude and latitude extent.
        Longitudes are normalized to -180 to 180 for the check.

        Parameters
        ----------
        extent : Tuple[float, float, float, float]
            A tuple representing (lon0, lon1, lat0, lat1) for the bounding box.

        Returns
        -------
        Optional[Feature]
            The Feature object itself if its location is within the extent, otherwise None.
            Returns None if 'LON(E)' or 'LAT(N)' attributes are missing.
        """
        try:
            lon0, lon1, lat0, lat1 = extent
            lat, lon = self.get_latlon()

            # Normalize longitude of the feature to -180 to 180 range for consistent comparison
            # This is important if `LON(E)` from file is 0-360
            if lon >= 180:
                lon -= 360

            if lat0 <= lat <= lat1 and lon0 <= lon <= lon1:
                return self
            else:
                return None
        except KeyError:
            logging.warning(f"Feature (ID: {getattr(self, 'ID', 'N/A')}) missing 'LON(E)' or 'LAT(N)' attributes for subsetting.")
            return None


class FeatureCollection(Sized, Iterable):
    """
    A collection of Feature objects, providing methods for convenient
    manipulation and conversion to DataFrame.
    """
    _features: List[Feature]

    def __init__(self, features: Optional[List[Feature]] = None):
        """
        Initializes a FeatureCollection.

        Parameters
        ----------
        features : Optional[List[Feature]], optional
            A list of Feature objects to initialize the collection, by default None.
        """
        self._features = features if features is not None else []

    def __len__(self) -> int:
        """Returns the number of features in the collection."""
        return len(self._features)

    def __getitem__(self, index: Union[int, slice]) -> Union[Feature, 'FeatureCollection']:
        """
        Allows indexing and slicing of the feature collection.

        Parameters
        ----------
        index : Union[int, slice]
            The index or slice to apply.

        Returns
        -------
        Union[Feature, FeatureCollection]
            A single Feature if an integer index is used, or a new FeatureCollection
            if a slice is used.
        """
        if isinstance(index, int):
            return self._features[index]
        elif isinstance(index, slice):
            return FeatureCollection(self._features[index])
        else:
            raise TypeError("FeatureCollection indices must be integers or slices.")

    def __iter__(self) -> Iterator[Feature]:
        """Allows iteration over the features in the collection."""
        return iter(self._features)

    def __repr__(self) -> str:
        """Provides a string representation of the FeatureCollection."""
        return f"FeatureCollection with {len(self)} features."

    def to_dataframe(self) -> pd.DataFrame:
        """
        Converts the FeatureCollection into a single Pandas DataFrame.

        Returns
        -------
        pd.DataFrame
            A DataFrame where each row represents a Feature.
        """
        if not self._features:
            # Return an empty DataFrame with expected columns, derived from na_values keys
            # and additional columns usually added to track data.
            default_cols = list(na_values.keys()) + ['color', 'location', 'contains_zmin', 'valid_time']
            return pd.DataFrame(columns=default_cols)
        return pd.DataFrame([feature.to_series() for feature in self._features])

    def filter_by_id(self, ids: Iterable[int]) -> 'FeatureCollection':
        """
        Filters the collection to include only features whose ID is in the given list.

        Parameters
        ----------
        ids : Iterable[int]
            An iterable of integer IDs to filter by.

        Returns
        -------
        FeatureCollection
            A new FeatureCollection containing only the filtered features.
        """
        if not self._features:
            return FeatureCollection([])

        # Get IDs from the current features (convert to set for efficient lookup)
        current_ids = {feature.ID for feature in self._features if hasattr(feature, 'ID')}

        # Check if any requested ID is present in the current collection
        requested_ids_set = set(ids)
        if not current_ids.intersection(requested_ids_set):
            logging.warning(
                f"No requested IDs {list(ids)} found in current collection. "
                "Returning empty FeatureCollection."
            )
            return FeatureCollection([])

        filtered_features = [feature for feature in self._features if hasattr(feature, 'ID') and feature.ID in ids]
        return FeatureCollection(filtered_features)

    def add_feature(self, feature: Feature):
        """Adds a single Feature object to the collection."""
        self._features.append(feature)

    def extend_features(self, features: Iterable[Feature]):
        """Adds multiple Feature objects to the collection."""
        self._features.extend(features)

    def subset_lonlat(self, extent: Tuple[float, float, float, float]) -> 'FeatureCollection':
        """
        Subsets the FeatureCollection to include only features whose locations
        fall within the specified longitude and latitude extent.
        Longitudes of features are normalized to -180 to 180 for the check.

        Parameters
        ----------
        extent : Tuple[float, float, float, float]
            A tuple representing (lon0, lon1, lat0, lat1) for the bounding box.
            Longitudes should be in the -180 to 180 range.

        Returns
        -------
        FeatureCollection
            A new FeatureCollection containing only the features within the extent.
        """
        lon0, lon1, lat0, lat1 = extent

        # Ensure extent longitudes are in the -180 to 180 range for consistency
        # (though the method itself handles conversion of feature longitudes)
        lon0_norm = ((lon0 + 180) % 360) - 180
        lon1_norm = ((lon1 + 180) % 360) - 180

        filtered_features: List[Feature] = []
        for feature in self._features:
            # The feature's own subset_lonlat method handles its longitude normalization
            if feature.subset_lonlat((lon0_norm, lon1_norm, lat0, lat1)) is not None:
                filtered_features.append(feature)

        logging.info(f"Subsetted FeatureCollection to extent {extent}. New count: {len(filtered_features)}")
        return FeatureCollection(filtered_features)


# Ignore specific UserWarnings from Cartopy, which can be noisy during plotting.
warnings.filterwarnings("ignore", category=UserWarning, module="cartopy")


# --- Changed functions ---

def label_id(ax: plt.Axes, features_input: Union[pd.DataFrame, Feature, FeatureCollection], **kwargs: Any) -> List[plt.Artist]:
    """
    Plot feature locations, draw lines from initial ('V') locations to current
    locations, and add ID markers.

    Args:
        ax (plt.Axes): The matplotlib axes object to plot on.
        features_input (Union[pd.DataFrame, Feature, FeatureCollection]):
            Can be a DataFrame, a single Feature, or a FeatureCollection.
            Will be converted to a DataFrame internally.
        **kwargs: Additional keyword arguments to pass to the plot functions,
                  e.g., 'markersize', 'linewidth'.

    Returns:
        List[plt.Artist]: A list of matplotlib artists (plot objects) created.
    """
    if isinstance(features_input, Feature):
        df = features_input.to_dataframe()
    elif isinstance(features_input, FeatureCollection):
        df = features_input.to_dataframe()
    elif isinstance(features_input, pd.DataFrame):
        df = features_input
    else:
        raise TypeError("features_input must be a DataFrame, Feature, or FeatureCollection.")

    markers: List[str] = [">", "<", "P", "D", "^", "v", "*", "s", "X", "d"]
    ts: List[plt.Artist] = []
    markersize: int = kwargs.pop("markersize", 10)
    linewidth: int = kwargs.pop("linewidth", 1)

    for _, row in df.iterrows():
        x = row["LON(E)"]
        y = row["LAT(N)"]
        # Ensure color column exists and is not NaN
        color = row.get("color", "gray") # Default to gray if 'color' column is missing or NaN
        if pd.isna(color):
            color = "gray"

        if "VLon(E)" in row and pd.notna(row["VLon(E)"]):
            vx = row["VLon(E)"]
            vy = row["VLat(N)"]
        else:
            # If 'V' location is not available (e.g., for analysis), use current location
            vx = x
            vy = y

        # Plot line and marker for ID
        logging.debug(f"marker and line from ({vx:.2f},{vy:.2f}) to ({x:.2f},{y:.2f}) for ID {row['ID']}")
        marker: str = markers[int(row["ID"] % len(markers))]
        (id_plot,) = ax.plot(
            [vx, x],
            [vy, y],
            marker=marker,
            markersize=markersize,
            linewidth=linewidth,
            color="k", # Line color black
            markeredgecolor="none",
            alpha=0.5,
            label=int(row["ID"]),
            transform=ccrs.Geodetic(), # Use Geodetic for straight lines on the globe
        )
        ts.append(id_plot)
        ts.append(
            ax.annotate(
                "",
                xytext=(vx, vy),
                xy=(x, y),
                arrowprops=dict(arrowstyle="->", color="k", linewidth=linewidth),
                transform=ccrs.PlateCarree(),  # because lon/lat are in geographic coords
            )
        )

    return ts

def tissot(ax: plt.Axes, features_input: Union[pd.DataFrame, Feature, FeatureCollection], **kwargs: Any) -> List[plt.Artist]:
    """
    Draw circles (Tissot's indicatrices) with 'Ro(km)' radius around
    (LON(E), LAT(N)) locations on a Cartopy axes.

    Parameters
    ----------
    ax : plt.Axes
        The matplotlib axes object to draw on.
    features_input : Union[pd.DataFrame, Feature, FeatureCollection]
        Can be a DataFrame, a single Feature, or a FeatureCollection.
        Will be converted to a DataFrame internally.
    **kwargs : Any
        Additional keyword arguments to pass to `ax.tissot()`,
        e.g., 'alpha', 'edgecolor', 'facecolor'.
        If 'color' is provided in kwargs, it overrides any 'color' column in df.

    Returns
    -------
    List[plt.Artist]
        A list of matplotlib artists (patches) created by `ax.tissot()`.
    """
    if isinstance(features_input, Feature):
        df = features_input.to_dataframe()
    elif isinstance(features_input, FeatureCollection):
        df = features_input.to_dataframe()
    elif isinstance(features_input, pd.DataFrame):
        df = features_input
    else:
        raise TypeError("features_input must be a DataFrame, Feature, or FeatureCollection.")

    ts: List[plt.Artist] = []
    if df.empty:
        return ts

    # If color is given in kwargs, override or create a 'color' column in the temporary DataFrame
    temp_df = df.copy()
    color_override = kwargs.pop("color", None)
    if color_override:
        temp_df["color"] = color_override
    elif "color" not in temp_df.columns:
        # Default color if not provided and not in DataFrame
        temp_df["color"] = "blue"

    for _, row in temp_df.iterrows():
        # Ensure Ro(km) is not NaN for tissot and is positive
        if pd.isna(row.get("Ro(km)")) or row["Ro(km)"] <= 0: # Use .get for robustness
            logging.debug(f"Skipping tissot for feature ID {row.get('ID', 'N/A')} due to invalid Ro(km): {row.get('Ro(km)')}")
            continue

        t = ax.tissot(
            rad_km=row["Ro(km)"],
            lons=row["LON(E)"],
            lats=row["LAT(N)"],
            color=row["color"],
            **kwargs,
        )
        ts.append(t)
    return ts


def animate(
    valid_time: Union[pd.Timestamp, np.datetime64, xarray.DataArray],
    ax: plt.Axes,
    itime: Union[pd.Timestamp, np.datetime64, xarray.DataArray],
    workdir: Path,
    forecast_length: int,
    ids: Optional[Iterable[int]] = None,
    isensemble: bool = False,
) -> Tuple[plt.Axes]:
    """
    Generates a single frame for an animation, plotting observed and forecast features.

    Parameters
    ----------
    valid_time : Union[pd.Timestamp, np.datetime64, xarray.DataArray]
        The valid time for the current animation frame.
    ax : plt.Axes
        The matplotlib axes object to draw on.
    itime : Union[pd.Timestamp, np.datetime64, xarray.DataArray]
        The initialization time of the forecast.
    workdir : Path
        The base working directory where forecast output is stored.
    forecast_length : int
        The forecast length in hours.
    ids : Optional[Iterable[int]], optional
        An iterable of feature IDs to filter and plot, by default None (plot all).
    isensemble : bool, optional
        True if plotting ensemble members, False for deterministic forecast, by default False.

    Returns
    -------
    Tuple[plt.Axes]
        A tuple containing the updated matplotlib Axes object.
    """
    logging.info(f"Animating for valid_time: {valid_time}")

    # Convert timestamps to pandas.Timestamp for consistent handling
    valid_time_ts = handleTimestamp(valid_time)
    itime_ts = handleTimestamp(itime)

    # Clear previous annotations to update the plot
    for coll in getattr(ax, "_annotations", []):
        try:
            if hasattr(coll, 'remove'):
                coll.remove()
                logging.debug(f"Removed old annotation: {coll}")
        except ValueError:
            # Handle cases where the artist might have already been removed by other means
            logging.debug(f"Could not remove {coll}, possibly already removed or not removable.")
    ax._annotations = []

    alpha = 0.8
    fhr = (valid_time_ts - itime_ts) / pd.Timedelta(hours=1)

    workdirs: List[Path] = [workdir / f"{itime_ts.strftime(fmt)}.F{forecast_length:03d}.C768"]
    if isensemble:  # all the members
        workdirs = list(workdir.glob(f"E{itime_ts.strftime(fmt)}.p??.*"))
        alpha = 0.2

    text_kw: Dict[str, Any] = dict(
        fontsize="xx-small",
        horizontalalignment="center",
        verticalalignment="center",
        transform=cartopy.crs.PlateCarree(),
        clip_on=True,
    )

    # --- Handle Observed Features ---
    # getobs now returns Union[Feature, FeatureCollection]
    obs_result: Union[Feature, FeatureCollection] = getobs(valid_time_ts, ids=ids)

    obs_result = obs_result.subset_lonlat(ax.get_extent(crs=ccrs.PlateCarree()))
    
    # Pass obs_result directly to tissot and label_id, they will handle conversion
    if obs_result: # Check if it's not an empty FeatureCollection or None
        ax._annotations.extend(tissot(ax, obs_result, alpha=0.4, color='k', facecolor="none"))
        ax._annotations.extend(label_id(ax, obs_result, **text_kw))


    # --- Handle Forecast Features ---
    for current_workdir in workdirs:
        logging.info(f"Processing forecast for workdir: {current_workdir.name}")
        # getfcst now returns Union[Feature, FeatureCollection]
        fcst_result: Union[Feature, FeatureCollection] = getfcst(itime_ts, valid_time_ts, current_workdir.parent,
                                                 isensemble=isensemble, ids=ids)

        ax.set_title(f"{valid_time_ts.strftime(fmt)} f{fhr:03.0f}")
        fcst_result = fcst_result.subset_lonlat(ax.get_extent(crs=ccrs.PlateCarree()))

        if fcst_result: # Check if it's not an empty FeatureCollection or None
            ax._annotations.extend(tissot(ax, fcst_result, alpha=alpha))

            # For label_id, we need to filter to visible features, which means converting to DataFrame first
            # The label_id function is now updated to take Union[pd.DataFrame, Feature, FeatureCollection]
            # but the filtering logic here still expects a DataFrame.
            # So, we convert the result to DataFrame and then apply the filter.
            fcst_df_for_labels = None
            if isinstance(fcst_result, Feature):
                fcst_df_for_labels = fcst_result.to_dataframe()
            elif isinstance(fcst_result, FeatureCollection):
                fcst_df_for_labels = fcst_result.to_dataframe()

            if fcst_df_for_labels is not None and not fcst_df_for_labels.empty:
                visible_df = fcst_df_for_labels[~fcst_df_for_labels.ID.isnull()]
                logging.info(f"{len(visible_df)} with ID")

                ax._annotations.extend(label_id(
                    ax,
                    visible_df, # Pass the filtered DataFrame
                    **text_kw,
                ))

    # Don't show duplicate labels in the legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='best')

    return (ax,)


def handleTimestamp(timestamp: Union[pd.Timestamp, np.datetime64, xarray.DataArray]) -> pd.Timestamp:
    """
    Converts various timestamp formats to a pandas.Timestamp object.

    Parameters
    ----------
    timestamp : Union[pd.Timestamp, np.datetime64, xarray.DataArray]
        The timestamp to convert. Can be a pandas Timestamp, numpy datetime64,
        or an xarray DataArray containing datetime64 values.

    Returns
    -------
    pd.Timestamp
        The converted pandas Timestamp object.

    Raises
    ------
    TypeError
        If the provided timestamp type is not supported.
    """
    # Case 1: If it's already a pandas Timestamp
    if isinstance(timestamp, pd.Timestamp):
        return timestamp
    # Case 2: If it's a numpy datetime64
    elif isinstance(timestamp, np.datetime64):
        return pd.to_datetime(timestamp)
    # Case 3: If it's an xarray DataArray (assuming it contains datetime64 values)
    elif isinstance(timestamp, xarray.DataArray):
        # If it's a single value DataArray, extract the scalar
        if timestamp.ndim == 0:
            return pd.to_datetime(timestamp.item())
        # If it's an array, convert the first element (or handle as needed for multiple times)
        # For this context (single time for a forecast/obs), assuming it's a scalar or 1-element array
        return pd.to_datetime(timestamp.values[0]) # Use [0] to get the scalar from a 0-dim array
    else:
        raise TypeError(f"Unsupported type for timestamp argument: {type(timestamp)}")


def getfcst(
    itime: Union[pd.Timestamp, np.datetime64, xarray.DataArray],
    valid_time: Union[pd.Timestamp, np.datetime64, xarray.DataArray],
    workdir: Path,
    isensemble: bool = False,
    ids: Optional[Iterable[int]] = None,
) -> Union[Feature, 'FeatureCollection']: # Type hint for FeatureCollection needs to be string if FeatureCollection is defined later
    """
    Loads forecast track data from CSV files and returns a FeatureCollection,
    or a single Feature if only one is found after filtering.

    Parameters
    ----------
    itime : Union[pd.Timestamp, np.datetime64, xarray.DataArray]
        The initialization time of the forecast.
    valid_time : Union[pd.Timestamp, np.datetime64, xarray.DataArray]
        The valid time for which to retrieve the forecast data.
    workdir : Path
        The base working directory where forecast output is stored.
    isensemble : bool, optional
        True if loading ensemble members, False for deterministic forecast, by default False.
    ids : Optional[Iterable[int]], optional
        An iterable of feature IDs to filter, by default None (return all).

    Returns
    -------
    Union[Feature, FeatureCollection]
        A collection of Feature objects, or a single Feature if only one is found.

    Raises
    ------
    AssertionError
        If `valid_time` is before `itime`.
    """
    itime_ts: pd.Timestamp = handleTimestamp(itime)
    valid_time_ts: pd.Timestamp = handleTimestamp(valid_time)
    assert itime_ts <= valid_time_ts, f"valid_time {valid_time_ts} before initialization time {itime_ts}"

    # Determine which forecast directories to look into
    workdirs: List[Path] = [workdir / itime_ts.strftime(fmt)]
    if isensemble:  # all the members
        # Adjust glob pattern to match ensemble member directories (e.g., E2020010100.p01.C768)
        workdirs = list(workdir.glob(f"E{itime_ts.strftime(fmt)}.p??.*"))
        if not workdirs:
            logging.warning(f"No ensemble member directories found for {itime_ts.strftime(fmt)} in {workdir}")
            return FeatureCollection([])

    fhr: float = (valid_time_ts - itime_ts) / pd.Timedelta(hours=1)

    all_fcst_dfs: List[pd.DataFrame] = []
    fname: str = (
        f"gfs.0p25.{itime_ts.strftime(fmt)}.f000.track"
        if itime_ts == valid_time_ts
        else f"diag_TroughsCutoffs.{itime_ts.strftime(fmt)}.f{fhr:03.0f}.track"
    )

    for current_workdir in workdirs:
        ifile: Path = current_workdir / fname
        if not ifile.exists():
            logging.warning(f"Forecast track file not found: {ifile}")
            continue # Skip to the next member/directory

        try:
            df_member = pd.read_csv(
                ifile,
                header=0,
                sep=r"\s+",
                na_values=na_values,
            )
            all_fcst_dfs.append(df_member)
        except pd.errors.EmptyDataError:
            logging.warning(f"Empty data or no columns in file: {ifile}")
            continue
        except Exception as e:
            logging.error(f"Error reading file {ifile}: {e}")
            continue

    if not all_fcst_dfs:
        logging.warning(f"No forecast data loaded for {itime_ts} at valid time {valid_time_ts}.")
        return FeatureCollection([])

    fcst_combined = pd.concat(all_fcst_dfs, ignore_index=True)

    fcst_combined["location"] = get_location(fcst_combined)
    fcst_combined["color"] = fcst_combined["location"].map(loc2color)
    fcst_combined["contains_zmin"] = get_contains_zmin(fcst_combined)
    fcst_combined["valid_time"] = valid_time_ts

    # Convert DataFrame rows to Feature objects
    features: List[Feature] = [Feature(row) for _, row in fcst_combined.iterrows()]
    feature_collection = FeatureCollection(features)

    if ids:
        feature_collection = feature_collection.filter_by_id(ids)

    # Return Feature if only one, otherwise FeatureCollection
    if len(feature_collection) == 1:
        return feature_collection[0]
    return feature_collection


def getobs(
    valid_time: Union[pd.Timestamp, np.datetime64, xarray.DataArray],
    ids: Optional[Iterable[int]] = None,
) -> Union[Feature, 'FeatureCollection']: # Type hint for FeatureCollection needs to be string if FeatureCollection is defined later
    """
    Loads observational track data from a CSV file and returns a FeatureCollection,
    or a single Feature if only one is found after filtering.

    Parameters
    ----------
    valid_time : Union[pd.Timestamp, np.datetime64, xarray.DataArray]
        The valid time for which to retrieve the observation data.
    ids : Optional[Iterable[int]], optional
        An iterable of feature IDs to filter, by default None (return all).

    Returns
    -------
    Union[Feature, FeatureCollection]
        A collection of Feature objects, or a single Feature if only one is found.

    Raises
    ------
    FileNotFoundError
        If the observation track file cannot be found.
    """
    valid_time_ts: pd.Timestamp = handleTimestamp(valid_time)
    obs_path: Path = Path(valid_time_ts.strftime(
        "/glade/u/home/klupo/work_new/postdoc/kasugaEA21/version9/HGT_500mb/"
        f"gfs.0p25.%Y%m%d%H.f000.track"
    ))

    if not obs_path.exists():
        logging.warning(f"Observation track file not found: {obs_path}. Returning empty collection.")
        return FeatureCollection([])

    try:
        obs_df = pd.read_csv(
            obs_path,
            header=0,
            sep=r"\s+",
            na_values=na_values,
        )
    except pd.errors.EmptyDataError:
        logging.warning(f"Empty data or no columns in observation file: {obs_path}")
        return FeatureCollection([])
    except Exception as e:
        logging.error(f"Error reading observation file {obs_path}: {e}")
        return FeatureCollection([])

    obs_df["color"] = "k"  # Observations typically colored black
    obs_df["contains_zmin"] = get_contains_zmin(obs_df)
    obs_df["valid_time"] = valid_time_ts

    # Convert DataFrame rows to Feature objects
    features: List[Feature] = [Feature(row) for _, row in obs_df.iterrows()]
    feature_collection = FeatureCollection(features)

    if ids:
        feature_collection = feature_collection.filter_by_id(ids)

    # Return Feature if only one, otherwise FeatureCollection
    if len(feature_collection) == 1:
        return feature_collection[0]
    return feature_collection


def rotate_field_xarray(
    field: xarray.DataArray,
    lat_a: float,
    lon_a: float,
    lat_b: float,
    lon_b: float,
) -> xarray.DataArray:
    """
    Rotate a field defined on a set of lat/lon points and interpolate it back to the original grid.
    This is useful for analyzing fields relative to a moving feature by effectively "recentering"
    the field.

    Parameters
    ----------
    field : xarray.DataArray
        Field values with 2D coordinates `XLAT` (latitude) and `XLONG` (longitude),
        and optionally a `valid_time` dimension.
    lat_a : float
        Latitude of the original reference point A.
    lon_a : float
        Longitude of the original reference point A.
    lat_b : float
        Latitude of the new position B.
    lon_b : float
        Longitude of the new position B.

    Returns
    -------
    xarray.DataArray
        Rotated and interpolated field with the same dimensions and coordinates as `field`.
    """
    if (lat_a, lon_a) == (lat_b, lon_b):
        logging.info("Start and end locations are equal. Returning the original array without rotation.")
        return field

    def latlon_to_cartesian(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
        """Vectorized conversion from lat/lon to 3D Cartesian coordinates."""
        lat_rad = np.radians(lat)
        lon_rad = np.radians(lon)
        x = np.cos(lat_rad) * np.cos(lon_rad)
        y = np.cos(lat_rad) * np.sin(lon_rad)
        z = np.sin(lat_rad)
        return np.stack([x, y, z], axis=-1)

    def cartesian_to_latlon(cartesian: np.ndarray) -> np.ndarray:
        """Vectorized conversion from 3D Cartesian coordinates to lat/lon."""
        x, y, z = cartesian.T
        lat = np.degrees(np.arcsin(z))
        lon = np.degrees(np.arctan2(y, x))
        return np.stack([lat, lon], axis=-1)

    def compute_rotation_matrix(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
        """Computes a 3x3 rotation matrix to rotate vector v1 to vector v2."""
        v1_norm = v1 / np.linalg.norm(v1)
        v2_norm = v2 / np.linalg.norm(v2)
        axis = np.cross(v1_norm, v2_norm)
        axis_norm = np.linalg.norm(axis)

        if axis_norm == 0:  # v1 and v2 are collinear
            return np.eye(3)
        axis = axis / axis_norm
        angle = np.arccos(np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)) # Ensure angle is valid

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

    # Ensure 'XLAT' and 'XLONG' exist as coordinates
    if "XLAT" not in field.coords or "XLONG" not in field.coords:
        # If 'XLAT' and 'XLONG' are data variables, make them coordinates
        if "XLAT" in field.data_vars and "XLONG" in field.data_vars:
            field = field.set_coords(["XLAT", "XLONG"])
        else:
            raise ValueError("Field must have 'XLAT' and 'XLONG' as coordinates or data variables.")


    # Extract latitude and longitude arrays
    lat_grid: np.ndarray = field["XLAT"].values
    lon_grid: np.ndarray = field["XLONG"].values
    lat_flat = lat_grid.ravel()
    lon_flat = lon_grid.ravel()

    # Calculate rotation
    A = latlon_to_cartesian(np.array([lat_a]), np.array([lon_a]))[0] # Convert scalar to array for function
    B = latlon_to_cartesian(np.array([lat_b]), np.array([lon_b]))[0]
    R = compute_rotation_matrix(A, B)

    # Rotate all points in one operation
    S_cartesian = latlon_to_cartesian(lat_flat, lon_flat)
    S_prime_cartesian = S_cartesian @ R.T
    S_prime = cartesian_to_latlon(S_prime_cartesian)
    lat_prime, lon_prime = S_prime[:, 0], S_prime[:, 1]

    def interpolate_field(field_slice: np.ndarray) -> np.ndarray:
        """Helper function for interpolation (used by apply_ufunc or directly)."""
        # griddata expects points as (y, x) or (lat, lon)
        return griddata(
            points=(lat_prime, lon_prime), # rotated source points
            values=field_slice.ravel(),    # values at rotated source points
            xi=(lat_flat, lon_flat),       # original target grid points
            method="linear",
        ).reshape(lat_grid.shape)

    # Handle interpolation for time-dependent or static fields
    if "valid_time" in field.dims:
        # Apply the interpolation to each time slice (if 'valid_time' is a dimension)
        interpolated_data = xarray.apply_ufunc(
            interpolate_field,
            field,
            input_core_dims=[["lat", "lon"]], # assuming 'lat' and 'lon' are the spatial dims
            output_core_dims=[["lat", "lon"]],
            exclude_dims=set(("lat", "lon")), # Exclude original spatial dims from iterating
            dask="parallelized",
            output_dtypes=[field.dtype],
        )
        interpolated_field = field.copy(data=interpolated_data) # Retain original coords and attrs
    else:
        # If no 'valid_time' dimension, apply interpolation directly
        interpolated_data = interpolate_field(field.values)
        interpolated_field = xarray.DataArray(
            interpolated_data,
            dims=field.dims, # Use original dimensions
            coords=field.coords, # Use original coordinates
            attrs=field.attrs,
        )

    # Update metadata to reflect the rotation
    interpolated_field.attrs.update(
        dict(rotated_from_lat=lat_a, rotated_from_lon=lon_a,
             rotated_to_lat=lat_b, rotated_to_lon=lon_b))
    return interpolated_field


def stack(ds: xarray.Dataset) -> xarray.Dataset:
    """
    Groups same variable from different pressure levels into new variables
    with a 'pfull' vertical dimension.

    Parameters
    ----------
    ds : xarray.Dataset
        The input xarray Dataset containing variables at discrete pressure levels
        (e.g., 't1000', 't850').

    Returns
    -------
    xarray.Dataset
        A new xarray Dataset with vertically stacked variables (e.g., 't' with a 'pfull' dimension).
    """
    potential_plevels: List[int] = [1000, 925, 850, 700, 600,
                                    500, 400, 300, 250, 200, 150, 100, 50, 10]

    # Create a copy to modify
    ds_modified = ds.copy()

    for statevar in ["omg", "q", "rh", "t", "u", "v", "z"]:
        # Handle state variables themselves (e.g., 't', 'q')
        potentialvars_state = [f"{statevar}{plev}" for plev in potential_plevels]
        this_vars_state = [var for var in potentialvars_state if var in ds_modified]
        if this_vars_state:
            stack_list: List[xarray.DataArray] = []
            for var in this_vars_state:
                plev = int(var.removeprefix(statevar))
                tmpvar = ds_modified[var].expand_dims(pfull=[plev])
                tmpvar.attrs["description"] = tmpvar.attrs.get("description", "").replace(f"{plev}-mb ", "")
                stack_list.append(tmpvar)
                ds_modified = ds_modified.drop_vars(var)
            ds_modified[statevar] = xarray.concat(stack_list, dim="pfull")
            logging.info(f"Stacked variable: {statevar}")

        # Handle tendency variables (e.g., 'dtemp3dt_mp', 'dtemp3dt_cnvgwd')
        for tend in [
            "cnvgwd", "deepcnv", "lw", "mp", "nonphys", "ogwd", "pbl", "rdamp", "shalcnv", "sw",
        ]:
            potentialvars_tend = [f"d{statevar}3dt{plev}_{tend}" for plev in potential_plevels]
            newvar_tend = f"d{statevar}3dt_{tend}"
            this_vars_tend = [var for var in potentialvars_tend if var in ds_modified]
            if not this_vars_tend:
                continue

            stack_list_tend: List[xarray.DataArray] = []
            for var in this_vars_tend:
                plev_str = var.removeprefix(f"d{statevar}3dt").removesuffix(f"_{tend}")
                try:
                    plev = int(plev_str)
                except ValueError:
                    logging.warning(f"Could not parse pressure level from variable name: {var}")
                    continue

                tmpvar = ds_modified[var].expand_dims(pfull=[plev])
                tmpvar.attrs["description"] = tmpvar.attrs.get("description", "").replace(f"{plev}-mb ", "")
                stack_list_tend.append(tmpvar)
                ds_modified = ds_modified.drop_vars(var)

            if stack_list_tend: # Only add if there were variables to stack
                ds_modified[newvar_tend] = xarray.concat(stack_list_tend, dim="pfull")
                logging.info(f"Stacked tendency variable: {newvar_tend}")
    return ds_modified
