import json
import os
import requests
import geopandas as gpd
from shapely.geometry import Polygon
import h3
from h3 import LatLngPoly
import shapely
from shapely.ops import transform
from shapely.geometry import Polygon
from pyproj import Transformer

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────
def to_2d(poly):
    """Return a 2‑D copy of a Shapely polygon (drops any Z values)."""
    if not poly.has_z:
        return poly
    exterior = [(x, y) for x, y, *_ in poly.exterior.coords]
    interiors = [
        [(x, y) for x, y, *_ in ring.coords] for ring in poly.interiors
    ]
    return Polygon(exterior, interiors)


def import_shapefile_to_shapely(path):
    """
    Load shapefile and return the unified geometry in EPSG:4326
    """
    gdf = gpd.read_file(path)

    if gdf.crs != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")

    # Dissolve all polygons into one geometry (if multiple)
    unified_geom = gdf.unary_union
    return unified_geom

def polygon_to_h3(unified_geom, resolution=9):
    h3_indices = set()
    shapelymapping = shapely.geometry.mapping(unified_geom)
    if shapelymapping['type'] == 'Polygon':
        # geo_json_geom already conforms to what h3 expects
        geo_geom = shapely.geometry.Polygon(unified_geom)
        lat_lng_coords = [(coord[1], coord[0]) for coord in geo_geom.exterior.coords]
        poly = LatLngPoly(lat_lng_coords)
        h3_indices.update(h3.h3shape_to_cells(poly, resolution))


    elif shapelymapping['type'] == 'MultiPolygon':
        # geo_json_geom already conforms to what h3 expects
        geo_geom = shapely.geometry.MultiPolygon(unified_geom)
        polygons = geo_geom.geoms
        h3_polys = set()
        for poly in polygons:
            lat_lng_coords = [(coord[1], coord[0]) for coord in geo_geom.exterior.coords]
            h3_poly = LatLngPoly(lat_lng_coords)
            h3_indices.update(h3.h3shape_to_cells(h3_poly, resolution))

    else:
        raise ValueError(f"Unsupported geometry type: {shapelymapping['type']}")

    return list(h3_indices)


def h3_to_geo(h3_id):
    geo_cells_h3 = h3.cells_to_geo(h3_id, tight=False)
    return geo_cells_h3

def save_h3_polygons(polygons, output_path='C:\\data_engineering\\pdal\\h3_polygons', output_filename="h3_tiles.geojson"):
    """
    Save the list of shapely polygons to a GeoJSON file
    """
    h3_path = os.path.join(output_path, output_filename)
    gdf = gpd.GeoDataFrame(geometry=polygons, crs="EPSG:4326")
    gdf.to_file(h3_path, driver="GeoJSON")

def shapely_to_esri_rings(polygon):
    def flatten_coords(coords):
        return [[coord[0], coord[1]] for coord in coords]
    rings = flatten_coords(polygon.exterior.coords)
    for interior in polygon.interiors:
        rings.append(flatten_coords(interior.coords))
    return rings

def get_usgs_datasets(user_aoi, layer_id):
    """
    Intersect AOI with the specified 3DEP Elevation Index layer and return GeoJSON features.

    Args:
        user_aoi (shapely.geometry.Polygon): The user's area of interest as a Shapely Polygon.
                                             Assumed to be in EPSG:4326 unless layer_id is 7.
        layer_id (int): The ID of the USGS 3DEP Elevation Index layer to query.
    """
    url = f"https://index.nationalmap.gov/arcgis/rest/services/3DEPElevationIndex/MapServer/{layer_id}/query"
    target_srid = "4326"
    aoi_transformed = user_aoi

    if layer_id == 7:
        target_srid = "3857"
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
        aoi_transformed = transform(transformer.transform, user_aoi)

    minx, miny, maxx, maxy = aoi_transformed.bounds
    geometry_str = f"{minx},{miny},{maxx},{maxy}"
    geometry_rings = shapely_to_esri_rings(aoi_transformed)

    # Envelope query
    params_envelope = {
        "where": "1=1",
        "geometry": geometry_str,
        "geometryType": "esriGeometryEnvelope",
        "inSR": target_srid,
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "*",
        "returnGeometry": "true",
        "f": "geojson",
    }

    # Polygon query
    polygon_geometry = {
        "rings": geometry_rings,
        "spatialReference": {"wkid": int(target_srid)}
    }

    params_polygon = {
        "where": "1=1",
        "geometry": json.dumps(polygon_geometry),
        "geometryType": "esriGeometryPolygon",
        "inSR": target_srid,
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "*",
        "returnGeometry": "true",
        "f": "geojson",
    }

    # Use GET for envelope
    features_envelope = requests.get(url, params=params_envelope)

    # Use POST for polygon
    # features_ring = requests.post(url, data=params_polygon)

    return features_envelope.json()