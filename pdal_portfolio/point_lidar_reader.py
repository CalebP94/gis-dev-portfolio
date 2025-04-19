"""
Updated 3DEP download + PDAL processing script
---------------------------------------------

•  Clips every point cloud to the exact AOI polygon (Z‑values stripped)
•  Runs the pipeline immediately via pdal.Pipeline().execute()
"""

import json
import os
import boto3
import requests
import pdal
import geopandas as gpd
import pyproj
from shapely.geometry import Polygon, shape
import h3
from h3 import LatLngPoly
import shapely
from shapely.ops import transform
from shapely.geometry import Polygon, MultiPolygon
from pyproj import Transformer
from shapely.geometry.polygon import orient
from pathlib import Path

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

def polygon_to_h3(unified_geom, resolution=11):
    h3_indices = set()
    geom_oriented = orient(unified_geom, sign=1.0)

    def coords_to_latlng(coords):
        x, y = coords.xy
        if len(x) >= 3:
            return [(y[i], x[i]) for i in range(len(x))]
        return None

    if isinstance(unified_geom, Polygon):
        outer = coords_to_latlng(geom_oriented.exterior.coords)
        if not outer:
            raise ValueError("Outer ring must have at least 3 points")

        poly = LatLngPoly(outer)
        h3_indices.update(h3.h3shape_to_cells(poly, resolution))

    elif isinstance(unified_geom, MultiPolygon):
        for part in unified_geom.geoms:
            outer = coords_to_latlng(part.exterior.coords)
            if not outer:
                continue  # skip bad shells

            poly = LatLngPoly([outer])
            h3_indices.update(h3.h3shape_to_cells(poly, resolution))

    else:
        raise ValueError("Expected Polygon or MultiPolygon")

    return list(h3_indices)


def h3_cell_production(h3_ids, output_path=r"C:\data_engineering\pdal\hexes_individual.shp"):

    features = []
    for h in h3_ids:
        ring = [(lon, lat) for lat, lon in h3.cell_to_boundary(h)]
        poly = Polygon(ring)

        if not poly.is_valid:    
            poly = poly.buffer(0)

        features.append({"geometry": poly, "h3_id": h})

    gdf = gpd.GeoDataFrame(features, crs="EPSG:4326")

    gdf.to_file(output_path, driver="ESRI Shapefile")   # or "ESRI Shapefile" if you must

    return gdf


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

    # Use GET for envelope
    features_envelope = requests.get(url, params=params_envelope)

    # Use POST for polygon
    # features_ring = requests.post(url, data=params_polygon)

    return features_envelope.json()



# ──────────────────────────────────────────────────────────────────────────────
# PDAL‑pipeline builders
# ──────────────────────────────────────────────────────────────────────────────
def build_pdal_pipeline(
    extent_epsg4326,
    usgs_3dep_dataset_names,
    pc_resolution,
    *,
    filterNoise=False,
    reclassify=False,
    savePointCloud=False,
    outCRS=4326,
    pc_outName="pc",
    pc_outType="laz",
):
    """Return a PDAL pipeline dict that fetches and crops point clouds."""
    aoi2d = to_2d(extent_epsg4326)
    catalog_url = "https://usgs-lidar-stac.s3-us-west-2.amazonaws.com/ept/catalog.json"
    response = requests.get(catalog_url)
    catalog = response.json()

    # Create a lowercase-to-original mapping
    ept_map = {}

    for link in catalog.get("links", []):
        href = link.get("href", "")
        if href.endswith(".json") and "/ept/" in href:
            original = href.split("/")[-1].replace(".json", "")
            lowercase = original.lower()
            lowercase = lowercase.replace("-", "_")
            ept_map[lowercase] = original

    stages = []
    for name in usgs_3dep_dataset_names:
        lower_name = name.lower()
        if lower_name in ept_map:
            name = ept_map[lower_name]
        else:
            return None
        url = f"https://s3-us-west-2.amazonaws.com/usgs-lidar-public/{name}/ept.json"
        #use to understand the size of the point cloud, but causes slowdown (significant)
        # point_count = get_point_count(extent_epsg4326=aoi2d, workunit=name)
        
        r = requests.head(url, allow_redirects=True)
        if r.status_code != 200:
            return None
        stages.append(
            {
                "type": "readers.ept",
                "filename": f"{url}",
                "resolution": pc_resolution,
            }
        )

    # crop to AOI polygon
    stages.append({"type": "filters.crop", "polygon": aoi2d.wkt})

    # optional noise filters
    if filterNoise:
        stages.extend(
            [
                {"type": "filters.range", "limits": "Classification![7:7]"},
                {"type": "filters.range", "limits": "Classification![18:18]"},
            ]
        )

    # optional re‑classification to bare‑earth
    if reclassify:
        stages.extend(
            [
                {"type": "filters.assign", "value": "Classification = 0"},
                {"type": "filters.smrf"},
                {"type": "filters.range", "limits": "Classification[2:2]"},
            ]
        )

    # reprojection
    stages.append({"type": "filters.reprojection", "out_srs": f"EPSG:{outCRS}"})

    # optional writer
    if savePointCloud:
        if pc_outType not in {"las", "laz"}:
            raise ValueError("pc_outType must be 'las' or 'laz'")
        writer = {
            "type": "writers.las",
            "filename": f"{pc_outName}.{pc_outType}",
        }
        if pc_outType == "laz":
            writer["compression"] = "laszip"
        stages.append(writer)

    return {"pipeline": stages}


def make_DEM_pipeline(
    extent_epsg4326,
    usgs_3dep_dataset_name,
    pc_resolution,
    dem_resolution,
    *,
    filterNoise=True,
    reclassify=False,
    savePointCloud=False,
    outCRS=4326,
    pc_outName="pc",
    pc_outType="laz",
    demType="dtm",
    gridMethod="idw",
    dem_outName="dem",
    dem_outExt="tif",
    driver="GTiff",
):
    """Return a PDAL pipeline dict that outputs a DTM/DSM raster."""
    pipe = build_pdal_pipeline(
        extent_epsg4326,
        usgs_3dep_dataset_name,
        pc_resolution,
        filterNoise=filterNoise,
        reclassify=reclassify,
        savePointCloud=savePointCloud,
        outCRS=outCRS,
        pc_outName=pc_outName,
        pc_outType=pc_outType,
    )
    if not pipe:
        return None

    
    if demType == "dtm":
        pipe["pipeline"].append(
            {"type": "filters.range", "limits": "Classification[2:2]"}
        )
    elif demType != "dsm":
        raise ValueError("demType must be 'dtm' or 'dsm'")

    pipe["pipeline"].append(
        {
            "type": "writers.gdal",
            "filename": f"{dem_outName}.{dem_outExt}",
            "gdaldriver": driver,
            "nodata": -9999,
            "output_type": gridMethod,
            "resolution": float(dem_resolution),
            "gdalopts": (
                "COMPRESS=LZW,TILED=YES,blockxsize=256,"
                "blockysize=256,COPY_SRC_OVERVIEWS=YES"
            ),
        }
    )
    return pipe

def get_point_count(extent_epsg4326, workunit):

    url = f"https://s3-us-west-2.amazonaws.com/usgs-lidar-public/{workunit}/ept.json"
    aoi2d = to_2d(extent_epsg4326)

    pipeline_dict = {
        "pipeline": [
            {
                "type": "readers.ept",
                "filename": url,
                "resolution": 2.0
            },
            {
                "type": "filters.crop",
                "polygon": aoi2d.wkt
            },
            {
                "type": "filters.stats"
            }
        ]
    }

    try:
        pipeline = pdal.Pipeline(json.dumps(pipeline_dict))
        pipeline.execute()
        metadata = json.loads(pipeline.metadata)
        count = metadata["metadata"]["filters.stats"]["statistic"]["count"]
        return count
    except Exception as e:
        print(f"Failed to count points for {workunit}: {e}")
        return 0
    

def main(shapefile_aoi, h3_tiler, h3_layer_resolution):
    pass


# ──────────────────────────────────────────────────────────────────────────────
# Main — run one pipeline per 3DEP workunit
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # ── prep AOI & H3 tiles (unchanged) ─────────────────────────────────────────
    script_dir = os.path.dirname(__file__)
    shapefile_path = os.path.join(script_dir, "scrape_bounds.shp")
    user_aoi_4326 = import_shapefile_to_shapely(shapefile_path)
    h3_indices    = polygon_to_h3(user_aoi_4326)
    hex_gdf       = h3_cell_production(h3_indices)              

    # ── 3DEP index once, not in every loop ─────────────────────────────────────
    features      = get_usgs_datasets(user_aoi_4326, layer_id=8) 
    ds_gdf        = gpd.GeoDataFrame.from_features(features["features"],
                                                crs="EPSG:4326")

    print(ds_gdf)
    # ── MAIN LOOP: one DEM per (hex, dataset) ───────────────────────────────────
    run_commands = []
    for _, hex_row in hex_gdf.iterrows():
        h3_id   = hex_row.h3_id
        hex_geom = hex_row.geometry                  
        # keep only datasets that intersect this hex
        intersects = ds_gdf[ds_gdf.intersects(hex_geom)]
        if intersects.empty:
            continue                                 
        
        for _, ds_row in intersects.iterrows():
            workunit = ds_row.workunit 
            out_root = Path(f"C:\\data_engineering\\pdal\\{workunit}")
            os.makedirs(out_root,exist_ok=True)
            # original property name
            try:

                pipeline_dict = make_DEM_pipeline(
                    extent_epsg4326     = hex_geom,
                    usgs_3dep_dataset_name = [workunit],
                    pc_resolution       = 0.5,
                    dem_resolution      = 1.0,
                    filterNoise         = True,
                    reclassify          = False,
                    savePointCloud      = True,
                    outCRS              = 3857,
                    pc_outName          = out_root / f"{h3_id}_{workunit}_pc",
                    pc_outType          = "laz",
                    demType             = "dtm",
                    gridMethod          = "idw",
                    dem_outName         = out_root / f"{h3_id}_{workunit}_dem",
                    dem_outExt          = "tif",
                    driver              = "GTiff",
                )
                if pipeline_dict:
                    print(f"{h3_id}:{workunit}")
                    command_json = json.dumps(pipeline_dict)
                    run_commands.append(command_json)
                    pdal.Pipeline(json.dumps(pipeline_dict)).execute()
                else:
                    print(f"{workunit} not valid name for ept")
            except Exception as e:
                print(f"{e}")
