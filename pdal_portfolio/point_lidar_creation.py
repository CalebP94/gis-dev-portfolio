import pdal
import json
import numpy as np
import pandas as pd
import os
from shapely.geometry import Polygon
import pyproj
import boto3
import requests
import geopandas as gpd
import utils


def get_extent_from_ept_s3(bucket: str='usgs-lidar-public', key: str='bounds'):
    s3 = boto3.client("s3")
    response = s3.get_object(Bucket=bucket, Key=key)
    metadata = json.loads(response['Body'].read().decode('utf-8'))

    bounds = metadata["bounds"]  # [minx, miny, minz, maxx, maxy, maxz]
    minx, miny = bounds[0], bounds[1]
    maxx, maxy = bounds[3], bounds[4]

    # Transform from EPSG:3857 to EPSG:4326
    transformer = pyproj.Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
    sw_lon, sw_lat = transformer.transform(minx, miny)
    ne_lon, ne_lat = transformer.transform(maxx, maxy)

    return sw_lon, sw_lat, ne_lon, ne_lat



def import_shapefile_to_shapely(path):
    """
    Load shapefile and return the geometry as a shapely polygon (assumed EPSG:3857)
    """
    gdf = gpd.read_file(path)
    return gdf.loc['geometry']


def build_pdal_pipeline(extent_epsg3857, usgs_3dep_dataset_names, pc_resolution, filterNoise = False,
                        reclassify = False, savePointCloud = True, outCRS = 3857, pc_outName = 'filter_test', 
                        pc_outType = 'laz'):

    """
    Build pdal pipeline for requesting, processing, and saving point cloud data. Each processing step is a 'stage' 
    in the final pdal pipeline. Each stages is appended to the 'pointcloud_pipeline' object to produce the final pipeline.
    
    Parameters:
    extent_epsg3857 (shapely polygon): Polygon for user-defined AOI in Web Mercator projection (EPS:3857)Polygon is generated 
                            either through the 'handle_draw' methor or by inputing their own shapefile.
    usgs_3dep_dataset_names (str): List of name of the 3DEP dataset(s) that the data will be obtained. This parameter is set 
                                determined through intersecttino of the 3DEP and AOI polys.
    pc_resolution (float): The desired resolution of the pointcloud based on the following definition:
        
                        Source: https://pdal.io/stages/readers.ept.html#readers-ept
                            A point resolution limit to select, expressed as a grid cell edge length. 
                            Units correspond to resource coordinate system units. For example, 
                            for a coordinate system expressed in meters, a resolution value of 0.1 
                            will select points up to a ground resolution of 100 points per square meter.
                            The resulting resolution may not be exactly this value: the minimum possible 
                            resolution that is at least as precise as the requested resolution will be selected. 
                            Therefore the result may be a bit more precise than requested.
                            
    filterNoise (bool): Option to remove points from USGS Class 7 (Low Noise) and Class 18 (High Noise).
    reclassify (bool): Option to remove USGS classes and run SMRF to classify ground points only. Default == False.
    savePointCloud (bool): Option to save (or not) the point cloud data. If savePointCloud == False, 
           the pc_outName and pc_outType parameters are not used and can be any value.
    outCRS (int): Output coordinate reference systemt (CRS), specified by ESPG code (e.g., 3857 - Web Mercator)
    pc_outName (str): Desired name of file on user's local file system. If savePointcloud = False, 
                  pc_outName can be in value.
    pc_outType (str):  Desired file extension. Input must be either 'las' or 'laz'. If savePointcloud = False, 
                  pc_outName can be in value. If a different file type is requested,the user will get error.
    
    Returns:
        pointcloud_pipeline (dict): Dictionary of processing stages in sequential order that define PDAL pipeline.

    Raises: 
        Exception: If user passes in argument that is not 'las' or 'laz'.
    """
    
    #this is the basic pipeline which only accesses the 3DEP data
    readers = []

    url = "https://s3-us-west-2.amazonaws.com/usgs-lidar-public/{}/ept.json".format(usgs_3dep_dataset_names)
    reader = {
        "type": "readers.ept",
        "filename": str(url),
        "polygon": str(extent_epsg3857),
        "requests": 3,
        "resolution": pc_resolution
    }
    readers.append(reader)
        
    pointcloud_pipeline = {
            "pipeline":
                readers
    }
    
    if filterNoise == True:
        
        # Filter stage for class 7
        filter_stage_class7 = {
            "type": "filters.range",
            "limits": "Classification![7:7]"
        }

        # Filter stage for class 18
        filter_stage_class18 = {
            "type": "filters.range",
            "limits": "Classification![18:18]"
        }

        # Append both filter stages to the pipeline separately
        pointcloud_pipeline['pipeline'].append(filter_stage_class7)
        pointcloud_pipeline['pipeline'].append(filter_stage_class18)
    
    if reclassify == True:
        
        remove_classes_stage = {
            "type":"filters.assign",
            "value":"Classification = 0"
        }
        
        classify_ground_stage = {
            "type":"filters.smrf"
        }
        
        reclass_stage = {
            "type":"filters.range",
            "limits":"Classification[2:2]"
        }

        pointcloud_pipeline['pipeline'].append(remove_classes_stage)
        pointcloud_pipeline['pipeline'].append(classify_ground_stage)
        pointcloud_pipeline['pipeline'].append(reclass_stage)
        
    reprojection_stage = {
        "type":"filters.reprojection",
        "out_srs":"EPSG:{}".format(outCRS)
    }
    
    pointcloud_pipeline['pipeline'].append(reprojection_stage)
    
    if savePointCloud == True:
        
        if pc_outType == 'las':
            savePC_stage = {
                "type": "writers.las",
                "filename": str(pc_outName)+'.'+ str(pc_outType),
            }
        elif pc_outType == 'laz':    
            savePC_stage = {
                "type": "writers.las",
                "compression": "laszip",
                "filename": str(pc_outName)+'.'+ str(pc_outType),
            }
        else:
            raise Exception("pc_outType must be 'las' or 'laz'.")

        pointcloud_pipeline['pipeline'].append(savePC_stage)
        
    return pointcloud_pipeline

def make_DEM_pipeline(extent_epsg3857, usgs_3dep_dataset_name, pc_resolution, dem_resolution,
                      filterNoise = True, reclassify = False, savePointCloud = False, outCRS = 3857,
                      pc_outName = 'filter_test', pc_outType = 'laz', demType = 'dtm', gridMethod = 'idw', 
                      dem_outName = 'dem_test', dem_outExt = 'tif', driver = "GTiff"):
    
    """
    Build pdal pipeline for creating a digital elevation model (DEM) product from the requested point cloud data. The 
    user must specify whether a digital terrain (bare earth) model (DTM) or digital surface model (DSM) will be created, 
    the output DTM/DSM resolution, and the gridding method desired. 

    The `build_pdal_pipeline() method is used to request the data from the Amazon Web Services ept bucket, and the 
    user may define any processing steps (filtering, reclassifying, reprojecting). The user must also specify whether 
    the point cloud should be saved or not. Saving the point cloud is not necessary for the generation of the DEM. 

    Parameters:
        extent_epsg3857 (shapely polygon): User-defined AOI in Web Mercator projection (EPS:3857). Polygon is generated 
                                           either through the 'handle_draw' methor or by inputing their own shapefile.
                                           This parameter is set automatically when the user-defined AOI is chosen.
        usgs_3dep_dataset_names (list): List of name of the 3DEP dataset(s) that the data will be obtained. This parameter is set 
                                        determined through intersecttino of the 3DEP and AOI polys.
        pc_resolution (float): The desired resolution of the pointcloud based on the following definition:

                        Source: https://pdal.io/stages/readers.ept.html#readers-ept
                            A point resolution limit to select, expressed as a grid cell edge length. 
                            Units correspond to resource coordinate system units. For example, 
                            for a coordinate system expressed in meters, a resolution value of 0.1 
                            will select points up to a ground resolution of 100 points per square meter.
                            The resulting resolution may not be exactly this value: the minimum possible 
                            resolution that is at least as precise as the requested resolution will be selected. 
                            Therefore the result may be a bit more precise than requested.

        pc_outName (str): Desired name of file on user's local file system. If savePointcloud = False, 
                          pc_outName can be in value.
        pc_outType (str): Desired file extension. Input must be either 'las' or 'laz'. If savePointcloud = False, 
                          pc_outName can be in value. If a different file type is requested,the user will get error.
    
        dem_resolution (float): Desired grid size (in meters) for output raster DEM 
        filterNoise (bool): Option to remove points from USGS Class 7 (Low Noise) and Class 18 (High Noise).
        reclassify (bool): Option to remove USGS classes and run SMRF to classify ground points only. Default == False.
        savePointCloud (bool): Option to save (or not) the point cloud data. If savePointCloud == False, the pc_outName 
                               and pc_outType parameters are not used and can be any value.

        outCRS (int): Output coordinate reference systemt (CRS), specified by ESPG code (e.g., 3857 - Web Mercator)
        pc_outName (str): Desired name of file on user's local file system. If savePointcloud = False, 
                          pc_outName can be in value.
        pc_outType (str): Desired file extension. Input must be either 'las' or 'laz'. If a different file type is requested,
                    the user will get error stating "Extension must be 'las' or 'laz'". If savePointcloud = False, 
                    pc_outName can be in value.
        demType (str): Type of DEM produced. Input must 'dtm' (digital terrain model) or 'dsm' (digital surface model).
        gridMethod (str): Method used. Options are 'min', 'mean', 'max', 'idw'.
        dem_outName (str): Desired name of DEM file on user's local file system.
        dem_outExt (str): DEM file extension. Default is TIF.
        driver (str): File format. Default is GTIFF
    
    Returns:
        dem_pipeline (dict): Dictionary of processing stages in sequential order that define PDAL pipeline.
    Raises: 
        Exception: If user passes in argument that is not 'las' or 'laz'.
        Exception: If user passes in argument that is not 'dtm' or 'dsm'

    """

    dem_pipeline = build_pdal_pipeline(extent_epsg3857, usgs_3dep_dataset_name, pc_resolution,
                                              filterNoise, reclassify, savePointCloud, outCRS, pc_outName, pc_outType)
    
    if demType == 'dsm':
        dem_stage = {
                "type":"writers.gdal",
                "filename":str(dem_outName)+ '.' + str(dem_outExt),
                "gdaldriver":driver,
                "nodata":-9999,
                "output_type":gridMethod,
                "resolution":float(dem_resolution),
                "gdalopts":"COMPRESS=LZW,TILED=YES,blockxsize=256,blockysize=256,COPY_SRC_OVERVIEWS=YES"
        }
    
    elif demType == 'dtm':
        groundfilter_stage = {
                "type":"filters.range",
                "limits":"Classification[2:2]"
        }

        dem_pipeline['pipeline'].append(groundfilter_stage)

        dem_stage = {
                "type":"writers.gdal",
                "filename":str(dem_outName)+ '.' + str(dem_outExt),
                "gdaldriver":driver,
                "nodata":-9999,
                "output_type":gridMethod,
                "resolution":float(dem_resolution),
                "gdalopts":"COMPRESS=LZW,TILED=YES,blockxsize=256,blockysize=256,COPY_SRC_OVERVIEWS=YES"
        }
    
    else:
        raise Exception("demType must be 'dsm' or 'dtm'.")
        
    dem_pipeline['pipeline'].append(dem_stage)
    
    return dem_pipeline


def get_usgs_datasets(user_aoi_3857):
    """
    Queries USGS 3DEP Elevation Index with AOI bounding box and returns dataset names
    """
    # Step 1: Get bounds of AOI in EPSG:3857
    bounds = user_aoi_3857.bounds  # (minx, miny, maxx, maxy)
    geometry_str = f"{bounds[0]},{bounds[1]},{bounds[2]},{bounds[3]}"

    # Step 2: Set up query parameters
    url = "https://index.nationalmap.gov/arcgis/rest/services/3DEPElevationIndex/MapServer/11/query"
    params = {
        "where": "1=1",
        "geometry": geometry_str,
        "geometryType": "esriGeometryEnvelope",
        "inSR": "3857",
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "*",
        "returnGeometry": "true",
        "f": "geojson"
    }

    # Step 3: Send request
    response = requests.get(url, params=params)
    return response.json()




if __name__ == '__main__':

    shapefile_path = r"C:\\data_engineering\\pdal\\scrape_bounds.shp"
    user_aoi_3857 = utils.import_shapefile_to_shapely(shapefile_path)
    data = utils.get_usgs_datasets(user_aoi_3857, layer_id=8)
    for feature in data['features']:
        workunit = f"{feature['properties']['workunit']}"
        pipeline = make_DEM_pipeline(
            extent_epsg3857=user_aoi_3857,
            usgs_3dep_dataset_name=f"{workunit}",  # Replace with your dataset name
            pc_resolution=0.5,
            dem_resolution=1.0,
            filterNoise=True,
            reclassify=False,
            savePointCloud=True,
            outCRS=3857,
            pc_outName=f"C:\\data_engineering\\pdal\\{workunit}_pc",
            pc_outType="laz",
            demType="dtm",
            gridMethod="idw",
            dem_outName=f"C:\\data_engineering\\pdal\\{workunit}_dem",
            dem_outExt="tif",
            driver="GTiff"
        )


