# generate basins


import os
import time
import glob

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import geopandas as gpd
import pandas as pd
import rioxarray as rxr

import richdem as rd

import random

import multiprocessing as mp

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'processed_data/')

DEM_resolution = 90 # EarthEnv DEM is 90m resolution
basin_threshold = int(1E6 / (90 * 90)) # min number of cells comprising a basin

processed_data_dir = os.path.join(DATA_DIR, 'processed_dem')
derived_basins_dir = os.path.join(DATA_DIR, 'derived_basins')

processed_regions = [e for e in os.listdir(derived_basins_dir) if not e.endswith('.geojson')]


def retrieve_polygon(fname):
    shp_num = fname.split('_')[-1].split('.')[0]
    return (shp_num, gpd.read_file(fname))


def retrieve_raster(fpath):
    rds = rxr.open_rasterio(fpath, masked=True, mask_and_scale=True)
    crs = rds.rio.crs
    affine = rds.rio.transform(recalc=False)
    return rds, crs, affine


def clip_raster_to_basin(basin_polygon, raster):
    
    crs = raster.rio.crs.to_epsg()
    if not crs:
        crs = raster.rio.crs.to_wkt()
    
    basin_polygon = basin_polygon.to_crs(crs)
    bounds = tuple(basin_polygon.bounds.values[0])

    # trimmed_box = box(*bounds).bounds
    try:
        subset_raster = raster.copy().rio.clip_box(*bounds)
        # bounds_area = (bounds[2] - bounds[0]) * (bounds[3] - bounds[1]) / 1E6
        # poly_area = basin_polygon.geometry.area.values[0] / 1E6
        clipped_raster = subset_raster.rio.clip(basin_polygon.geometry, basin_polygon.crs, all_touched=True)
        return clipped_raster, True
    except Exception as e:
        print(e)
        return None, False


def process_basin_elevation(clipped_raster):
    # evaluate masked raster data
    values = clipped_raster.data.flatten()
    mean_val = np.nanmean(values)
    median_val = np.nanmedian(values)
    min_val = np.nanmin(values)
    max_val = np.nanmax(values)
    return mean_val, median_val, min_val, max_val


def calculate_gravelius_and_perimeter(polygon):
    
    p = polygon.to_crs('EPSG:3005')
    perimeter = p.geometry.length.values[0]
    area = p.geometry.area.values[0] 
    if area == 0:
        return np.nan, perimeter
    else:
        perimeter_equivalent_circle = np.sqrt(4 * np.pi * area)
        gravelius = perimeter / perimeter_equivalent_circle
    
    return gravelius, perimeter



def calculate_slope_and_aspect(clipped_raster):   
    
    wkt = clipped_raster.rio.crs.to_wkt()
    affine = clipped_raster.rio.transform()
    
    rdem_clipped = rd.rdarray(
        clipped_raster.data[0], 
        no_data=clipped_raster.rio.nodata, 
        # projection=wkt, 
    )
    rdem_clipped.geotransform = affine.to_gdal()
    rdem_clipped.projection = wkt

    slope = rd.TerrainAttribute(rdem_clipped, attrib='slope_degrees')
    aspect = rd.TerrainAttribute(rdem_clipped, attrib='aspect')
    # print('slope, aspect')
    # print(slope, aspect)
    return np.nanmean(slope), np.nanmean(aspect)



def main():
    methods = ['RND', 'NBR', 'ACC']

    ppt_sample_size = 10
    take_sample = False

    for region in processed_regions:
        for method in methods:
            t_start = time.time()
            print(f'Processing {region}.')

            basin_path = os.path.join(DATA_DIR, f'derived_basins/{region}/{region}_derived_basins_{method}.geojson')
            basin_df = gpd.read_file(basin_path)

            if take_sample:
                sample_idx = random.choices(basin_df.index, k=ppt_sample_size)
                basin_df = basin_df.loc[sample_idx, :].copy().sort_index()
                print(f'   ...processing {len(basin_df)} basins.')

            dem_path = os.path.join(DATA_DIR, f'processed_dem/EENV_DEM/{region}_EENV_DEM_3005.tif')
            region_raster, _, _  = retrieve_raster(dem_path)

            output_data_folder = os.path.join(DATA_DIR, f'basin_properties/{region}/')
            if not os.path.exists(output_data_folder):
                os.makedirs(output_data_folder)
            
            output_fpath = os.path.join(output_data_folder, f'{region}_basin_properties_{method}.csv')
            
            for i, row in basin_df.iterrows():
                
                basin_polygon = basin_df[basin_df.index == i]

                clipped_raster, raster_loaded = clip_raster_to_basin(basin_polygon, region_raster)
                
                mean_el, median_el, min_el, max_el = process_basin_elevation(clipped_raster)
                basin_df.loc[i, 'Drainage_Area_km2'] = row['geometry'].area / 1E6
                basin_df.loc[i, 'mean_el'] = mean_el
                # median elevation, but label it the same as hysets column (Elevation_m)
                basin_df.loc[i, 'Elevation_m'] = median_el
                # basin_df.loc[i, 'min_el'] = min_el
                # basin_df.loc[i, 'max_el'] = max_el
                
                gravelius, perimeter = calculate_gravelius_and_perimeter(basin_polygon)
                
                basin_df.loc[i, 'Gravelius'] = gravelius
                basin_df.loc[i, 'Perimeter'] = perimeter

                slope, aspect = calculate_slope_and_aspect(clipped_raster)
                basin_df.loc[i, 'Slope_deg'] = slope
                basin_df.loc[i, 'Aspect_deg'] = aspect

            # drop the geometry column to save disk and write time
            basin_df = basin_df[[c for c in basin_df.columns if c != 'geometry']]
            basin_df.to_csv(output_fpath)

            t_end = time.time()
            unit_time = (t_end-t_start) / len(basin_df)
            print(f'   {t_end-t_start:.1f}s to delineate {len(basin_df)} basins.  ({unit_time:.2f}s per basin.)')
            print('')

        # print(sdafdasd)


if __name__ == '__main__':
    t0 = time.time()
    main()
    t1 = time.time()
    print('')
    print('###################################################')
    print('')
    print(f'Script completed in {t1-t0:.1}s.')
    print('__________________________________________________')