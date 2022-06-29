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

import random

import multiprocessing as mp

from whitebox.whitebox_tools import WhiteboxTools

wbt = WhiteboxTools()
wbt.verbose = False
# wbt.rust_backtrace = True

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'processed_data/')

DEM_resolution = 90 # EarthEnv DEM is 90m resolution
basin_threshold = int(1E6 / (90 * 90)) # min number of cells comprising a basin

processed_data_dir = os.path.join(DATA_DIR, 'processed_dem')

region_files = os.listdir(os.path.join(processed_data_dir, 'EENV_DEM'))
region_codes = sorted(list(set([e.split('_')[0] for e in region_files])))

region_codes = [
    #.07-.08 .26 .34-.35 .42-.44
    # '08P', '08O', '07G', '07U',
    #                      1.33
    #'08G', '08H', '08E', '08A',
    # '08A',
    '08D', '09A', '08F', '08B', '08C',
    'ERockies', '08N', 'Peace', 
    'Fraser', 'Liard'
    ]

def delineate_basin(data):

    try:

        wbt.watershed(
            data['fdir_path'], 
            data['ppt_path'], 
            data['basin_raster_path'], 
            esri_pntr=False,
        )

        wbt.raster_to_vector_polygons(
            data['basin_raster_path'], 
            data['basin_polygon_path'], 
        )

    except Exception as ex:
        print(ex)
    
    try:
        os.remove(data['basin_raster_path'])
        return None
    except Exception as ex:
        return data


def unnest_basins(region, fdir_path, basin_output_path, method):
    pour_pt_folder = os.path.join(DATA_DIR, f'pour_points/{region}/{method}')
    ppt_shp_files = [e for e in os.listdir(pour_pt_folder) if e.endswith('.shp')]
    assert len(ppt_shp_files) == 1
    ppt_path = os.path.join(pour_pt_folder, ppt_shp_files[0])

    basins_raster_path = os.path.join(basin_output_path, f'{region}_basins.tif')

    wbt.unnest_basins(
        fdir_path, 
        ppt_path, 
        basins_raster_path,
        esri_pntr=False, 
        # callback=default_callback
    )

def retrieve_polygon(fname):
    shp_num = fname.split('_')[-1].split('.')[0]
    return (shp_num, gpd.read_file(fname))


def check_polygon_df(region):
    poly_path = os.path.join(DATA_DIR, f'derived_basins/{region}_derived_basin_sample.geojson')
    if os.path.exists(poly_path):
        return True
    else:
        return False


def clean_up_basin_polygons(output_polygon_folder):
    files = os.listdir(output_polygon_folder)
    for f in files:
        os.remove(os.path.join(output_polygon_folder, f))

    os.rmdir(output_polygon_folder)
    

def convert_results_to_geojson(region, polygon_folder, method):
    
    all_files = [f for f in os.listdir(polygon_folder) if f.endswith('.shp')]
    all_paths = [os.path.join(polygon_folder, f) for f in all_files]
    print(f'   ...converting {len(all_paths)} basins to geojson.')

    with mp.Pool() as p:
        results = p.map(retrieve_polygon, all_paths)

    ids = [e[0] for e in results]
    all_gdfs = [e[1] for e in results]
    
    merged = pd.concat(all_gdfs)

    merged['basin_id'] = ids

    gdf = gpd.GeoDataFrame(merged[['basin_id']], geometry=merged['geometry'].values, crs='EPSG:3005')
    n_basins = len(gdf)
    fname = f'derived_basins/{region}/{region}_derived_basins_{method}.geojson'
    gdf.to_file(os.path.join(DATA_DIR, fname), driver='GeoJSON')

    print(f'   ...{len(gdf)} basins merged to geodataframe.')
    return len(gdf)


def create_input_array(region, fdir_path, polygon_folder, basin_folder, ppt_folder):
    """create a list of dict items featuring the information needed 
    as input in to the watershed delineation function.

    Args:
        region (_type_): _description_
        fdir_path (_type_): _description_
        polygon_folder (_type_): _description_
        basin_folder (_type_): _description_
        ppt_folder (_type_): _description_
        make_gdf (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
    """
    all_ppt_files = os.listdir(ppt_folder)
    ppt_files = [e for e in all_ppt_files if e.endswith('.shp')]
    input_array = []
    for f in ppt_files:
        ppt_num = int(f.split('_')[-1].split('.')[0])
        basin_raster_filename = f'basin_{ppt_num:07d}.tif'
        basin_polygon_filename = f'basin_{ppt_num:07d}.shp'
        basin_raster_path = os.path.join(basin_folder, basin_raster_filename)
        basin_polygon_path = os.path.join(polygon_folder, basin_polygon_filename)

        ppt_path = os.path.join(ppt_folder, f)

        d = {
            'region': region,
            'fdir_path': fdir_path,
            'ppt_path': ppt_path,
            'ppt_num': ppt_num,
            'basin_raster_path': basin_raster_path,
            'basin_polygon_path': basin_polygon_path,
            }

        input_array.append(d)

    return input_array


# def load_ppt_file(region, method, ppt_folder):
#     ppt_file_prefix = f'{region}_ppt_{method}'
#     # ppt_path = os.path.join(ppt_folder, ppt_file_prefix)
#     all_ppt_files = os.listdir(ppt_folder)
#     shp_files = [e for e in all_ppt_files if e.endswith('.shp')]
#     ppt_file = [e for e in shp_files if e.startswith(ppt_file_prefix)]

#     assert len(ppt_file) == 1

#     ppt_path = os.path.join(ppt_folder, ppt_file[0])
#     ppt_gdf = gpd.read_file(ppt_path)

#     print(ppt_gdf)

#     return ppt_gdf


def main():

    methods = ['RAND', 'CONF', 'GRAD']

    use_unnest = False
    make_gdf = True

    # ppt_sample_size = 100
    for region in region_codes:
        for method in methods:
            t_start = time.time()
            print(f'Processing {region}.')
            ppt_folder = os.path.join(DATA_DIR, f'pour_points/{region}/{method}')

            # ppt_gdf = load_ppt_file(region, method, ppt_folder)

            fdir_file = f'EENV_DEM/{region}_EENV_DEM_3005_temp_fdir.tif'
            fdir_path = os.path.join(processed_data_dir, fdir_file)

            temp_basin_raster_folder = os.path.join(DATA_DIR, f'derived_basins/{region}/basin_rasters/')
            output_polygon_folder = os.path.join(DATA_DIR, f'derived_basins/{region}/basin_polygons/')

            for f in [temp_basin_raster_folder, output_polygon_folder]:
                if not os.path.exists(f):
                    os.makedirs(f)

            if use_unnest:
                unnest_basins(region, fdir_path, temp_basin_raster_folder)
            else:
                # write the individual pour points to separate shp files in parallel
                # first create a new folder to store the files if it doesn't exist
                ppt_folder = os.path.join(DATA_DIR, f'pour_points/{region}/{method}/')

                basin_input_array = create_input_array(region, fdir_path, output_polygon_folder, temp_basin_raster_folder, ppt_folder)

                t1 = time.time()

                with mp.Pool(5) as p:
                    p.map(delineate_basin, basin_input_array)

                n_basins_created = len(basin_input_array)
                t2 = time.time()
                print(f'{n_basins_created} pour point shapefiles written in {t2-t1:.1e}s')

            for f in [temp_basin_raster_folder]:
                os.rmdir(f)

            convert_results_to_geojson(region, output_polygon_folder, method)

            clean_up_basin_polygons(output_polygon_folder)

            t_end = time.time()
            print(t_end - t_start)
            unit_time = (t_end-t_start) / n_basins_created
            print(f'   {t_end-t_start:.1f}s to delineate {n_basins_created} basins.  ({unit_time:.2f}s per basin.)')
            print('')
            print('')
            print('')
            print('#################')
            print('')
        print(asdf)


if __name__ == '__main__':
    t0 = time.time()
    main()
    t1 = time.time()
    print('')
    print('###################################################')
    print('')
    print(f'Script completed in {t1-t0:.1}s.')
    print('__________________________________________________')