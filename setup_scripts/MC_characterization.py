# generate basins


import os
import time
import glob
import json

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import geopandas as gpd
import pandas as pd
import rioxarray as rxr

import random

import multiprocessing as mp

from shapely.geometry import Point

from whitebox.whitebox_tools import WhiteboxTools

wbt = WhiteboxTools()
wbt.verbose = False


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'processed_data/')
EXT_MEDIA = '/media/danbot/Samsung_T5/geospatial_data/basin_generator/'

# DATA_DIR = EXT_MEDIA

PPT_DIR = EXT_MEDIA

DEM_resolution = 90 # EarthEnv DEM is 90m resolution
basin_threshold = int(1E6 / (90 * 90)) # min number of cells comprising a basin

dem_folder = os.path.join(EXT_MEDIA, 'processed_dem/EENV_DEM')
region_files = os.listdir(dem_folder)
region_codes = sorted(list(set([e.split('_')[0] for e in region_files])))

basin_tracker_filename = 'delineated_basin_tracker.json'
basin_tracker_fpath = os.path.join(DATA_DIR, basin_tracker_filename)

basin_tracker = {}
if os.path.exists(basin_tracker_fpath):
    with open(basin_tracker_fpath) as f:
        basin_tracker = json.load(f)
        # print(basin_tracker)

region_codes = [
    # '08P','08O', '07G', 
    # '07U', '07O', 
    '07U',
    # '08G', '08H', '08E', '08A', 
    # '08D', '09A', '08F', '08B', 
    # '08C', 'ERockies', '08N', 
    # 'Peace', 
    # 'Fraser', 
    # 'Liard',
    ]

def retrieve_raster(region, raster_type):
    filename = f'{region}_EENV_DEM_3005_{raster_type}.tif'
    fpath = os.path.join(dem_folder, filename)
    raster = rxr.open_rasterio(fpath, mask_and_scale=True)
    crs = raster.rio.crs
    affine = raster.rio.transform(recalc=False)
    return raster, crs, affine


def get_region_area(region):
    polygon_path = os.path.join(DATA_DIR, 'merged_basin_groups/region_polygons/')
    poly_files = os.listdir(polygon_path)
    file = [e for e in poly_files if e.startswith(region)]
    if len(file) != 1:
        raise Exception; 'Region shape file not found.'
    fpath = os.path.join(polygon_path, file[0])
    gdf = gpd.read_file(fpath)
    gdf = gdf.to_crs(3005)
    return gdf['geometry'].area.values[0] / 1E6


def retrieve_polygon(fname):
    shp_num = fname.split('_')[-1].split('.')[0]
    return (shp_num, gpd.read_file(fname))


def check_polygon_df(region):
    poly_path = os.path.join(DATA_DIR, f'derived_basins/{region}_derived_basin_sample.geojson')
    if os.path.exists(poly_path):
        return True
    else:
        return False


def create_pour_point_gdf(stream, acc, pts, crs, n_chunks=2):
    """Break apart the list of stream pixels to avoid memory 
    allocation issue when indexing large rasters.

    Args:
        stream (_type_): _description_
        confluences (_type_): _description_
        n_chunks (int, optional): _description_. Defaults to 2.

    Returns:
        _type_: _description_
    """
    
    n_chunks = int(10 * np.log(len(pts)))

    # print(f'{len(pts)} points, breaking into {n_chunks} chunks.')

    conf_chunks = np.array_split(np.asarray(pts), indices_or_sections=n_chunks)

    point_array = []
    acc_array, id_array = [], []
    for chunk in conf_chunks:

        xis = [int(c[0]) for c in chunk]
        yis = [int(c[1]) for c in chunk]
        acc_array += [acc.data[0][c[0], c[1]] for c in chunk]

        ppts = stream[0, xis, yis]
        coords = tuple(map(tuple, zip(ppts.coords['x'].values, ppts.coords['y'].values)))
        point_array += [Point(p) for p in coords]

    df = pd.DataFrame()
    df['num_acc_cells'] = acc_array
    df['pt_id'] = list(range(len(df)))

    gdf = gpd.GeoDataFrame(df, geometry=point_array, crs=crs)
    # print(f'Created dataframe with {len(gdf)} pour points')
    return gdf


def retrieve_stream_px(region):
    # retrieve the stream raster 
    stream_data, crs, affine = retrieve_raster(region, 'pruned_stream')
    acc_data, crs, affine = retrieve_raster(region, 'accum')

    S = stream_data.data[0]      
    # get a list of stream pixel indices  
    stream_px = np.argwhere( S == 1 )

    return stream_px, stream_data, acc_data, crs


# def create_RAND_sample(stream_px, stream_data, acc_data, sample_size, crs):
#     # select a random sample of stream pixels
#     sample_ppts = random.choices(stream_px, k=sample_size)
#     n_chunks = int(np.log(len(stream_px)))
#     # create a pour point dataframe
#     ppts_gdf = create_pour_point_gdf(stream_data, acc_data, sample_ppts, crs, n_chunks=n_chunks)
#     return ppts_gdf


def create_simulation_ppt_file(region, ppt_gdf, basin_tracker, temp_filepath, filesize, sample_size):

    temp_gdf = ppt_gdf.copy()

    if region not in basin_tracker.keys():
        basin_tracker[region] = []
    else:
        pass
        # filter out any pour point indices for which a basin 
        # has already been delineated
        # existing_ppt_indices = basin_tracker[region]
        # temp_gdf = temp_gdf[~temp_gdf['raster_idx'].isin(existing_ppt_indices)].copy()
        # n_existing = len(temp_gdf) - len(ppt_gdf)

        # Figure out how to map pour points (raster index)
        # to the basin output from the unnest_basins() function
        # so that repeated basin delineations can be avoided.
        # we can know it was tracked, but sharing results between
        # simulations will take some more thought...

        # print('#############################')
        # print('')
        # print(f'Skipping {n_existing} pour points with an existing basin polygon.')
        # print('')
        # print('#############################')
        # print('')

        # existing_ppt_indices += list(temp_gdf['raster_idx'].values)

        # basin_tracker[region] = list(set(existing_ppt_indices))
        # with open(basin_tracker_fpath, 'w') as json_file:
        #     # print('Updating basin tracker file.')
        #     json.dump(basin_tracker, json_file)

    # divide the dataframe into chunks for batch processing
    # save the pour point dataframe to temporary files
    # and limit temporary raster files to 10GB / batch
    batch_limit = 1E5
    n_batches = int(filesize * sample_size / batch_limit) + 1
    print(f' Running {n_batches} batch(es) on {filesize:.1f}MB raster.')
    batch_paths = []
    n = 0
    temp_gdf['FID'] = temp_gdf['pt_id'].copy()
    if sample_size * filesize < batch_limit:
        temp_fpath = temp_filepath.split('.')[0] + f'_{n:04}.shp'
        temp_gdf.to_file(temp_fpath)
        batch_paths.append(temp_fpath)
    else:                    
        batches = np.array_split(temp_gdf.index, indices_or_sections=n_batches)
        for batch in batches:
            batch_gdf = temp_gdf[temp_gdf.index.isin(batch)].copy()            
            temp_fpath = temp_filepath.split('.')[0] + f'_{n:04}.shp'
            batch_gdf.to_file(temp_fpath)
            n += 1
            batch_paths.append(temp_fpath)    
    return batch_paths


def batch_basin_delineation(fdir_path, ppt_path, temp_raster_path):

    wbt.unnest_basins(
        fdir_path, 
        ppt_path, 
        temp_raster_path,
        esri_pntr=False, 
        # callback=default_callback
    )


def raster_to_vector_basins(temp_raster_folder, polygon_folder):
    batch_rasters = sorted(os.listdir(temp_raster_folder))
    n = 0
    for f in batch_rasters:
        fpath = os.path.join(temp_raster_folder, f)
        output_polygon_path = os.path.join(polygon_folder, f'temp_polygons_{n:04}.shp')
        # this function creates rasters of ordered sets of 
        # non-overlapping basins
        wbt.raster_to_vector_polygons(
            fpath,
            output_polygon_path,
        )
        n += 1


def save_gdf(output_polygon_folder, region, method, n_sim, all_polygons):
    crs = all_polygons[0].crs
    print(f'The crs is {crs}')
    print('')
    print(all_polygons[0].head())
    all_gdfs = gpd.GeoDataFrame(pd.concat(all_polygons), crs=crs)
    print(len(all_gdfs))
    print('')
    print(asdfasdfsd)
    gdf_fpath = os.path.join(output_polygon_folder, f'{region}_basins_{method}_{n_sim:04}.geojson')
    all_gdfs.to_file(gdf_fpath, driver='GeoJSON')


def reformat_basin_polygons_and_cleanup(temp_raster_path, temp_polygon_path, crs, min_basin_area):

    raster_files = os.listdir(temp_raster_path)
    for f in raster_files:
        # print(f'removing {f}')
        os.remove(os.path.join(temp_raster_path, f))

    gdfs = []
    polygon_files = [e for e in os.listdir(temp_polygon_path) if e.endswith('.shp')]

    for f in polygon_files:
        
        fp = os.path.join(temp_polygon_path, f)
        gdf = gpd.read_file(fp)
        # gdf = gdf.explode(index_parts=False)
        gdf['area'] = gdf.geometry.area / 1E6
        gdf = gdf[gdf['area'] > min_basin_area]
        
        gdfs.append(gdf)
    
    temp_polygon_files = os.listdir(temp_polygon_path) 

    for f in temp_polygon_files:
        os.remove(os.path.join(temp_polygon_path, f))

    all_gdfs = gpd.GeoDataFrame(pd.concat(gdfs), crs=crs)
    all_gdfs.reset_index(inplace=True)

    return all_gdfs
    

def main():

    methods = ['RAND', 'CONF', 'GRAB', 'GRAP']
    # methods = ['CONF', 'GRAD']
    
    n_simulations = 50
    # set the sample size based on area 
    # density of stations: 1 station per 100 km^2
    simulated_station_density = 100
    revision = 'RD'
    min_basin_area = 2.0 # km^2

    for region in region_codes:
        # if region == 'Fraser':
        #     n_simulations = 8
        # else:
        #     methods = ['RAND', 'CONF', 'GRAD']
        #     n_simulations = 50

        # define the sample size
        region_area = get_region_area(region)
        sample_size = int(region_area / simulated_station_density)

        fdir_file = f'{region}_EENV_DEM_3005_temp_fdir.tif'
        fdir_path = os.path.join(dem_folder, fdir_file)

        filesize = (os.path.getsize(fdir_path) >> 20 )

        temp_basin_raster_folder = os.path.join(DATA_DIR, f'derived_basins/{region}/basin_rasters/')
        temp_polygon_folder = os.path.join(DATA_DIR, f'derived_basins/temp/')
        
        ppt_folder = os.path.join(PPT_DIR, f'pour_points/{region}/')
        if not os.path.exists(ppt_folder):
            os.makedirs(ppt_folder)
        temp_ppt_filepath = os.path.join(ppt_folder, 'temp_ppt_sample.shp')

        if 'media' in EXT_MEDIA:
            output_polygon_folder = os.path.join(EXT_MEDIA, f'{revision}/derived_basins/{region}')
        else:
            output_polygon_folder = os.path.join(EXT_MEDIA, f'derived_basins/{region}')

        if not os.path.exists(output_polygon_folder):
            os.makedirs(output_polygon_folder)

        for f in [temp_basin_raster_folder, temp_polygon_folder]:
            if not os.path.exists(f):
                os.makedirs(f)

        for method in methods:
            t_start = time.time()

            print(f'Processing {region}: {n_simulations} simulations of {sample_size} basin samples using {method} method.')

            existing_sims = sorted([e for e in os.listdir(output_polygon_folder) if method in e])
            
            # if there are existing simulations, continue numbering
            last_sim = len(existing_sims) + 1

            ppt_fname = f'{region}_pour_pts_{method}.geojson'
            ta = time.time()
            all_ppts = gpd.read_file(os.path.join(ppt_folder, ppt_fname))
            tb = time.time()
            print(f'   ...time to read ppt file: {tb-ta:.1f}s')

            for n_sim in range(last_sim, last_sim + n_simulations, 1):

                t0 = time.time()
                            
                # select a random sample from the dataframe
                ppt_gdf = all_ppts.sample(n=sample_size)
                t_batch1 = time.time()
                print(f'Time to read sample: {t_batch1-t0:.1f}')

                simulation_fpaths = create_simulation_ppt_file(region, ppt_gdf, basin_tracker, temp_ppt_filepath, filesize, sample_size)

                print(simulation_fpaths)
                print('')
                print(simulation_fpaths[0])
                print(asfsd)

                all_polygons = []
                for ppt_path in simulation_fpaths:

                    t_ppt = time.time()
                    
                    temp_fname = f'{region}_{method}_{n_sim:04}.tif'
                    temp_basin_raster_path = os.path.join(temp_basin_raster_folder, temp_fname)

                    tdel = time.time()
                    
                    batch_basin_delineation(fdir_path, ppt_path, temp_basin_raster_path)

                    trp = time.time()
                    raster_to_vector_basins(temp_basin_raster_folder, temp_polygon_folder)

                    pclean = time.time()

                    polygons = reformat_basin_polygons_and_cleanup(temp_basin_raster_folder, temp_polygon_folder, output_polygon_folder, ppt_gdf.crs, n_sim, min_basin_area)

                    all_polygons.append(polygons)
                    
                    t_batch2 = time.time()
                    if n_sim % 10 == 0:    
                        t_clean = t_batch2 - pclean
                        trtv = pclean - trp
                        bbdt = trp - tdel
                        print(f'    ...Batch basin delineation: {bbdt:.3f}')
                        print(f'    ...raster to vector basins: {trtv:.3f}')
                        print(f'    ...reformat to geojson: {t_clean:.3f}')

                    # run batch feature extraction
                    # store the dropped indexes and add 
                    # their associated properties 
                    # (can't, they aren't ordered)

                save_gdf(output_polygon_folder, region, method, n_sim, all_polygons)
                
                if n_sim % 10 == 0:
                    t_end = time.time()
                    unit_time = (t_end - t_start) / (n_sim - last_sim)
                    print(f'   ...simulation {n_sim} processed in {t_end-t_start:.1f}s ({unit_time:.1f})s/sim..')
            

if __name__ == '__main__':
    t0 = time.time()
    main()
    t1 = time.time()
    print('')
    print('###################################################')
    print('')
    print(f'Script completed in {t1-t0:.2f}s.')
    print('__________________________________________________')