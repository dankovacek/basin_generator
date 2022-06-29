# generate basins

from concurrent.futures import process
import os
import sys

import random

import time

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as rxr

from shapely.geometry import Polygon, Point
from shapely.validation import make_valid

from multiprocessing import Pool

from whitebox.whitebox_tools import WhiteboxTools

wbt = WhiteboxTools()
wbt.verbose = False

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'processed_data/')

DEM_resolution = 90 # EarthEnv DEM is 90m resolution
basin_threshold = int(1E6 / (90 * 90)) # min number of cells comprising a basin

processed_data_dir = os.path.join(DATA_DIR, 'processed_dem')

region_files = os.listdir(os.path.join(processed_data_dir, 'EENV_DEM'))

region_codes = sorted(list(set([e.split('_')[0] for e in region_files])))


def retrieve_raster(region, raster_type):
    filename = f'{region}_EENV_DEM_3005_{raster_type}.tif'
    raster_path = os.path.join(processed_data_dir, f'EENV_DEM/{filename}')
    raster = rxr.open_rasterio(raster_path, mask_and_scale=True)
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


def random_stream_point_selection(S, ppt_sample_size):
    """
    Randomly select stream pixels.
    """
    acc, crs, affine = retrieve_raster(region, 'accum')
    A = acc.data[0]
    stream_px = np.argwhere( S == 1 )
    stream_sample = random.choices(stream_px, k=ppt_sample_size)
    ppts = []
    c = 0
    for i, j in stream_sample:
        cell_acc = A[i, j]
        ppts.append((i, j, cell_acc, c))
        c += 1

    pct_tracked = len(ppts) / len(stream_px) * 100
    print(f'Tracked {c} randomly selected points.')
    print(f'{len(ppts):.1e}/{len(stream_px):.1e} stream pixels are points of interest ({pct_tracked:.0f}%).')
        
    return ppts, pct_tracked


def find_stream_confluences_by_fdir(S, ppt_sample_size):
    """Assume confluences are stream pixels connected adjacent or diagonally with more than two other stream pixels.

    Args:
        stream (2d matrix): Binary matrix of stream pixels (1), non-stream pixels (0)

    Returns:
        _type_: list of pixel indices corresponding to confluences
    """
    # retrieve the flow direction raster
    fdir, crs, affine = retrieve_raster(region, 'temp_fdir')
    acc, crs, affine = retrieve_raster(region, 'accum')

    A = acc.data[0]

    F = fdir.data[0]

    # the matrix is very sparse, retrieve just the indices
    # corresponding to (all) stream pixels
    stream_px = np.argwhere( S == 1 )

    confluences = []
    c = 0
    for (i, j) in stream_px:
        
        # create a boolean matrix for stream cells
        S_W = np.where((S[i-1:i+2, j-1:j+2] == 1), 1, 0)
        # retrieve D8 pointer from a 3x3 window 
        # around the target pixel
        F_W = F[i-1:i+2, j-1:j+2].round(2) 
        # mask the flow direction by element-wise multiplication
        # using the boolean stream matrix
        W = np.multiply(F_W, S_W)

        # get the accumulation window to find the drainage area
        # of each inflow cell (as well as the target cell)
        A_W = np.multiply(A[i-1:i+2, j-1:j+2], S_W)        

        inflow_directions = np.array([
            [4, 8, 16],
            [2, np.nan, 32],
            [1, 128, 64]
        ])

        # find flow direction cells
        # flowing towards the center cell.
        inflow_cells = np.argwhere(W == inflow_directions)
        num_matches = np.sum(W == inflow_directions)
    
        if num_matches > 1:
            # add the target cell
            cell_acc = A_W[1, 1]
            confluences.append((i, j, cell_acc, c))
            # add the inflow cell indices
            inflow_cells = np.argwhere((W == inflow_directions))
            for ci, cj in inflow_cells:
                ix = ci + i - 1
                jx = cj + j - 1
                cell_acc = A_W[ci, cj]  
                confluences.append((ix, jx, cell_acc, c))
            c += 1

    # if len(confluences) > ppt_sample_size:
    #     return random.choices(confluences, k=ppt_sample_size)

    pct_found = len(confluences) / len(stream_px) * 100
    print(f'{c} confluences found.')
    print(f'{len(confluences):.1e}/{len(stream_px):.1e} stream pixels are points of interest ({pct_found:.0f}%).')

    return confluences, pct_found


def find_stream_confluences_by_acc_gradient(region, S, threshold, ppt_sample_size):
    """_summary_

    Args:
        region (_type_): _description_
        idxs (_type_): _description_
        threshold (_type_): _description_

    Returns:
        _type_: _description_
    """
    acc, crs, affine = retrieve_raster(region, 'accum')

    stream_px = np.argwhere( S == 1 )

    px_of_interest = []
    c = 0
    for i, j in stream_px:
        # get the accumulation value of the target cell
        cell_acc = acc.data[0][i, j]

        # create a boolean matrix for stream cells
        S_W = np.where((S[i-1:i+2, j-1:j+2] == 1), True, np.nan)
        
        # retrieve flow accumulation from a 3x3 window 
        # around the target pixel
        A_W = acc.data[0][i-1:i+2, j-1:j+2].round(2) 

        # mask the flow direction by element-wise multiplication
        # using the boolean stream matrix
        A = np.multiply(A_W, S_W)
        # print(f'threshold: {threshold}')
        # print(A)
        # print(asdfs)

        # num_stream_cells = np.count_nonzero(~np.isnan(W))

        # calculate the acc gradient w.r.t. the target cell
        acc_grad = cell_acc - A
        jump_px = np.argwhere( acc_grad > threshold )

        num_jumps = len(jump_px)

        if (num_jumps > 0):
            px_of_interest.append((i, j, cell_acc, c))
            c += 1

    # if len(px_of_interest) > ppt_sample_size:
    #     return random.choices(px_of_interest, k=ppt_sample_size)

    pct_found = len(px_of_interest) / len(stream_px) * 100
    print(f'{len(px_of_interest):.1e}/{len(stream_px):.1e} stream pixels are points of interest ({pct_found:.0f}%).')

    return px_of_interest, pct_found


def create_pour_point_gdf(stream, confluences, crs, n_chunks=2):
    """Break apart the list of stream pixels to avoid memory 
    allocation issue when indexing large rasters.

    Args:
        stream (_type_): _description_
        confluences (_type_): _description_
        n_chunks (int, optional): _description_. Defaults to 2.

    Returns:
        _type_: _description_
    """
    
    n_chunks = int(10 * np.log(len(confluences)))

    conf_chunks = np.array_split(np.asarray(confluences), indices_or_sections=n_chunks)

    point_array = []
    acc_array, id_array = [], []
    idx_array = []
    for chunk in conf_chunks:

        xis = [int(c[0]) for c in chunk]
        yis = [int(c[1]) for c in chunk]
        acc_array += [int(c[2]) for c in chunk]
        id_array += [int(c[3]) for c in chunk]
        idx_array += [str(f'{c[0]},{c[1]}') for c in chunk]

        ppts = stream[0, xis, yis]
        coords = tuple(map(tuple, zip(ppts.coords['x'].values, ppts.coords['y'].values)))
        point_array += [Point(p) for p in coords]

    df = pd.DataFrame()
    df['num_acc_cells'] = acc_array
    df['pt_id'] = id_array
    df['raster_idx'] = idx_array

    gdf = gpd.GeoDataFrame(df, geometry=point_array, crs=crs)
    print(f'{len(gdf)} pour points created.')   

    return gdf


def save_gdf(data):
    gdf, fpath = data
    gdf.to_file(fpath)


def create_confluence_point_vector(region, stream, crs, confluences, n_sample, method, separate_output_files=False):

    ppt_gdf = create_pour_point_gdf(stream, confluences, crs, method)

    output_folder = os.path.join(DATA_DIR, f'pour_points/{region}/')

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if separate_output_files:
        # if you want to create individual pour point files (shp)
        gdfs = []
        for i, _ in ppt_gdf.iterrows():
            ppt = ppt_gdf[ppt_gdf.index == i].copy()
            acc_val = ppt['num_acc_cells'].values[0]
            pt_id = ppt['pt_id'].values[0]
            gdf_path = os.path.join(output_folder, f'{region}_ppt_{method}_{i}_acc{acc_val}_id{pt_id}.shp')
            gdfs.append((ppt, gdf_path))

        with Pool() as p:
            p.map(save_gdf, gdfs)
    
    output_folder =  os.path.join(DATA_DIR, f'pour_points/{region}')
    output_path = os.path.join(output_folder, f'{region}_pour_pts_{method}.geojson')
    # output to a single gdf geojson
    ppt_gdf.to_file(output_path, driver='GeoJSON')



methods = ['CONF', 'GRAD']
# methods = ['CONF', 'GRAD']
# methods = ['RAND']

separate_output_files = False

cell_tracking_info = {}

for region in sorted(region_codes):
    # region = '07G'
    # region = 'Liard'
    # print('')

    stream, crs, affine = retrieve_raster(region, 'pruned_stream')

    S = stream.data[0]

    cell_tracking_info= {}

    tracking_df = pd.DataFrame()
    tracked_pcts = []
    for method in methods:
        print('')
        print(f'Processing {region} using {method} method.')

        region_area_km2 = get_region_area(region)

        ppt_sample_size = int(region_area_km2 / 10)
        # ppt_sample_size = 10

        print(f'Generating {ppt_sample_size} pour points for a region of {region_area_km2:.1f} km^2.')
        
        output_path = os.path.join(DATA_DIR, f'pour_points/')

        existing_files = os.listdir(output_path)

        if len(existing_files) == 5 * ppt_sample_size:
            print(f'   ...{ppt_sample_size} samples already exist.  Skipping {region}.')
        else:
            stream_px = np.where(S)

            t0 = time.time()

            if method == 'RAND':
                ppts, pct_cells_tracked = random_stream_point_selection(S, ppt_sample_size)
            elif method == 'CONF':
                ppts, pct_cells_tracked = find_stream_confluences_by_fdir(S, ppt_sample_size)
            elif method == 'GRAD':
                ppts, pct_cells_tracked = find_stream_confluences_by_acc_gradient(region, S, basin_threshold, ppt_sample_size)
            else:
                raise Exception(f'"{method}" method does not exist.  Variable "method" must be set to one of RND, NBR, or ACC.')

            t1 = time.time()
            print(f'Time to find pour point sample: {t1-t0:.1f}s')

            pntr_path = create_confluence_point_vector(region, stream, crs, ppts, ppt_sample_size, method, separate_output_files)
    
            tracking_df[method] = [len(ppts), pct_cells_tracked]

    tracking_df.index = ['num_cells_tracked', 'pct_cells_tracked']
    if not os.path.exists(os.path.join(output_path, f'ppt_stats/')):
        os.makedirs(os.path.join(output_path, f'ppt_stats/'))
    tracking_df.to_csv(os.path.join(output_path, f'ppt_stats/{region}_ppt_stats.csv'))

    print(asdfsad)

    
    
