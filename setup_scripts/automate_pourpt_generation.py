# generate basins

from concurrent.futures import process
import os
import sys

import random

import time
from turtle import clear

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

# if you are using an external SSD to store files (good idea)
EXT_MEDIA = '/media/danbot/Samsung_T5/geospatial_data/basin_generator/'
DATA_DIR = EXT_MEDIA

output_path = os.path.join(DATA_DIR, f'pour_points/')
if not os.path.exists(output_path):
    os.mkdir(output_path)

DEM_resolution = 90 # EarthEnv DEM is approximately 90m resolution

# this should correspond with the threshold accumulation
# set in "derive_flow_accumulation.py"
min_basin_area = 5 # km^2
# min number of cells comprising a basin
basin_threshold = int(min_basin_area * 1E6 / (90 * 90)) 

processed_data_dir = os.path.join(DATA_DIR, 'processed_dem')

region_files = os.listdir(os.path.join(processed_data_dir, 'EENV_DEM'))

region_codes = sorted(list(set([e.split('_')[0] for e in region_files])))

inflow_direction_matrix = np.array([
    [4, 8, 16],
    [2, 0, 32],
    [1, 128, 64]
])

d8_matrix = np.array([
    [64, 128, 1],
    [32, 0, 2],
    [16, 8, 4]
])

def retrieve_raster(region, raster_type):
    filename = f'{region}_EENV_DEM_3005_{raster_type}.tif'
    raster_path = os.path.join(processed_data_dir, f'EENV_DEM/{filename}')
    raster = rxr.open_rasterio(raster_path, mask_and_scale=True)
    crs = raster.rio.crs
    affine = raster.rio.transform(recalc=False)
    return raster, crs, affine


def get_region_area(region):
    polygon_path = os.path.join(BASE_DIR, 'processed_data/merged_basin_groups/region_polygons/')
    poly_files = os.listdir(polygon_path)
    file = [e for e in poly_files if e.startswith(region)]
    if len(file) != 1:
        raise Exception; 'Region shape file not found.'
    fpath = os.path.join(polygon_path, file[0])
    gdf = gpd.read_file(fpath)
    gdf = gdf.to_crs(3005)
    return gdf['geometry'].area.values[0] / 1E6


def random_stream_point_selection(stream_px, A):
    """
    Here we simply return all of the stream pixels.
    The random selection will occur later, but it is 
    faster to read the file once and do iterative 
    random selections instead of creating a temporary 
    gdf in the random selection step for every simulation.
    """
    ppts = []
    c = 0
    for i, j in stream_px:
        cell_acc = A[i, j]
        ppts.append((i, j, cell_acc, c))
        c += 1

    pct_tracked = len(ppts) / len(stream_px) * 100
    print(f'Tracked {c} randomly selected points.')
    print(f'{len(ppts):.3e}/{len(stream_px):.3e} stream pixels are points of interest ({pct_tracked:.0f}%).')
        
    return ppts, pct_tracked


def create_flow_direction_mask(S, F, i, j, fdir, filter_stream=True):    
    if filter_stream:
        # create a boolean matrix for stream cells
        S_w = np.where(S[i - 1:i+2, j-1:j+2] == 1, 1, 0)
    else:
        S_w = np.ones_like(fdir)

    # retrieve D8 pointer from a 3x3 window 
    # around the target pixel
    F_w = F[i-1:i+2, j-1:j+2]
    
    # mask the flow direction by element-wise multiplication
    # using the boolean stream matrix
    F_wf = np.multiply(F_w, S_w)

    # get the accumulation window to find the drainage area
    # of each inflow cell (as well as the target cell)
    # A_w = np.multiply(A[i-1:i+2, j-1:j+2], S_w)
    return np.where(F_wf == fdir, True, False)


def find_stream_confluences_by_fdir(stream_px, S, A, F):
    """Assume confluences are stream pixels connected adjacent or diagonally with more than two other stream pixels.

    Args:
        stream (2d matrix): Binary matrix of stream pixels (1), non-stream pixels (0)

    Returns:
        _type_: list of pixel indices corresponding to confluences
    """
    confluences = []
    c = 0
    for (i, j) in stream_px:
        # find flow direction cells
        # flowing towards the center cell.
        F_in = create_flow_direction_mask(S, F, i, j, inflow_direction_matrix)
        
        # set the focus cell as true so it gets added
        F_in[1, 1] = True
        # # num_matches = np.sum(W_w == inflow_direction_matrix)
        inflow_cells = ()
    
        if np.sum(F_in) > 2:
            c += 1
            # get the focus cell accumulation
            cell_acc = A[i, j]
            # get inflow cell indices
            inflow_cells = np.argwhere(F_in)
            for ci, cj in inflow_cells:
                ix = ci + i - 1
                jx = cj + j - 1
                cell_acc = A[ix, jx]  
                confluences.append((ix, jx, cell_acc, c))

    pct_found = len(confluences) / len(stream_px) * 100
    print(f'{c} confluences found.')
    print(f'{len(confluences):.3e}/{len(stream_px):.3e} stream pixels are points of interest ({pct_found:.0f}%).')

    return confluences, pct_found


def find_gradient_stream_cells(stream_px, S, A, F, threshold, method):
    """_summary_

    Args:
        region (_type_): _description_
        idxs (_type_): _description_
        threshold (_type_): _description_

    Returns:
        _type_: _description_
    """
    px_of_interest = []
    tracked_indices = []
    c = 0
    for i, j in stream_px:
        # get the accumulation value of the target cell
        focal_cell_acc = A[i, j]
                
        # retrieve flow accumulation from a 3x3 window 
        # around the target pixel
        A_w = A[i-1:i+2, j-1:j+2]

        # create a boolean matrix for stream cells
        S_w = np.where((S[i-1:i+2, j-1:j+2] == 1), True, False)
        
        # calculate the acc gradient w.r.t. the focal cell
        # and filter out non-stream cells
        A_g = np.abs(focal_cell_acc - A_w)

        # create a boolean matrix for cells that flow into the focal cell
        F_w = create_flow_direction_mask(S, F, i, j, inflow_direction_matrix, filter_stream=True)

        # flip the bit on the cell the focal cell points at
        # to count it in the gradient calculation. 
        focal_cell_direction = np.argwhere(F[i, j]==d8_matrix)
        dir_idx = focal_cell_direction[0]
        F_w[dir_idx[0],dir_idx[1]] = True

        # mask the flow direction by element-wise multiplication
        # using the boolean flow direction matrix.
        # this yields all stream cells flowing into focal cell
        A_gf = np.multiply(A_g, F_w)

        # GRAB - AB stands for absolute -- find all cells where
        # the accumulation difference is greater than a constant
        # number of cells (minimum area change of interest).
        threshold_crit = (A_gf > threshold).any()

        # proportion criteria not needed for the GRAB method
        proportion_crit = True
        if method == 'GRAP':
            # AP stands for area proportional -- find cells where
            # the difference is greater than 10% of the accumulation
            # of the focal cell
            proportion_threshold = 0.5 * focal_cell_acc
            proportion_crit = (A_gf > proportion_threshold).any()

        if threshold_crit & proportion_crit:            
            if [i, j] not in tracked_indices:
                px_of_interest.append((i, j, focal_cell_acc, c))
                tracked_indices.append([i, j])
                c += 1

    pct_found = len(px_of_interest) / len(stream_px) * 100
    print(f'{len(px_of_interest):.3e}/{len(stream_px):.3e} stream pixels are points of interest ({pct_found:.0f}%).')

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


def create_confluence_point_vector(region, stream, crs, confluences, method, separate_output_files=False):
    ta = time.time()
    ppt_gdf = create_pour_point_gdf(stream, confluences, crs, method)
    tb = time.time()
    print(f'   ...time to create ppt gdf: {tb-ta:.1f}')

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
    
    t0 = time.time()
    output_folder =  os.path.join(DATA_DIR, f'pour_points/{region}')
    output_path = os.path.join(output_folder, f'{region}_pour_pts_{method}.geojson')
    # output to a single gdf geojson
    ppt_gdf.to_file(output_path, driver='GeoJSON')
    t1 = time.time()
    print(f'   ...time to write ppt file: {t1-t0:.1f}')


methods = ['RAND', 'CONF', 'GRAB', 'GRAP']
methods = ['CONF', 'GRAB', 'GRAP']
# methods = ['GRAP']

separate_output_files = False

cell_tracking_info = {}

region_codes = ['07G']

for region in sorted(region_codes):

    region_area_km2 = get_region_area(region)

    ppt_sample_size = int(region_area_km2 / 100)
    # ppt_sample_size = 10

    print(f'Generating {ppt_sample_size} pour points for a region of {region_area_km2:.1f} km^2 to yield 1 station per 100 km^2.')

    rt0 = time.time()
    
    stream, crs, affine = retrieve_raster(region, 'pruned_stream')

    fdir, crs, affine = retrieve_raster(region, 'temp_fdir')

    acc, crs, affine = retrieve_raster(region, 'accum')

    rt1 = time.time()
    print(f'   ...time to load resources: {rt1-rt0:.1f}s.')

    S = stream.data[0]
    F = fdir.data[0]
    A = acc.data[0]

    stream_px = np.argwhere(S == 1)

    tracking_df = pd.DataFrame()
    cell_tracking_info= {}
    tracked_pcts = []
    for method in methods:
        print('')
        print(f'Processing {region} using {method} method.')

        t0 = time.time()

        if method == 'RAND':
            ppts, pct_cells_tracked = random_stream_point_selection(stream_px, A, ppt_sample_size)
        elif method == 'CONF':
            ppts, pct_cells_tracked = find_stream_confluences_by_fdir(stream_px, S, A, F)
        elif method.startswith('GRA'):
            ppts, pct_cells_tracked = find_gradient_stream_cells(stream_px, S, A, F, basin_threshold, method)
        else:
            raise Exception(f'"{method}" method does not exist.  Variable "method" must be set to one of RND, NBR, or ACC.')

        t1 = time.time()
        print(f'Time to find pour point sample: {t1-t0:.1f}s')

        pntr_path = create_confluence_point_vector(region, stream, crs, ppts, method, separate_output_files)

        tracking_df[method] = [len(ppts), pct_cells_tracked]

    tracking_df.index = ['num_cells_tracked', 'pct_cells_tracked']
    if not os.path.exists(os.path.join(output_path, f'ppt_stats/')):
        os.makedirs(os.path.join(output_path, f'ppt_stats/'))
    tracking_df.to_csv(os.path.join(output_path, f'ppt_stats/{region}_ppt_stats.csv'))

    
    
