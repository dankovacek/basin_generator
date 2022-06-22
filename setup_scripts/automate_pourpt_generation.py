# generate basins

from concurrent.futures import process
import os
import sys

import random

import time

import warnings
warnings.filterwarnings('ignore')

import numpy as np
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


def retrieve_raster(region, type):
    filename = f'{region}_EENV_DEM_3005_{type}.tif'
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


def random_stream_point_selection(stream, ppt_sample_size):
    """
    Randomly select stream pixels.
    """
    D = stream.data[0].copy()
    stream_px = np.argwhere( D == 1 )
    return random.choices(stream_px, k=ppt_sample_size)
    

def find_stream_confluences_by_neighbors(stream, ppt_sample_size):
    """Assume confluences are stream pixels connected adjacent or diagonally with more than two other stream pixels.

    Args:
        stream (2d matrix): Binary matrix of stream pixels (1), non-stream pixels (0)

    Returns:
        _type_: list of pixel indices corresponding to confluences
    """
    D = stream.data[0].copy()
    # the matrix is very sparse, retrieve just the indices
    # corresponding to (all) stream pixels
    stream_px = np.argwhere( D == 1 )
    # print(non_nan_px[0])
    confluences = []
    c = 0
    for (i, j) in stream_px:
        
        W = D[i-1:i+2, j-1:j+2].round(2)
        # get the number of stream cells.
        # The stream_px array is the indices of all stream pixels.
        # The number of neighboring cells is one less than the 
        # total numeric cells.
        # 0 neighboring cells: pit, not a stream (don't track)
        # 1: headwater / outlet (track)
        # 2: mid-stream (don't track)
        # 3+: confluence
        num_stream_cells = np.count_nonzero(~np.isnan(W))

        num_neighbors = num_stream_cells - 1
    
        if (num_neighbors > 2):
            c += 1
            confluences.append((i, j))

    if len(confluences) > ppt_sample_size:
        return random.choices(confluences, k=ppt_sample_size)



    print(f'{len(confluences):.1e}/{len(stream_px):.1e} stream pixels are points of interest.')

    return confluences


def find_stream_confluences_by_acc_gradient(region, stream, threshold, ppt_sample_size):
    """_summary_

    Args:
        region (_type_): _description_
        idxs (_type_): _description_
        threshold (_type_): _description_

    Returns:
        _type_: _description_
    """
    acc, crs, affine = retrieve_raster(region, 'accum')

    A = acc.data[0].copy()

    # A = np.where(A > threshold, A, np.nan)

    D = stream.data[0].copy()
    stream_px = np.argwhere( D == 1 )

    px_of_interest = []
    
    for i, j in stream_px:
        # get the accumulation value of the target cell
        cell_acc = A[i, j]
        # get the accumulation values of a window 
        # of cells surrounding the target cell
        W = A[i-1:i+2, j-1:j+2].round(2)

        # num_stream_cells = np.count_nonzero(~np.isnan(W))

        # calculate the acc gradient w.r.t. the target cell
        acc_grad = cell_acc - W
        jump_px = np.argwhere( acc_grad > threshold )

        num_jumps = len(jump_px)

        if (num_jumps > 0):
            px_of_interest.append((i, j))

    if len(px_of_interest) > ppt_sample_size:
        return random.choices(px_of_interest, k=ppt_sample_size)

    return px_of_interest


def generate_pour_points(stream, confluences, n_chunks=2):
    
    n_chunks = int(10 * np.log(len(confluences)))

    conf_chunks = np.array_split(np.asarray(confluences), indices_or_sections=n_chunks)

    point_array = []
    for chunk in conf_chunks:
        xis = [c[0] for c in chunk]
        yis = [c[1] for c in chunk]

        ppts = stream[0, xis, yis]
        coords = tuple(map(tuple, zip(ppts.coords['x'].values, ppts.coords['y'].values)))
        point_array += [Point(p) for p in coords]
    
    return point_array


def save_gdf(data):
    gdf, fpath = data
    gdf.to_file(fpath)


def create_confluence_point_vector(region, stream, crs, confluences, n_sample, method, separate_output_files=False):

    ppts = generate_pour_points(stream, confluences)

    ppts = random.choices(ppts, k=n_sample)

    gdf = gpd.GeoDataFrame(geometry=ppts, crs=crs)
    print(f'{len(gdf)} pour points created.')

    output_folder = os.path.join(DATA_DIR, f'pour_points/{region}/{method}')

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if separate_output_files:
        # if you want to create individual pour point files (shp)
        gdfs = []
        for i, _ in gdf.iterrows():
            ppt = gdf[gdf.index == i].copy()
            gdf_path = os.path.join(output_folder, f'{region}_ppt_{method}_{i}.shp')
            gdfs.append((ppt, gdf_path))

        with Pool() as p:
            p.map(save_gdf, gdfs)
    
    output_folder =  os.path.join(DATA_DIR, f'pour_points/{region}')
    output_path = os.path.join(output_folder, f'{region}_pour_pts_{method}_N{n_sample}_RA.geojson')
    # output to a single gdf geojson
    gdf.to_file(output_path, driver='GeoJSON')



method = 'RND'
method = 'NBR'
method = 'ACC'

separate_output_files = True

for region in ['07G']:#sorted(region_codes):
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
        stream, crs, affine = retrieve_raster(region, 'pruned_stream')

        stream_px = np.where(stream)

        t0 = time.time()

        if method == 'RND':
            ppts = random_stream_point_selection(stream, ppt_sample_size)
        elif method == 'NBR':
            ppts = find_stream_confluences_by_neighbors(stream, ppt_sample_size)
        elif method == 'ACC':
            ppts = find_stream_confluences_by_acc_gradient(region, stream, basin_threshold, ppt_sample_size)
        else:
            raise Exception(f'"{method}" method does not exist.  Variable "method" must be set to one of RND, NBR, or ACC.')

        t1 = time.time()
        print(f'Time to find pour point sample: {t1-t0:.1f}s')

        pntr_path = create_confluence_point_vector(region, stream, crs, ppts, ppt_sample_size, method, separate_output_files)

        print(asfdsadf)
