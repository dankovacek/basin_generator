# generate basins


import os
import time
import glob
import math

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import geopandas as gpd
import pandas as pd
import rioxarray as rxr

import richdem as rd

import random

import multiprocessing as mp

from numba import jit

REV = 'RC'

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'processed_data/')

DATA_DIR = '/media/danbot/Samsung_T5/geospatial_data/basin_generator/'

DEM_resolution = 90 # EarthEnv DEM is 90m resolution
basin_threshold = int(1E6 / (90 * 90)) # min number of cells comprising a basin

processed_dem_dir = os.path.join(DATA_DIR, 'processed_dem/EENV_DEM/')
derived_basins_dir = os.path.join(DATA_DIR, f'{REV}/derived_basins')

processed_regions = sorted(list(set([e.split('_')[0] for e in os.listdir(processed_dem_dir)])))

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

@jit(nopython=True)
def process_slope_and_aspect(E, el_px, resolution, shape):
    dx, dy = resolution
    S = np.empty_like(E)
    S[:] = np.nan # track slope (in degrees)
    tot_p, tot_q = 0, 0
    for i, j in el_px:
        if (i == 0) | (j == 0) | (i == shape[0]) | (j == shape[1]):
            continue
            
        E_w = E[i - 1:i+2, j-1:j+2]

        a = E_w[0,0]
        b = E_w[1,0]
        c = E_w[2,0]
        d = E_w[0,1]
        f = E_w[2,1]
        g = E_w[0,2]
        h = E_w[1,2]
        # skip i and j because they're already used
        k = E_w[2,2]        

        all_vals = np.array([a, b, c, d, f, g, h, k])

        val_check = np.isfinite(all_vals)

        if np.all(val_check):
            p = ((c + 2*f + k) - (a + 2*d + g)) / (8 * abs(dx))
            q = ((c + 2*b + a) - (k + 2*h + g)) / (8 * abs(dy))
            cell_slope = np.sqrt(p*p + q*q)
            S[i, j] = (180 / np.pi) * math.atan(cell_slope)
            tot_p += p
            tot_q += q

    return S, tot_q, tot_p



def calculate_slope_and_aspect(clipped_raster):  
    """Calculate mean basin slope and aspect 
    according to Hill (1981).

    Args:
        clipped_raster (array): dem raster

    Returns:
        slope, aspect: scalar mean values
    """
    
    wkt = clipped_raster.rio.crs.to_wkt()
    affine = clipped_raster.rio.transform()
    resolution = clipped_raster.rio.resolution()
    raster_shape = clipped_raster.rio.width, clipped_raster.rio.height    

    rdem_clipped = rd.rdarray(
        clipped_raster.data[0], 
        no_data=clipped_raster.rio.nodata, 
        # projection=wkt, 
    )

    rdem_clipped.geotransform = affine.to_gdal()
    rdem_clipped.projection = wkt

    # ts0 = time.time()
    # slope = rd.TerrainAttribute(rdem_clipped, attrib='slope_degrees')
    # aspect_deg = rd.TerrainAttribute(rdem_clipped, attrib='aspect')
    # ts2 = time.time()

    el_px = np.argwhere(np.isfinite(rdem_clipped))
    S, tot_q, tot_p = process_slope_and_aspect(rdem_clipped, el_px, 
    resolution, raster_shape)

    mean_slope_deg = np.nanmean(S)
    # should be within a hundredth of a degree or so.
    # print(f'my slope: {mean_slope_deg:.4f}, rdem: {np.nanmean(slope):.4f}')

    mean_aspect_rad =  math.atan2(tot_q, -tot_p)
    mean_aspect_deg = (180.0 / np.pi) * mean_aspect_rad

    if(mean_aspect_deg < 0):
        mean_aspect_deg = 90 - mean_aspect_deg
    elif(mean_aspect_deg > 90.0):
        mean_aspect_deg = 360.0 - mean_aspect_deg + 90.0
    else:
        mean_aspect_deg = 90.0 - mean_aspect_deg

    # ts4 = time.time()

    # print(f'richdem slope + aspect: {ts2-ts0:.2f}')
    # print(f' custom slope + aspect: {ts4 - ts2:.2f}')

    return mean_slope_deg, mean_aspect_deg


def get_simulation_sample_filenames(region, method):
    all_files = os.listdir(f'{derived_basins_dir}/{region}/')
    return [e for e in all_files if f'basins_{method}' in e]


def main():
    methods = ['RAND', 'CONF', 'GRAD']

    ppt_sample_size = 10
    take_subsample = False

    for region in sorted(processed_regions):
        print(f'Processing {region} attributes.')
        for method in methods:
            print(f'   ...starting {method} method extraction.')            

            output_data_folder = os.path.join(BASE_DIR, f'processed_data/basin_properties/{REV}/{region}/')

            if not os.path.exists(output_data_folder):
                os.makedirs(output_data_folder)

            existing_files = os.listdir(output_data_folder)

            sample_files = get_simulation_sample_filenames(region, method)

            processed_sims = [e.split('_')[3] + '_' + e.split('_')[4].split('.')[0] for e in existing_files]

            processed_files = [f'{region}_basins_{s}.geojson' for s in processed_sims]

            # filtered_samples = [e for e in sample_files if e not in processed_files]
            filtered_samples = sample_files

            if len(filtered_samples) == 0:
                continue
            
            sim_times = []
            n = 0
            for sample_file in filtered_samples:

                simulation_no = sample_file.split('.')[0].split('_')[-1]

                sample_path = os.path.join(derived_basins_dir, region, sample_file)

                t_start = time.time()
                print(f'Processing {region}.')

                basin_df = gpd.read_file(sample_path)
                basin_df['Flags'] = ''

                basin_df[['mean_el', 'min_el', 'Elevation_m', 'Gravelius', 'Perimeter', 'Slope_deg', 'Aspect_deg']] = None

                if take_subsample:
                    sample_idx = random.choices(basin_df.index, k=ppt_sample_size)
                    basin_df = basin_df.loc[sample_idx, :].copy().sort_index()
                    print(f'   ...processing {len(basin_df)} basins.')

                dem_path = os.path.join(processed_dem_dir, f'{region}_EENV_DEM_3005.tif')
                region_raster, _, _  = retrieve_raster(dem_path)
                
                output_fpath = os.path.join(output_data_folder, f'{region}_basin_properties_{method}_{simulation_no}.csv')
                
                for i, row in basin_df.iterrows():
                    
                    basin_polygon = basin_df[basin_df.index == i]

                    clipped_raster, raster_loaded = clip_raster_to_basin(basin_polygon, region_raster)
                    
                    if not raster_loaded:
                        continue
                                       
                    mean_el, median_el, min_el, max_el = process_basin_elevation(clipped_raster)
                    
                    basin_df.loc[i, 'Drainage_Area_km2'] = row['geometry'].area / 1E6
                    basin_df.loc[i, 'mean_el'] = mean_el
                    # median elevation, but label it the same as hysets column (Elevation_m)
                    basin_df.loc[i, 'min_el'] = min_el
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

                n += 1

                t_end = time.time()
                unit_time = (t_end-t_start) / len(basin_df)
                sim_times.append(t_end - t_start)
                remaining_files = len(filtered_samples) - n
                t_remaining = np.mean(sim_times) * remaining_files
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