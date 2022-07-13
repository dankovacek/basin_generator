import os
import sys
import zipfile

import time

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import geopandas as gpd
import rioxarray as rxr

from shapely.geometry import Polygon

from shapely.validation import make_valid

# specify the DEM source
DEM_source = 'EENV_DEM'

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_DIR = os.path.join(BASE_DIR, 'data/')

# this is a custom path because I want my files saved to an external disk
DATA_DIR = '/media/danbot/Samsung_T5/geospatial_data/basin_generator/data/'

mask_dir = os.path.join(BASE_DIR, f'processed_data/merged_basin_groups/')

EENV_DEM_DIR = os.path.join(DATA_DIR, DEM_source)

assert os.path.exists(EENV_DEM_DIR)

# this command builds the dem mosaic "virtual raster"
mosaic_path = os.path.join(DATA_DIR, f'{DEM_source}_mosaic_4326.vrt')

vrt_command = f"gdalbuildvrt -resolution highest -a_srs EPSG:4326 {mosaic_path} {EENV_DEM_DIR}/EarthEnv-DEM90_*.bil"

os.system(vrt_command)

print(asfsd)

t0 = time.time()

epsg_code = 4326
vrt_file = f'{DEM_source}_mosaic_4326.vrt'
if DEM_source == 'USGS_3DEP':
    epsg_code = 4269
    vrt_file = f'processed_dem/{DEM_source}_mosaic_{epsg_code}.vrt'

dem_mosaic_file = os.path.join(EENV_DEM_DIR, vrt_file)

def get_crs_and_resolution(fname):
    raster = rxr.open_rasterio(fname)
    crs = raster.rio.crs.to_epsg()
    res = raster.rio.resolution()   
    return crs, res


dem_crs, (w_res, h_res) = get_crs_and_resolution(dem_mosaic_file)


def check_mask_validity(mask_path):

    mask = gpd.read_file(mask_path)
    gtype = mask.geometry.values[0].geom_type
    if gtype == 'GeometryCollection':
        raise Exception; 'geometry should not be GeometryCollection'

    if mask.geometry.is_valid.values[0]:
        print(f'   ...mask is valid.')
        # mask.to_file(mask_path, driver='GeoJSON')
    else:
        file = mask_path.split('/')[-1]
        print(f'   ...{file} mask is invalid:')
        mask['geometry'] = mask.geometry.apply(lambda m: make_valid(m))
        mask = mask.dissolve()
        # drop invalid geometries
        mask = mask[mask.geometry.is_valid]
        mask['area'] = mask.geometry.area
        # reproject to 4326 to correspond with DEM tile mosaic
        mask = mask.to_crs(4326)
        
        fixed = all(mask.geometry.is_valid)
        
        if fixed:
            mask.to_file(mask_path, driver='GeoJSON')
            print(f'   ...invalid mask corrected: {fixed}')
            # mask.to_file('DEM_treatment_test/08NM071_original.geojson')
            # print(asdfsad)
            # mask = mask.to_crs(3005)
            # mask = mask.buffer(1000, resolution=16, single_sided=True)
            # m = '08NM071_buffered.geojson'
            # mask.to_file(f'DEM_treatment_test/{m}')
        else:
            print(f'   ...invalid mask could not be corrected')

if not os.path.exists(mask_dir + 'region_polygons'):
    with zipfile.ZipFile(mask_dir + 'region_polygons.zip', 'r') as zip_ref:
        zip_ref.extractall(mask_dir + 'region_polygons')

mask_dir = os.path.join(mask_dir, 'region_polygons/')
all_masks = [e for e in os.listdir(mask_dir) if e.endswith('.geojson')]

all_codes = [e.split('_')[0] for e in all_masks]
# all_codes = [e for e in all_codes if e in ['07U', 'ERockies', '08N']]

i = 0
for code in all_codes:

    file = [e for e in all_masks if code in e][0]
    
    fpath = os.path.join(mask_dir, file)

    if '_' not in file:
        splitter = '.'
    else:
        splitter = '_'
    grp_code = file.split(splitter)[0]
    print(f'Starting polygon merge on {grp_code}.')

    
    named_layer = file.split('.')[0]

    mask_check = check_mask_validity(fpath)
    

    # set the output initial path and reprojected path
    out_path = f'{EENV_DEM_DIR}{grp_code}_{DEM_source}_{epsg_code}.tif'
    out_path_reprojected = f'{EENV_DEM_DIR}{grp_code}_{DEM_source}_3005.tif'

    # if you want to modify the resulting DEM resolution,
    # you'll need to replace 1 with some other factor here
    rfactor = 1
    trw = abs(w_res*rfactor)
    trh = abs(h_res*rfactor)

    if not os.path.exists(out_path_reprojected):

        if rfactor == 1:
            command = f'gdalwarp -s_srs epsg:{epsg_code} -cutline {fpath} -cl {named_layer} -crop_to_cutline -multi -of gtiff {dem_mosaic_file} {out_path} -wo NUM_THREADS=ALL_CPUS'
        else:
            command = f'gdalwarp -s_srs epsg:4326 -cutline {fpath} -cl {named_layer} -ot Float32 -crop_to_cutline -tr {trw} {trh} -multi -of gtiff {dem_mosaic_file} {out_path} -wo CUTLINE_ALL_TOUCHED=TRUE -wo NUM_THREADS=ALL_CPUS'
        print('')
        print('__________________')
        print(command)
        print('')
        try:
            os.system(command)
        except Exception as e:
            raise Exception; e
    else:
        fname = out_path_reprojected.split('/')[-1]
        print(f'   ...{fname} exists, skipping dem cutline operation..')

    
    # check # pixels low res        
    if not os.path.exists(out_path_reprojected):
        # reproject to epsg 3005
        lr = rxr.open_rasterio(out_path, masked=True, default='dem')
        lr = lr.rio.reproject(3005)
        lr.rio.to_raster(out_path_reprojected)
        lr_shape = lr.rio.shape
        n_pix = lr_shape[0] * lr_shape[0]
        print(f'   ...img has {n_pix:.2e} pixels')
        os.remove(out_path)
    else:
        fname = out_path_reprojected.split('/')[-1]
        print(f'   ...{fname} exists, skipping dem reprojection..')
    
    t1 = time.time()
    print(f'      {i+1}/{len(all_masks)} Completed tile merge: {out_path_reprojected} created in {t1-t0:.1f}s.')
    print('')
    print('')
    i += 1
    