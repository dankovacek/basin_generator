# derive flow accumulation network from DEM
import os

import warnings
warnings.filterwarnings('ignore')


from whitebox.whitebox_tools import WhiteboxTools

wbt = WhiteboxTools()
wbt.verbose = False

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'processed_data/')

processed_data_dir = os.path.join(DATA_DIR, 'processed_dem')

region_files = os.listdir(os.path.join(processed_data_dir, 'EENV_DEM'))
region_codes = [e.split('_')[0] for e in region_files]

DEM_resolution = 90 # EENV DEM90 is approximately 90m resolution

# determine the threshold number of cells 
# corresponding to the minimum drainage area
minimum_basin_size = 5 # km^2
threshold = int(minimum_basin_size * 1E6 / (DEM_resolution**2))

for region in sorted(region_codes):
    print(f'processing flow direction and accumulation for {region}')
    dem_file = f'EENV_DEM/{region}_EENV_DEM_3005.tif'
    dem_path = os.path.join(processed_data_dir, dem_file)
    out_dem_file = f'EENV_DEM/{region}_EENV_DEM_3005_temp_dem.tif'
    out_pntr_file = f'EENV_DEM/{region}_EENV_DEM_3005_temp_fdir.tif'
    out_accum_file = f'EENV_DEM/{region}_EENV_DEM_3005_accum.tif'
    out_stream_file = f'EENV_DEM/{region}_EENV_DEM_3005_stream.tif'

    accum_path = os.path.join(processed_data_dir, out_accum_file)
    if not os.path.exists(accum_path):
        wbt.flow_accumulation_full_workflow(
            dem_path, 
            os.path.join(processed_data_dir, out_dem_file), 
            os.path.join(processed_data_dir, out_pntr_file), 
            os.path.join(processed_data_dir, out_accum_file), 
            out_type='cells', 
            log=False, 
            clip=False, 
            esri_pntr=False, 
        )

    stream_path = os.path.join(processed_data_dir, out_stream_file)
    if not os.path.exists(stream_path):
        wbt.extract_streams(
            accum_path, 
            stream_path, 
            threshold, 
            zero_background=False, 
        )
    out_pruned_stream_file = f'EENV_DEM/{region}_EENV_DEM_3005_pruned_stream.tif'
    pruned_stream_path = os.path.join(processed_data_dir, out_pruned_stream_file)
    if not os.path.exists(pruned_stream_path):
        wbt.remove_short_streams(
            os.path.join(processed_data_dir, out_pntr_file),
            stream_path, 
            pruned_stream_path, 
            100, # min tributary length in map units (metres) 
            esri_pntr=False, 
        )