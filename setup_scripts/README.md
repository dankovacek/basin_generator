Automated Basin Delineation
===========================

This repository provides an automated pipeline for generating large
samples of basins from DEM. First, DEM files are collected from an
open-source repository, then assembled into a raster tile mosaic. The
study region (British Columbia) is broken into hydraulically “complete”
sub-regions, defined by boundaries featuring only outflow. The polygons
describing these sub-regions are used to create clipped rasters, which
are then hydraulically enforced (fill depressions and resolve flats) to
create flow direction, flow accumulation, and stream rasters. The stream
raster is used as a binary mask to identify a set of pour points to use
for the final step, basin delineation.

The resulting collection of basins is an estimate of the population, but
the basin delineation process raises interesting questions about the
representativeness of the sample. In the DEM preprocessing steps, we use
topography to identify stream networks. Stream networks are represented
by cells that meet some criteria of minimum flow accumulation. That is,
surface runoff will only accumulate if there is enough area upstream to
contribute runoff. There is no single number to represent this minimum
threshold, so we make an estimate and possibly even calibrate it with
satellite imagery, orthophoto, or ground-survey of stream headwaters.

1.  How many basins do we need for a representative sample? We’ve seen
    smaller samples do OK in approximating the population. Could be used
    to reduce the sampling burden.

Data Acquisition and Preprocessing
----------------------------------

Set up Computing Environment
----------------------------

**Note, this code was tested on Ubuntu Linux.** These instructions are
intended to minimize setup time and compute resources.

<!-- >`sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable`   -->

> `sudo apt update`  
> `sudo apt upgrade`

Install dependencies:
&gt;`sudo apt-get install gdal-bin libgdal-dev libproj15 libproj-dev openmpi-bin libopenmpi-dev libboost-iostreams-dev parallel unzip dos2unix zip`

**Clone the repository (from the root directory)**

> `git clone https://github.com/dankovacek/basin_generator`

Change directories to the `basin_generator` folder:  
&gt;`cd basin_generator`

### Install Python package manager (pip)

If not automatically installed, install Python and virtualenv (assuming
Python3 is installed by default on a linux distribution):

> `sudo apt install python3.8-venv pip`

Create Python 3.8+ virtual environment at the root level directory:

> `python3 -m venv env/`

Activate the virual environment:  
&gt;`source env/bin/activate`

Install Python packages:  
&gt;`pip install -r requirements.txt`

Data Acquisition and Processing
-------------------------------

### EarthEnv DEM90 (Robinson, Regetz, and Guralnick 2014)

HYSETS used EarthEnv DEM90 (~90m resolutoin) data to derive basins and
physical / topographical attributes. The tiles are available at
[earthenv.org/DEM](https://earthenv.org/DEM). We can download the tiles
and merge them into a virtual raster with gdal using the following
script in `setup_scripts`:  
&gt;`python get_EarthEnv_DEM.py`

Links to invidivual DEM tiles look like the following:  
&gt;`http://mirrors.iplantcollaborative.org/earthenv_dem_data/EarthEnv-DEM90/EarthEnv-DEM90_N55W125.tar.gz`

The resulting .vrt mosaic should look like below (red outline added to
emphasize land mass bounds.):

![DEM Mosaic of BC and administrative boundary
regions](../img/wsc-ssda-processed.png)

DEM Preparation
---------------

This step represents the heavy lifting where large regions of DEM such
as the Liard, Peace, and and Fraser River basins are processed into flow
direction and flow accumulation.
[Whiteboxtools](python%20process_dem_by_basin.py) was used here for the
step of delineating a large set of basins.

<!-- Note: the breach [depression function](https://jblindsay.github.io/ghrg/Whitebox/Help/BreachDepressions.html) run on the DEM is a bottleneck step.   -->

Create the individual region DEM files using the provided region
polygons and the DEM tile mosaic created in the previous step:

> `cd setup_scripts/`  
> `python create_complete_region_DEMS.py`

Process the region DEMs to create rasters representing flow direction,
flow accumulation, and stream network:  
&gt;`python derive_flow_accumulation.py`

Using the stream raster, generate pour points at headwaters and
confluences:  
&gt;`python automate_pourpt_generation.py`

Generate a basin for each of the pour points:  
&gt;`setup_scripts/python pysheds_derive_basin_polygons.py`

Additional Notes
----------------

<!-- Automate citation formatting for the README document.

>`pandoc -t markdown_strict -citeproc README-draft.md -o README.md --bibliography bib/bibliography.bib` -->

Robinson, Natalie, James Regetz, and Robert P Guralnick. 2014.
“EarthEnv-Dem90: A Nearly-Global, Void-Free, Multi-Scale Smoothed, 90m
Digital Elevation Model from Fused Aster and Srtm Data.” *ISPRS Journal
of Photogrammetry and Remote Sensing* 87: 57–67.
