"""
Created on Thu Dec 17 16:14:03 2020

This script generates a SCRIP file from a cartesian CISM grid. Such a grid file could be generated using the notebook
"CISM_grid_create.ipynb".
A SCRIP file could be generated directly from a lat-lon file from the CISM tools. However, the corner cells defined 
by Lat0-Lon0 are missing the edge corners. Hence we are generating the file from the cartesian coordinate definition. 

This script assumes the CISM grid file variable convention x1-y1 for cell centers and x0-y0 for cell corners. 

@author: gunterl
"""


import sys, os
sys.path.insert(1, '../lib/')

from locate_vectorize import locate
import shutil
import numpy as np
from netCDF4 import Dataset
from pathlib import Path
import constants



                        #################
                        # Configuration #
                        #################

IS = 'GrIS'  # ice sheet for which the SCRIP file will be generated


# Defining the file to write the CISM grid NetCDF file
pathlocal = os.getcwd()
srcPath = f"{pathlocal}/cism_grid/"
inputFile  = f"xy_CISM3_grid_{IS}_04000.nc"

dstPath = f"{srcPath}SCRIP_FILES/"
outputFile = f"SCRIPgrid_{IS}_4km_epsg3413_c160524.nc"

# Creating the path where the Lat-Lon file will be written if it does not exist
Path(dstPath).mkdir(parents=True, exist_ok=True)





                        #############
                        # Main body #
                        #############

srcFile = srcPath + inputFile
outFile = dstPath + outputFile

Rearth = constants.Rearth/1000.  # converting Earth's radius from m to km
ecc = constants.ecc
PI = constants.pi



if IS not in ['AIS', 'GrIS']:
    sys.exit(f"Chose an ice sheet option in the list ['GrIS', 'AIS']")


if IS in ['GrIS']:
    
    print('Reading Greenland info')
    slat   = 70.                   # standard latitude
    itrans = 1                     # Transformation from I,J to lat-lon
    ihem   = 1                     # Northern hemisphere
    xydist = 0


if IS in ['AIS']:
    
    print('Reading Antarctic info')
    slat   = 71.                   # standard latitude
    itrans = 1                     # Transformation from I,J to lat-lon
    ihem   = 2                     # Southern hemisphere
    xydist = 0


# Illustration of a CISM grid cell


#       center NW     center NE
#           .-----------.        
#           |           |
#           |           |
#           |     .     |
#           | (corner)  |
#           |           |
#           .-----------.
#       center SW     center SE

########################################
# Reading the BedMachine original file #
########################################

ncidin = Dataset(srcFile, 'r')

xcenter = ncidin["x1"][:]
ycenter = ncidin["y1"][:]

create_x0 = 0
try:
    xcorner_temp = ncidin["x0"][:]
    ycorner_temp = ncidin["y0"][:]
except:
    create_x0 = 1
    print('need to create x0')


# coordinate dimensions
nx = len(xcenter)
ny = len(ycenter)

dx = xcenter[1]-xcenter[0]
dy = ycenter[1]-ycenter[0]


if create_x0:
    print('Creating x0')
    xcorner_temp = np.zeros(nx)
    ycorner_temp = np.zeros(ny)
    xcorner_temp = xcenter[0:-1] + dx/2.
    ycorner_temp = ycenter[0:-1] + dy/2.

# Extending the corner point so that each center points are included within a 
# defined cell. 
# Note: the initial corner is the lower left of CISM cartesian grid.

xcorner = np.zeros(nx+1)
ycorner = np.zeros(ny+1)

xcorner[1:-1] = xcorner_temp[:]
ycorner[1:-1] = ycorner_temp[:]

xcorner[0]  = xcorner_temp[0]-dx
xcorner[-1] = xcorner_temp[-1]+dx

ycorner[0]  = ycorner_temp[0]-dy
ycorner[-1] = ycorner_temp[-1]+dy

nxc = len(xcorner)
nyc = len(ycorner)


# Illustration of 4 CISM grid cells after extension


#         corner      corner      corner  
#           .-----------.-----------.        
#           |           |           |
#           |           |           |
#           |     .     |     .     |
#           |  center   |  center   |
#           |           |           |
#    corner .-----------.-----------. corner         
#           |           |           |
#           |           |           |
#           |     .     |     .     |
#           |  center   |  center   |
#           |           |           |
#           .-----------.-----------. 
#        corner       corner      corner





############################ 
# Creating the output file #
############################


print('Creating the output file')

# coordinate dimensions
nx = len(xcenter)
ny = len(ycenter)
nxc = len(xcorner)
nyc = len(ycorner)

grid_size = nx*ny

# Creating the netcdf output file
ncfile = Dataset(outFile, 'w')

# Create dimensions
ncfile.createDimension('grid_size', grid_size)
ncfile.createDimension('grid_corners', 4)
ncfile.createDimension('grid_rank', 2)


grid_dims            = ncfile.createVariable('grid_dims','i4',('grid_rank',))
grid_dims.long_name = "grid dimensions"
grid_dims.units     = "-"

grid_imask            = ncfile.createVariable('grid_imask','i4',('grid_size',))
grid_imask.long_name = "grid_imask"
grid_imask.units     = "-"

lat_center           = ncfile.createVariable('grid_center_lat','f8',('grid_size',))
lat_center.long_name = "grid_center_lat"
lat_center.units     = "degrees"

lon_center           = ncfile.createVariable('grid_center_lon','f8',('grid_size',))
lon_center.long_name = "grid_center_lat"
lon_center.units     = "degrees"

lat_corner           = ncfile.createVariable('grid_corner_lat','f8',('grid_size','grid_corners'))
lat_corner.long_name = "grid_center_lat"
lat_corner.units     = "degrees"

lon_corner           = ncfile.createVariable('grid_corner_lon','f8',('grid_size','grid_corners'))
lon_corner.long_name = "grid_center_lat"
lon_corner.units     = "degrees"

#grid_area            = ncfile.createVariable('grid_area','f4',('grid_size',))
#grid_area.long_name = "grid_area"
#grid_area.units     = "radians^2"


# Global attributes
ncfile.author = 'Gunter Leguy'
if IS in ['AIS']:
    ncfile.title = 'CISM EPSG:3031 grid'
if IS in ['GrIS']:
    ncfile.title = 'CISM EPSG:3413 grid'
ncfile.earth_radius = Rearth 
ncfile.eccentricity = ecc
ncfile.standard_parallel = -(-1)**itrans * slat
ncfile.false_easting = 0. ;
ncfile.false_northing = 0. ;
ncfile.note = "Used the same projection properties as those used in BedMachine2 from Morlighem M. et al., (2019)"


grid_dims[0] = int(nx)
grid_dims[1] = int(ny)

grid_imask[:] = 1

# Creating the cartesian grids
Xgridcenter, Ygridcenter = np.meshgrid(xcenter/1000.,ycenter/1000., sparse=False, indexing='ij')
Xgridcorner, Ygridcorner = np.meshgrid(xcorner/1000.,ycorner/1000., sparse=False, indexing='ij')


nx = len(xcenter)
ny = len(ycenter)
nxc = len(xcorner)
nyc = len(ycorner)

# Creating the lat-lon grid
LAT2d = np.zeros((nx,ny))
LON2d = np.zeros((nx,ny))
X2d = np.zeros((nx,ny))
Y2d = np.zeros((nx,ny))

xc = np.zeros(nx*ny)
yc = np.zeros(nx*ny)

for i in range(ny):

    xincenter = Xgridcenter[:,i] 
    yincenter = Ygridcenter[:,i]  

    latcenter,loncenter = locate(xincenter, yincenter, 0, 0, itrans, ihem, xydist, Rearth, ecc, slat)
    lon_center[i*nx:(i+1)*nx] = loncenter
    lat_center[i*nx:(i+1)*nx] = latcenter


# Ordering of cell corners

#           3           2
#           .-----------.        
#           |           |
#           |           |
#           |           |
#           |           |
#           |           |
#           .-----------.
#          0            1    


    xincorner = Xgridcorner[:,i]  
    yincorner = Ygridcorner[:,i]  

    latcorner,loncorner = locate(xincorner, yincorner, 0, 0, itrans, ihem, xydist, Rearth, ecc, slat)


    lon_corner[i*nx:(i+1)*nx,0] = loncorner[0:-1]
    lat_corner[i*nx:(i+1)*nx,0] = latcorner[0:-1]

    lon_corner[i*nx:(i+1)*nx,1] = loncorner[1::]
    lat_corner[i*nx:(i+1)*nx,1] = latcorner[1::]

    if i > 0:
        # For each iteration, the top row of corner cell is the bottom row of the next iteration.
        # The top row of corner cells needs to be filled in separately. 

        lon_corner[(i-1)*nx:i*nx,2] = lon_corner[i*nx:(i+1)*nx,1]
        lat_corner[(i-1)*nx:i*nx,2] = lat_corner[i*nx:(i+1)*nx,1]

        lon_corner[(i-1)*nx:i*nx,3] = lon_corner[i*nx:(i+1)*nx,0]
        lat_corner[(i-1)*nx:i*nx,3] = lat_corner[i*nx:(i+1)*nx,0]


# Filling in the top row of corner coordinates
xincorner = Xgridcorner[:,nyc-1]  
yincorner = Ygridcorner[:,nyc-1]  

latcorner,loncorner = locate(xincorner, yincorner, 0, 0, itrans, ihem, xydist, Rearth, ecc, slat)


lon_corner[-nx::,3] = loncorner[0:-1]
lat_corner[-nx::,3] = latcorner[0:-1]

lon_corner[-nx::,2] = loncorner[1::]
lat_corner[-nx::,2] = latcorner[1::]


ncfile.close()





