
"""
Created on Thu May 16 11:57:00 2019

@author: gunterl
"""

from netCDF4 import Dataset
import numpy as np
import os, sys




def locate(x, y, lat, lon, itrans, ihem, xydist, rearth, ecc, slat):

# This function transforms I,J coordinates of an SSM/I grid cell to latitutde 
# and longitude coordinate. The programs can perform the inverse transformation 
# as well. This function interfaces with 2 other functions MAPXY and MAPLL.
#
# Note: This code is adapted from the fortran version of the subroutine 
#       in the NSIDC website Written by  V.J.Troisi - January, 1990
#       Updated by  N.A.Sandoval - November, 1995 
#       explanations: https://nsidc.org/data/polar-stereo/ps_grids.html and
#                     https://nsidc.org/data/polar-stereo/tools.html
#       obtained here: "ftp://sidads.colorado.edu/pub/DATASETS/seaice/polar-stereo/tools/"
#
#       
#
# input:
#   x,y    : Reals used when trans = 1. These entries describe 
#            the lower left coodinates (km) of a cell in an SSM/I grid.
#   lat,lon: Reals used when trans = 2. These reals describe the latitude and 
#            longitude (degrees) in an SSM/I grid which the function locate 
#            will transform to an I,J grid cell position.
#            Note: All latitudes and longitudes must be entered as positive numbers.
#   itrans : Integer describing the type of transformation performed by locate.
#            (1=1=I,J-to-Lat,Lon; 2=Lat,Lon-to-I,J)
#   ihem   : Integer describing one of the two polar regions (1=North , 2=South)
#   xydist : Distance from the origin of the grid in the cartesian plane.
#            The x-y coordinates for the edge of the lower left pixel
#   rearth : Earth's radius (km)
#   ecc    : eccentricity of the ellipsoid
#   slat   : Standard latitude for the SSM/I grids (typically 70 or 71)


    ##################
    #### Constants ###
    ##################
    PI = np.pi    # pi = 3.14...
    
    
    
    ##################
    #### Main body ###
    ##################
    
    
    # Define the sign and meridian offset (delta) for the SSM/I grids.
    if (ihem == 1):
        # Northern hemisphere (NH)
        SGN   = 1.0
        delta = 45.
    else:
        # Southern hemisphere (SH)
        SGN   = -1.0
        delta = 0.0
       

    # Beginning the translation:
    #
    # Case 1, form x-y to lat-lon
    if (itrans == 1):
        # In this case the lat-lon function entry can be ignored and will be 
        # rewritten in the transformation process. 

#         Convert I,J pairs to x and y distances from origin.
#         note: For some image display programs, the grid will be flipped in
#         the Y direction. Change j for i to be consistent. 
#        x=((j-1)*res-(xydist(1,1)-res/2.)
#        kk=numy(ihem,gtype)-(i-1)           
#        y=((kk-1)*res-(xydist(2,1)-res/2.)

        # Transform x and y distances to lat-lon.
        (lat,lon) = mapxy(x,y,slat,SGN,ecc,rearth)

        # Transform radians to positive degrees.
        lon=lon*180./PI
        lat=lat*180./PI
        lon=lon-delta   


        # Convert longitude to positive degrees
        np.where(lon <= 0.0, lon+360, lon)
        np.where(lon >= 360.0, lon-360, lon)
        
        return (lat,lon)
        
        
    else:
        # Case 2, form lat-lon to x-y 

        # Read the latitude and longitude pair and transform to cartesian where
        # that pair is located.
        # Transform degrees to radians:
        lat = np.abs(lat)*PI/180.
        lon = (lon+delta)*PI/180.

        # Transform latitude and longitude to x and y distances from origin
        (x,y) = mapll(lat,lon,slat,SGN,ecc,rearth)

        return (x,y)




def mapxy(x,y,slat,SGN,ecc,rearth):

# This function converts from Polar Stereographic (X,Y) coordinates to 
# geodetic latitude and longitude for the polar regions. The equations are
# from Snyder, J. P., 1982,  Map Projections Used by the U.S. Geological 
# Survey, Geological Survey Bulletin 1532, U.S. Government Printing Office.
# See JPL Technical Memorandum 3349-85-101 for further details.     
#
# Note: This function is adapted from the one revised by V. J. Troisi - January 1990
#       explanations: https://nsidc.org/data/polar-stereo/ps_grids.html and
#                     https://nsidc.org/data/polar-stereo/tools.html
#       obtained here: "ftp://sidads.colorado.edu/pub/DATASETS/seaice/polar-stereo/tools/"
#
#
# Inputs:
#   x,y    : real Polar stereographic X,Y coordinates (km)
#   slat   : Standard latitude for the SSM/I grids (typically 70 or 71)
#   SGN    : Sign for SSM/I grid (1 = NH, -1 = SH)
#   ecc    : eccentricity of the ellipsoid
#   rearth : Earth's radius (km)

    ##################
    #### Constants ###
    ##################
    PI = np.pi      # pi = 3.14...
    E2 = ecc**2     # square of eccentricity


    ##################
    #### Main body ###
    ##################

    SL = slat*PI/180.           # Converting from deg to rad
    RHO = np.sqrt(x**2+y**2)    # Radius computation
    
    
    
    # Case when not at the pole or vicinity
    CM = np.cos(SL)/np.sqrt(1.0-E2*(np.sin(SL)**2))
    T  = np.tan((PI/4.0)-(SL/(2.0)))/((1.0-ecc*np.sin(SL))/(1.0 + ecc*np.sin(SL)))**(ecc/2.0)

    if (np.abs(slat-90.) <= 1.e-5):
        T = RHO*np.sqrt((1.+ecc)**(1.+ecc)*(1.-ecc)**(1.-ecc))/2./rearth
    else:
        T = RHO*T/(rearth*CM)

    CHI = (PI/2.0)-2.0*np.arctan(T)          # Conformal latitude
    # Approximation version of the implicit equation:
    # lat = 2 arctan(tan(PI/4. + CHI/2)*((1+ecc*sin(lat))/(1-ecc*sin(lat)))^(ecc/2)) - PI/2
    lat = (CHI+((E2/2.0)+(5.0*E2**2.0/24.0)+(E2**3.0/12.0))*np.sin(2*CHI) +
          ((7.0*E2**2.0/48.0)+(29.0*E2**3/240.0))*np.sin(4.0*CHI) +
          (7.0*E2**3.0/120.0)*np.sin(6.0*CHI))
    lat = SGN*lat
    lon = np.arctan2(x*SGN,-y*SGN)
    lon = SGN*lon
        
    mask = RHO>0.1    # condition for coordinates to be away enough from pole.
                
    # Case at the pole or vicinity
    lat[~mask] = 90.*SGN
    lon[~mask] = 0.0

      
    return (lat,lon)




def mapll(lat,lon,slat,SGN,ecc,rearth):

# This function converts from geodetic latitude and longitude to Polar 
# Stereographic (X,Y) coordinates to for the polar regions. The equations are
# from Snyder, J. P., 1982,  Map Projections Used by the U.S. Geological 
# Survey, Geological Survey Bulletin 1532, U.S. Government Printing Office.
# See JPL Technical Memorandum 3349-85-101 for further details.     
#
# Note: This function is adapted from the one revised by Xiaoming Li - October 1996 
#       (Corrected equation for RHO)    
#       explanations: https://nsidc.org/data/polar-stereo/ps_grids.html and
#                     https://nsidc.org/data/polar-stereo/tools.html
#       obtained here: "ftp://sidads.colorado.edu/pub/DATASETS/seaice/polar-stereo/tools/"
#
#
# Inputs:
#   lat   : real geodetic latitude (deg, +90 to -90)
#   lon   : real geodetic longitude (deg, 0 to 360)
#   slat  : Standard latitude for the SSM/I grids (typically 70 or 71)
#   SGN   : Sign for SSM/I grid (1 = NH, -1 = SH)
#   ecc   : eccentricity of the ellipsoid
#   rearth: Earth's radius (km)

    ##################
    #### Constants ###
    ##################
    PI = np.pi      # pi = 3.14...
    E2 = ecc**2     # square of eccentricity


    ##################
    #### Main body ###
    ##################

    # Compute X and Y in grid coordinates

    T = np.tan(PI/4.-lat/2.)/((1.-ecc*np.sin(lat))/(1.+ecc*np.sin(lat)))**(ecc/2.)
    if (np.abs(90. - slat) < 1.e-5):
        RHO = 2.*rearth*T/((1.+ecc)**(1.+ecc)*(1.-ecc)**(1.-ecc))**(1/2.)
    else:
        SL  = slat*PI/180.
        TC  = np.tan(PI/4.-SL/2.)/((1.-ecc*np.sin(SL))/(1.+ecc*np.sin(SL)))**(ecc/2.)
        MC  = np.cos(SL)/np.sqrt(1.0-E2*(np.sin(SL)**2))
        RHO = rearth*MC*T/TC

    y =-RHO*SGN*np.cos(SGN*lon)
    x = RHO*SGN*np.sin(SGN*lon)

    # Condition if x,y are pole coordinates.    
    mask = (np.abs(lat) < PI/2.)

    x[~mask] = 0.0
    y[~mask] = 0.0


    return (x,y)


