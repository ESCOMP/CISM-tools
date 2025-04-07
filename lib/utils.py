
"""
Created on April 02 2025

@author: gunterl
"""

from __future__ import annotations

import constants
from locate_vectorize import locate
import numpy as np
import os, sys




def create_latlon_coord(x, y, IS, R=constants.Rearth, ecc=constants.ecc):
    """
    This function uses stereographic projection to project CISM cartesian to lat-lon coordinates.
    It uses the assumption that CISM is defined on a regular and uniform grid. 
    Inputs:
        x   : 1d array of CISM x-axis coordinates (m)
        y   : 1d array of CISM y-axis coordinates (m)
        IS  : string selecting the ice sheet grid we'd like to project. Options are 'GrIS' or 'AIS'
        R   : Earth's radius (m) (set by default to the value defined in constants.py)
        ecc : Earth's eccentricity (set by default to the value defined in constants.py)        

    Outputs:
        lat2d : 2d-array of CISM latitude coordinates (deg)
        lon2d : 2d-array of CISM longitude coordinates (deg)        
    """

    if IS not in ['GrIS', 'AIS']:
        sys.exit(f"The ice sheet {IS} is not yet supported in this function. Choose an option in ['GrIS', 'AIS']")

    if IS in ['GrIS']:
        # Adding the Lat-Lon field to the output file
        itrans = 1      # option to project I-J to Lat-Lon
        ihem   = 1      # option for Southern hemisphere
        xydist = 0.     # Distance from the origin of the grid in the Cartesian plane
        slat   = 70.    # Standard latitude for the SSMM/I grids
    
        # Creating a meshgrid from the cartesian coordinates
        x2d, y2d = np.meshgrid(x, y)
        lat2d, lon2d = locate(x2d/1000., y2d/1000., 0, 0, itrans, ihem, xydist, R/1000., ecc, slat)
        # Location of the South pole
        ind = np.where(np.abs(lat2d) > 100)[:]
        if len(ind[0]>0):    
            lat2d[ind[0], ind[1]] = -90
            
    if IS in ['AIS']:
        # Adding the Lat-Lon field to the output file
        itrans = 1      # option to project I-J to Lat-Lon
        ihem   = 2      # option for Southern hemisphere
        xydist = 0.     # Distance from the origin of the grid in the Cartesian plane
        slat   = 71.    # Standard latitude for the SSMM/I grids
    
        # Creating a meshgrid from the cartesian coordinates
        x2d, y2d = np.meshgrid(x, y)
        lat2d, lon2d = locate(x2d/1000., y2d/1000., 0, 0, itrans, ihem, xydist, R/1000., ecc, slat)
        # Location of the South pole
        ind = np.where(np.abs(lat_2d) > 100)[:]
        if len(ind[0]>0):    
            lat2d[ind[0], ind[1]] = -90
    
    return lat2d, lon2d











