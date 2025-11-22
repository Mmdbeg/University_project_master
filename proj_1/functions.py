import numpy as np
import math
import pandas as pd
import pyproj as pj
from pygeodesy.geoids import GeoidPGM
from pygeodesy.ellipsoidalKarney import LatLon  
import matplotlib.pyplot as plt

# Constants for WGS84
a = 6378137.0              # semi-major axis
f = 1 / 298.257223563
b = a * (1 - f)
e = math.sqrt(1 - (b / a)**2)
k0 = 0.9996  # scale factor

def utm_to_latlon(x, y, zone, northern=True):
    # UTM zone central meridian (in degrees)
    lon0 = (zone - 1) * 6 - 180 + 3  
    lon0_rad = math.radians(lon0)
    
    # Remove 500,000 meter offset for longitude
    x -= 500000.0
    
    # If in southern hemisphere, remove 10,000,000 meter offset
    if not northern:
        y -= 10000000.0
    
    # Calculate the meridional arc
    M = y / k0
    
    # Compute the footprint latitude
    mu = M / (a * (1 - e**2 / 4 - 3*e**4 / 64 - 5*e**6 / 256))
    e1 = (1 - math.sqrt(1 - e**2)) / (1 + math.sqrt(1 - e**2))
    
    J1 = (3*e1/2 - 27*e1**3/32)
    J2 = (21*e1**2/16 - 55*e1**4/32)
    J3 = (151*e1**3/96)
    J4 = (1097*e1**4/512)
    
    fp = mu + J1*math.sin(2*mu) + J2*math.sin(4*mu) + J3*math.sin(6*mu) + J4*math.sin(8*mu)
    
    # Compute latitude and longitude
    e2 = e**2 / (1 - e**2)
    C1 = e2 * math.cos(fp)**2
    T1 = math.tan(fp)**2
    N1 = a / math.sqrt(1 - e**2 * math.sin(fp)**2)
    R1 = N1 * (1 - e**2) / (1 - e**2 * math.sin(fp)**2)
    D = x / (N1 * k0)
    
    Q1 = N1 * math.tan(fp) / R1
    Q2 = D**2 / 2
    Q3 = (5 + 3*T1 + 10*C1 - 4*C1**2 - 9*e2) * D**4 / 24
    Q4 = (61 + 90*T1 + 298*C1 + 45*T1**2 - 252*e2 - 3*C1**2) * D**6 / 720
    lat = fp - Q1 * (Q2 - Q3 + Q4)
    
    Q5 = D
    Q6 = (1 + 2*T1 + C1) * D**3 / 6
    Q7 = (5 - 2*C1 + 28*T1 - 3*C1**2 + 8*e2 + 24*T1**2) * D**5 / 120
    lon = lon0_rad + (Q5 - Q6 + Q7) / math.cos(fp)
    
    # Convert from radians to degrees
    lat_deg = math.degrees(lat)
    lon_deg = math.degrees(lon)
    
    return lat_deg, lon_deg



def index(df,i,col):

    a = df.loc[i,col]
    return a



def wgs84_to_itrf2000(X_wgs, Y_wgs, Z_wgs):
    # Translation (meters)
    Tx, Ty, Tz = -0.551 , 0.373 , 0.817

    # Rotations in arc-seconds → radians
    Rx = math.radians(0.001063 / 3600)
    Ry = math.radians(-0.009047 / 3600)
    Rz = math.radians(0.011414 / 3600)

    # Scale factor
    Scale_inv = 1 - 0.004874e-6  # approximate inverse

    # Rotation matrix
    R = np.array([
        [1, Rz, Ry],
        [Rz, 1, Rx],
        [Ry, Rx, 1]
    ])

    # Original coordinates
    XYZ_wgs = np.array([X_wgs, Y_wgs, Z_wgs])

    # Apply transformation
    XYZ_itrf = np.array([Tx, Ty, Tz]) + Scale_inv * R @ XYZ_wgs

    return XYZ_itrf[0], XYZ_itrf[1], XYZ_itrf[2]


# ITRF transformation functions-----------------------------------------------------------------------------------------
def itrf2000_to_itrf2005(X2000, Y2000, Z2000):
    """
    Convert coordinates from ITRF2000 to ITRF2005 using 7-parameter Helmert transformation
    with parameters updated for 5-year difference.
    
    Units:
    - Translations: meters
    - Rotations: arc-seconds
    - Scale: ppm
    """
    y=5 # ref epoch 2000
    # Translation (meters)
    Tx = -0.1e-3 + (0.2e-3 * y)  # meters
    Ty =  0.8e-3 + (-0.1e-3 * y)
    Tz =  5.8e-3 + (1.8e-3 * y)

    # Scale factor (ppm → scale)
    D = -0.40e-9 + (-0.08e-9 * y)  # dimensionless

    # Rotations in arc-seconds → radians
    rx = math.radians(0 / 3600)
    ry = math.radians(0 / 3600)
    rz = math.radians(0 / 3600)

    # Rotation matrix (small-angle approximation)
    R = np.array([
        [1,  rz, -ry],
        [-rz, 1,  rx],
        [ry, -rx, 1]
    ])

    # Original coordinates
    XYZ2000 = np.array([X2000, Y2000, Z2000])

    # Apply Helmert transformation
    XYZ2005 = np.array([Tx, Ty, Tz]) + (1 + D) * (R @ XYZ2000)

    return XYZ2005[0], XYZ2005[1], XYZ2005[2]


def itrf2005_to_itrf2008(X2005, Y2005, Z2005):
    """
    Convert coordinates from ITRF2005 to ITRF2008 using 7-parameter Helmert transformation
    with parameters updated for 8-year difference.
    
    Units:
    - Translations: meters
    - Rotations: arc-seconds
    - Scale: ppm
    """
    y=8 # ref epoch 2000
    # Translation (meters)
    Tx = 2e-3 + (-0.3e-3 * y)  # meters
    Ty =  0.9e-3 + (0 * y)
    Tz =  4.7e-3 + (0 * y)

    # Scale factor (ppm → scale)
    D = -0.94e-9 + (0 * y)  # dimensionless

    # Rotations in arc-seconds → radians
    rx = math.radians(0 / 3600)
    ry = math.radians(0 / 3600)
    rz = math.radians(0 / 3600)

    # Rotation matrix (small-angle approximation)
    R = np.array([
        [1,  rz, -ry],
        [-rz, 1,  rx],
        [ry, -rx, 1]
    ])

    # Original coordinates
    XYZ2005 = np.array([X2005, Y2005, Z2005])

    # Apply Helmert transformation
    XYZ2008 = np.array([Tx, Ty, Tz]) + (1 + D) * (R @ XYZ2005)

    return XYZ2008[0], XYZ2008[1], XYZ2008[2]


def itrf2008_to_itrf2014(X2008, Y2008, Z2008):
    """
    Convert coordinates from ITRF2008 to ITRF2014 using 7-parameter Helmert transformation
    with parameters updated for 4-year difference.
    
    Units:
    - Translations: meters
    - Rotations: arc-seconds
    - Scale: ppm
    """
    y=4 # ref epoch = 2010
    # Translation (meters)
    Tx = -1.6e-3 + (0 * y)  # meters
    Ty =  -1.9e-3 + (0 * y)
    Tz =  -2.4e-3 + (0.1e-3 * y)

    # Scale factor (ppm → scale)
    D = 0.02e-9 + (-0.03e-9* y)  # dimensionless

    # Rotations in arc-seconds → radians
    rx = math.radians(0 / 3600)
    ry = math.radians(0 / 3600)
    rz = math.radians(0 / 3600)

    # Rotation matrix (small-angle approximation)
    R = np.array([
        [1,  rz, -ry],
        [-rz, 1,  rx],
        [ry, -rx, 1]
    ])

    # Original coordinates
    XYZ2008 = np.array([X2008, Y2008, Z2008])

    # Apply Helmert transformation
    XYZ2014 = np.array([Tx, Ty, Tz]) + (1 + D) * (R @ XYZ2008)

    return XYZ2014[0], XYZ2014[1], XYZ2014[2]


def itrf2014_to_itrf2020(X2014, Y2014, Z2014):
    """
    Convert coordinates from ITRF2014 to ITRF2020 using 7-parameter Helmert transformation
    with parameters updated for 5-year difference.
    
    Units:
    - Translations: meters
    - Rotations: arc-seconds
    - Scale: ppm
    """
    dy = 5 # ref epoch 2015 
    # Translation (meters)
    Tx = 1.4e-3 + (0 * dy)  # meters
    Ty =  0.9e-3 + (0.1e-3 * dy)
    Tz =  -1.4e-3 + (-0.2e-3 * dy)

    # Scale factor (ppm → scale)
    D = 0.42e-9 + (0* dy)  # dimensionless

    # Rotations in arc-seconds → radians
    rx = math.radians(0 / 3600)
    ry = math.radians(0 / 3600)
    rz = math.radians(0 / 3600)

    # Rotation matrix (small-angle approximation)
    R = np.array([
        [1,  rz, -ry],
        [-rz, 1,  rx],
        [ry, -rx, 1]
    ])

    # Original coordinates
    XYZ2014 = np.array([X2014, Y2014, Z2014])

    # Apply Helmert transformation
    XYZ2020 = np.array([Tx, Ty, Tz]) + (1 + D) * (R @ XYZ2014)

    return XYZ2020[0], XYZ2020[1], XYZ2020[2]


def itrf2000_to_itrf2005_in2022(X2000, Y2000, Z2000):
    """
    Convert coordinates from ITRF2000 to ITRF2005 using 7-parameter Helmert transformation
    with parameters updated for 22-year difference.
    
    Units:
    - Translations: meters
    - Rotations: arc-seconds
    - Scale: ppm
    """
    y=22 # ref epoch 2000 until 2022
    # Translation (meters)
    Tx = -0.1e-3 + (0.2e-3 * y)  # meters
    Ty =  0.8e-3 + (-0.1e-3 * y)
    Tz =  5.8e-3 + (1.8e-3 * y)

    # Scale factor (ppm → scale)
    D = -0.40e-9 + (-0.08e-9 * y)  # dimensionless

    # Rotations in arc-seconds → radians
    rx = math.radians(0 / 3600)
    ry = math.radians(0 / 3600)
    rz = math.radians(0 / 3600)

    # Rotation matrix (small-angle approximation)
    R = np.array([
        [1,  rz, -ry],
        [-rz, 1,  rx],
        [ry, -rx, 1]
    ])

    # Original coordinates
    XYZ2000 = np.array([X2000, Y2000, Z2000])

    # Apply Helmert transformation
    XYZ2005 = np.array([Tx, Ty, Tz]) + (1 + D) * (R @ XYZ2000)

    return XYZ2005[0], XYZ2005[1], XYZ2005[2]



def itrf2000_to_itrf2020_in2022(X2000, Y2000, Z2000):
    """
    Convert coordinates from ITRF2000 to ITRF2020 using 7-parameter Helmert transformation
    with parameters updated for 5+2 -year difference.
    
    Units:
    - Translations: meters
    - Rotations: arc-seconds
    - Scale: ppb
    """
    dy = 7 # ref epoch 2015 UNTIL 2022 
    # Translation (meters)
    Tx = 0.2e-3 + (-0.1e-3 * dy)  # meters
    Ty =  -0.8e-3 + (0 * dy)
    Tz =  34.2e-3 + (1.7e-3 * dy)

    # Scale factor (ppm → scale)
    D = -2.25e-9 + (-0.11e-9* dy)  # dimensionless

    # Rotations in arc-seconds → radians
    rx = math.radians(0 / 3600)
    ry = math.radians(0 / 3600)
    rz = math.radians(0 / 3600)

    # Rotation matrix (small-angle approximation)
    R = np.array([
        [1,  rz, -ry],
        [-rz, 1,  rx],
        [ry, -rx, 1]
    ])

    # Original coordinates
    XYZ2000 = np.array([X2000, Y2000, Z2000])

    # Apply Helmert transformation
    XYZ2020 = np.array([Tx, Ty, Tz]) + (1 + D) * (R @ XYZ2000)

    return XYZ2020[0], XYZ2020[1], XYZ2020[2]



#-----------------------------------------------------------------------------------------


