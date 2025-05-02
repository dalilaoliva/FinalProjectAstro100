import numpy as np
import math
import matplotlib.pyplot as plt
import os
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import glob
import warnings
warnings.filterwarnings("ignore")
from matplotlib.gridspec import GridSpec
from astropy.io import fits
from scipy.signal import find_peaks
from numpy.polynomial import Polynomial
from scipy.optimize import curve_fit
from scipy import constants
import pandas as pd
import numpy as np
from scipy.constants import G

data = np.genfromtxt('M104.txt', delimiter=',')

# finding the velocity
def findingRedshift(lambda_observed, lambda_rest):
    return (float(lambda_observed) - float(lambda_rest))/(float(lambda_rest))

def findingVelocity(lambda_observed, lambda_rest): # meters/seconds
    return findingRedshift(lambda_observed, lambda_rest)*constants.c


def listOfRedfshits(data):
    redshift = []
    velocities = []
    for i in data:
        redshift.append(findingRedshift(i, 6563))
        velocities.append(findingVelocity(i, 6563))
    return redshift, velocities

redshifts, velocities = listOfRedfshits(data)
print(redshifts)

print("redshift:", np.average(redshifts))

errors = np.genfromtxt('M104Error.txt', delimiter=',')
_, velocitiesErr = listOfRedfshits(errors)

def errorConstruction(velocitiesErr):
    bars = []
    for i in range(int(len(velocitiesErr)/2)):
        bars.append(velocitiesErr[i+1] - velocitiesErr[i])
    return bars

bars = np.abs(errorConstruction(velocitiesErr))

def plotVelocities(velocities, bars):
    # Generate position values (symmetric slit)
    distance_kpc = np.linspace(-4.6 / 1.085 , 4.6 * 1.085  , len(velocities)) 
    velocities = np.array(velocities) / 1000  # Convert to km/s
    bars = np.array(bars) / (3*1000)

    # Subtract center velocity
    zero_index = np.argmin(np.abs(distance_kpc))
    center_velocity = velocities[zero_index]
    corrected_velocities = velocities - center_velocity

    # Split left and right sides
    left_mask = distance_kpc < 0
    right_mask = distance_kpc >= 0

    # Left side
    x_left = np.abs(distance_kpc[left_mask])
    y_left = np.abs(corrected_velocities[left_mask])
    yerr_left = bars[left_mask]

    # Right side
    x_right = np.abs(distance_kpc[right_mask])
    y_right = np.abs(corrected_velocities[right_mask])
    yerr_right = bars[right_mask]

    # Polynomial fits
    poly_left = Polynomial.fit(x_left, y_left, deg=2).convert()
    poly_right = Polynomial.fit(x_right, y_right, deg=2).convert()

    # Fine x for smooth fits
    x_fit_left = np.linspace(x_left.min(), x_left.max(), 300)
    x_fit_right = np.linspace(x_right.min(), x_right.max(), 300)

    # Plot
    plt.figure(figsize=(10, 6))
    plt.errorbar(x_right, y_right, yerr=yerr_right, fmt='o', color='lightcoral', capsize=4, label='Right Side')
    plt.plot(x_fit_right, poly_right(x_fit_right), '-', color='red', label='Right Fit')
    
    plt.errorbar(x_left, y_left, yerr=yerr_left, fmt='o', color='xkcd:sky blue', capsize=4, label='Left Side')
    plt.plot(x_fit_left, poly_left(x_fit_left), '-', color='blue', label='Left Fit')

    plt.xlabel('Distance from Center (kpc)')
    plt.ylabel('Rotational Speed (km/s)')
    plt.title('Galaxy Rotation Curve: M104')
    plt.grid(True)
    plt.legend()
    plt.show()

    stable_mask_right = x_right > 3
    stable_mask_left = x_left > 3

    # Take average velocity in stable regions
    v_flat_right = np.mean(y_right[stable_mask_right])
    v_flat_left = np.mean(y_left[stable_mask_left])

    # You can also average both sides
    v_flat = (v_flat_right + v_flat_left) / 2
    v_flat = v_flat
    return v_flat

v = plotVelocities(velocities, bars)

def TotalMassCalculation(v, r):
    newr = r * 3.086*10**(19) #in m
    newv = v * 1000 # convert to m/s
    return ((newv)**2 * newr)/(G) 


def convertToSolarMasses(tm):
    return tm/(1.989*10**(30))

# print("Total stellar mass:", convertToSolarMasses(TotalMassCalculation(v, 4.6)))

# Specifically for UGC4972
def irregularGalaxyVel(data, errors):
    sumE = []
    for i in range(len(data)):
        sumE.append(np.abs(errors[i*2] - errors[(i+1)*2-1]))

    anotherV = np.average(sumE)/5890 *constants.c
    return anotherV

# print("Velocity:", irregularGalaxyVel(data, errors))
