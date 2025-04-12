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

def master_bias_fits(bias_files):
    print("Reading in the bias frames now")
    master_bias = []
    for fn in bias_files:
        data = fits.getdata(fn)  # default is extension 0
        master_bias.append(data[None])

    master_bias = np.concatenate(master_bias, axis=0)
    master_bias = np.median(master_bias, axis=0)
    print("Created the master bias")

    hdu = fits.PrimaryHDU(master_bias.astype(np.int16))
    hdu.writeto('masterBiasSpectroscopy.fits', overwrite=True)

    test_hdu = fits.open('masterBiasSpectroscopy.fits')
    test_hdu.info()

def master_flat_fits(flat_files):
    print("Reading in the flat frames now")
    master_flat = []
    master_bias = fits.getdata('masterBiasSpectroscopy.fits') 

    for fn in flat_files:
        data = fits.getdata(fn)
        data = data - master_bias
        master_flat.append(data[None])

    master_flat = np.concatenate(master_flat, axis=0)
    master_flat = np.median(master_flat, axis=0)
    master_flat = master_flat / np.median(master_flat)
    print("Created the master flat")

    hdu = fits.PrimaryHDU(master_flat.astype(np.int16))
    hdu.writeto('masterFlatSpectroscopy.fits', overwrite=True)

    test_hdu = fits.open('masterFlatSpectroscopy.fits')
    test_hdu.info()

def plotting(bias, flats, science):
    fig = plt.figure(figsize=(12,5))

    plt.subplot(411)
    plt.imshow(bias, aspect='auto', vmax=np.percentile(bias, 90))
    plt.colorbar()
    plt.text(2300, 30, 'Master bias', color='white')
    plt.ylabel('Spatial pixel')

    plt.subplot(412)
    plt.imshow(flats, aspect='auto', vmax=np.percentile(flats, 90))
    plt.colorbar()
    plt.text(2300, 30, 'Master flat', color='white')
    plt.ylabel('Spatial pixel')

    plt.subplot(413)
    data = fits.getdata(science)
    plt.imshow(data, aspect='auto', vmax=np.percentile(data, 90))
    plt.text(2200, 30, 'Science frame (raw)', color='white')
    plt.colorbar()
    plt.ylabel('Spatial pixel')

    plt.subplot(414)
    data = data - bias
    data = data / bias

    plt.imshow(np.log(data), aspect='auto', )
    plt.xlabel('Wavelength pixel')
    plt.ylabel('Spatial pixel')
    plt.text(2300, 30, 'Science frame', color='white')
    plt.colorbar()
    plt.show()

def smaller_frame(science):
    data = fits.getdata(science)
    fig = plt.figure(figsize=(12,5))
    
    plt.subplot(2, 1, 1)
    plt.imshow(np.log(data), aspect='auto', )
    plt.ylim(50, 100)
    plt.xlabel('Wavelength pixel')
    plt.ylabel('Spatial pixel')
    plt.colorbar()

    plt.subplot(2, 1, 2)
    plt.plot(data[80:90].mean(0))
    plt.xlabel('Wavelength pixel')
    plt.ylabel('Flux (arbitrary)')
    plt.ylim(1e2, 1e4)
    # plt.xlim(500, 2700)
    plt.semilogy()

    plt.show()

def pixel_to_wavelength(n_pixels, start_wavelength):
    dispersion = (9000-3700)/2725
    return start_wavelength + dispersion * np.arange(n_pixels)

def plot_1D_spectrum_with_wavelength(science_file, row_range=(75, 85)):
    data = fits.getdata(science_file)
    spectrum_1d = data[row_range[0]:row_range[1]].mean(axis=0)
    
    n_pixels = spectrum_1d.shape[0]
    wavelength = pixel_to_wavelength(n_pixels, 3700)

    plt.figure(figsize=(10, 4))
    plt.plot(wavelength, spectrum_1d, color='mediumblue')
    plt.xlabel('Wavelength (Å)')
    plt.ylabel('Flux (arbitrary units)')
    plt.title('Extracted 1D Spectrum')
    plt.yscale('log')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

bias_files = glob.glob('/Users/dalilaoliva/FinalProjectAstro100/FAST/drive-download-20250325T154020Z-001/*BIAS*.fits')
flat_files = glob.glob('/Users/dalilaoliva/FinalProjectAstro100/FAST/drive-download-20250325T153944Z-001/*FLAT*.fits')
science = '/Users/dalilaoliva/FinalProjectAstro100/FAST/0055.UGC5495.fits'

# master_bias_fits(bias_files)
# master_flat_fits(flat_files)

master_bias = fits.getdata('/Users/dalilaoliva/FinalProjectAstro100/masterBiasSpectroscopy.fits')
master_flat = fits.getdata('/Users/dalilaoliva/FinalProjectAstro100/masterFlatSpectroscopy.fits')

# plotting(master_bias, master_flat, science)
# smaller_frame(science)
# plot_1D_spectrum_with_wavelength(science)