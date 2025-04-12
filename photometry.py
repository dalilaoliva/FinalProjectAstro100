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
  hdus = [fits.PrimaryHDU()]
  for i in range(1,5):
    master_bias = []
    for fn in bias_files:
      data = fits.getdata(fn, i)
      master_bias.append(data[None])

    master_bias = np.concatenate(master_bias, axis=0)
    master_bias = np.median(master_bias, axis=0)
    print("Created the master bias")
    hdus.append(fits.ImageHDU(master_bias.astype(np.int16), name=f"IM{i}"))  # Store HDU

    hdulist = fits.HDUList(hdus)
    # Write to a new FITS file
    hdulist.writeto('masterBias.fits', overwrite=True)

  test_hdu = fits.open('masterBias.fits')
  test_hdu.info() 


def master_flat_fits(flat_files):
  print("Reading in the flat frames now")
  hdus = [fits.PrimaryHDU()]
  for i in range(1,5):
    master_flat = []
    master_bias = fits.open('masterBias.fits')[i].data
    for fn in flat_files:
      data = fits.getdata(fn)
      data = data - master_bias
      master_flat.append(data[None])

    master_flat = np.concatenate(master_flat, axis=0)
    master_flat = np.median(master_flat, axis=0)
    master_flat = master_flat / np.median(master_flat)
    print("Created the master flat")
    hdus.append(fits.ImageHDU(master_flat.astype(np.int16),name=f"IM{i}"))  # Store HDU
    hdulist = fits.HDUList(hdus)
    # Write to a new FITS file
    hdulist.writeto('masterFlat.fits', overwrite=True)

  test_hdu = fits.open('masterFlat.fits')
  test_hdu.info() 


def data_reduction(info, bias, flats):
  header = info[0].header
  filter_name = header.get("FILTER", "Unknown")
  print(bias[1].header)
  # print(filter_name)

  hdus = [fits.PrimaryHDU()]

  for i in range(1,5):
    filtered_bias = []

    for fn in bias:
       with fits.open(fn) as hdul:
        header = hdul[0].header  # Check primary HDU
        file_filter = header.get("FILTER", None)
        print(file_filter)

    #     if file_filter == filter_name:
    #         data = hdul[i].data
    #         filtered_bias.append(data[None])

    # if filtered_bias:  # Check if list is not empty before concatenating
    #         filtered_bias = np.stack(filtered_bias, axis=0)
    #         hdus.append(fits.ImageHDU(filtered_bias.astype(np.int16), name=f"IM{i}"))

    # if len(hdus) > 1:
    #     hdulist = fits.HDUList(hdus)
    #     hdulist.writeto(f"filteredBias_{filter_name}.fits", overwrite=True)
    # else:
    #     print(f"No matching bias frames found for filter {filter_name}.")

# test_hdu = fits.open('filterdBiasr.fits')
# test_hdu.info() 
  

def histogramCreator(info):
  for i in range(4):
    # data = fits.getdata(info[0],i+1)
    data = info[i+1].data  # Directly access the data from the HDU list
    plt.subplot(2,2, 4-i)
    plt.imshow(data, cmap='afmhot', vmax=np.percentile(data, 95), vmin=np.percentile(data,10))
  plt.show()


telescopeM81part1 = glob.glob("/Users/dalilaoliva/FinalProjectAstro100/KeplerCam/0020.M81.fits")
telescopeM81part2 = glob.glob("/Users/dalilaoliva/FinalProjectAstro100/KeplerCam/0025.M81.fits")

test_hdu = fits.open(telescopeM81part1[0])
test_hdu.info() 

bias_files = glob.glob('/Users/dalilaoliva/FinalProjectAstro100/KeplerCam/bias/*BIAS*.fits')
flat_files = glob.glob('/Users/dalilaoliva/FinalProjectAstro100/KeplerCam/flats/*FLAT*.fits')


master_bias_files = glob.glob('/Users/dalilaoliva/FinalProjectAstro100/KeplerCam/masterBias.fits')
master_flat_files = glob.glob('/Users/dalilaoliva/FinalProjectAstro100/KeplerCam/masterFlat.fits')
# master_bias_fits(bias_files)
# master_flat_fits(flat_files)
# histogramCreator(test_hdu)
# data_reduction(test_hdu, master_bias_files, master_flat_files)

header = master_bias_files[0].header
print(header)
