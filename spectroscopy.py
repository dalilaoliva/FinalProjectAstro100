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

top = 0
bottom = 161

bias_files = glob.glob('/Users/dalilaoliva/FinalProjectAstro100/FAST/drive-download-20250325T154020Z-001/*BIAS*.fits')
flat_files = glob.glob('/Users/dalilaoliva/FinalProjectAstro100/FAST/drive-download-20250325T153944Z-001/*FLAT*.fits')
science = fits.getdata('/Users/dalilaoliva/FinalProjectAstro100/0059.M104.fits')

# createTrimmedMasterBias(bias_files)
# createTrimmedMasterFlats(flat_files)
print(fits.getheader('/Users/dalilaoliva/FinalProjectAstro100/0083.NGC4157.fits'))

master_bias = fits.getdata('/Users/dalilaoliva/FinalProjectAstro100/masterBiasSpectroscopy.fits')
master_flat = fits.getdata('/Users/dalilaoliva/FinalProjectAstro100/masterFlatSpectroscopy.fits')
scienceCOMP = fits.getdata('/Users/dalilaoliva/FinalProjectAstro100/0060.COMP.fits')



# In this section, I trimmed and created the master flat and bias, as well as the normalized flat to find the polynomial fit of this one. 
def createTrimmedMasterBias(bias_files):
    print("Reading in the bias frames now")
    master_bias = []
    for fn in bias_files:
        trimmed_data = fits.getdata(fn)[:, :]  # Trim each bias file
        master_bias.append(trimmed_data)

    combined_bias = np.median(np.array(master_bias), axis=0)
    print("Created the master bias")

    hdu = fits.PrimaryHDU(combined_bias)
    hdu.writeto('masterBiasSpectroscopy.fits', overwrite=True)

    test_hdu = fits.open('masterBiasSpectroscopy.fits')
    test_hdu.info()

def createTrimmedMasterFlats(flat_files):
    print("Reading in the flat frames now")
    master_flat = []
    for fn in flat_files:
        trimmed_data = fits.getdata(fn)[:, :]  # Trim each bias file
        master_flat.append(trimmed_data)

    # Combine all flat data by taking the median across the array
    combined_flat = np.median(np.array(master_flat), axis=0)

    hdu = fits.PrimaryHDU(combined_flat)
    hdu.writeto('masterFlatSpectroscopy.fits', overwrite=True)

    test_hdu = fits.open('masterFlatSpectroscopy.fits')
    test_hdu.info()

def polynomialFitFlats(masterFlat):
    pixels = np.arange(0, masterFlat.shape[-1])

    # Calculate the mean flat for polynomial fitting
    mean_flat = np.mean(masterFlat, axis=0)

    # Iteratively find the best polynomial degree (without weighting)
    min_residual = np.inf
    best_degree = 0

    for degree in range(1, 10):  # Test polynomial degrees from 1 to 9
        p = np.polynomial.Polynomial.fit(pixels, mean_flat, degree)
        poly_fit = p(pixels)
        residual = np.sum((mean_flat - poly_fit) ** 2)  # Sum of squared residuals

        if residual < min_residual:
            min_residual = residual
            best_degree = degree

        # Use the best degree found
        p_best = np.polynomial.Polynomial.fit(pixels, mean_flat, best_degree)
        poly_fit_best = p_best(pixels)

        # Normalize the flat field
        normalized_flat = masterFlat / (mean_flat / poly_fit_best)

    print(p_best)
    hdu = fits.PrimaryHDU(normalized_flat)
    hdu.writeto('masterFlatNormalizedSpectroscopy.fits', overwrite=True)

    test_hdu = fits.open('masterFlatNormalizedSpectroscopy.fits')
    test_hdu.info()

def processScienceImage(science, masterBias, masterFlats):
    # science = science[0].data[:, :]

    # Subtract bias and apply flat-field correction
    bias_subtracted_image = science - masterBias
    flat_corrected_image = bias_subtracted_image / masterFlats

    return flat_corrected_image, science

# Visualization for combined bias and flat frames
def plot_bias_and_flat(combined_bias, combined_flat):
    plt.figure(figsize=(12, 4))

    plt.subplot(1, 2, 1)
    plt.imshow(combined_bias, cmap='gray', origin='lower', aspect='auto')
    plt.colorbar()
    plt.title('Combined Bias Frame')

    plt.subplot(1, 2, 2)
    plt.imshow(combined_flat, cmap='gray', origin='lower', aspect='auto')
    plt.colorbar()
    plt.title('Combined Flat Frame')

    plt.show()

# Function to visualize raw and corrected science images
def plot_science_images(raw_image, corrected_image):
    plt.figure(figsize=(12, 4))

    plt.subplot(1, 2, 1)
    plt.imshow(raw_image, cmap='gray', origin='lower', aspect='auto')
    plt.colorbar()
    plt.ylim(80,100)
    plt.title('Raw Science Image')

    plt.subplot(1, 2, 2)
    plt.imshow(corrected_image, cmap='gray', origin='lower', aspect='auto')
    plt.colorbar()
    plt.ylim(80,100)
    plt.title('Corrected Science Image')

    plt.show()



plot_bias_and_flat(master_bias, master_flat)
polynomialFitFlats(master_flat)

flatCorrectedScienceImage, science = processScienceImage(science, master_bias, master_flat)
master_flat_normalized = fits.getdata('/Users/dalilaoliva/FinalProjectAstro100/masterFlatNormalizedSpectroscopy.fits')

# Comp file parsing and calibration
def COMPfileCalibrated(scienceCOMP, master_bias, master_flat_normalized, ystart, yend):
    comp_data_trimmed = scienceCOMP[ystart:yend, :]

    # Normalize and subtract bias (using combined bias data)
    bias_subtracted_comp = comp_data_trimmed - master_bias[ystart:yend, :]

    # Trim both the bias-subtracted COMP data and the normalized flat to match the smaller size
    min_width = min(bias_subtracted_comp.shape[1], master_flat_normalized.shape[1])

    bias_subtracted_comp_trimmed = bias_subtracted_comp[:, :min_width]
    normalized_flat_trimmed = master_flat_normalized[ystart:yend, :min_width]

    # Apply flat correction (normalize by combined flat)
    calibrated_comp = bias_subtracted_comp_trimmed / normalized_flat_trimmed
    comp_spectrum_1d = np.sum(calibrated_comp, axis=0)
    return calibrated_comp, comp_spectrum_1d

def plot_calibratedCOMPFile(calibrated_comp):
    # Plot the raw, trimmed, and calibrated COMP image
    fig, axs = plt.subplots(1, 1, figsize=(12, 4))

    # Display the calibrated COMP file
    axs.imshow(calibrated_comp, cmap='gray', origin='lower', aspect='auto')
    axs.set_title('Calibrated COMP File')
    axs.set_xlabel('X Pixel')
    axs.set_ylabel('Y Pixel')

    plt.tight_layout()
    plt.show()

calibratedCOMP, calibratedCOMP1D = COMPfileCalibrated(scienceCOMP, master_bias, master_flat_normalized, top , bottom)
plot_calibratedCOMPFile(calibratedCOMP)


# 1D Comparison Lamp Spectrum and Polynomial Fit

def oneDcomparisonLampSpectrum(calibratedCOMP1D):
    peaks, properties = find_peaks(calibratedCOMP1D)

    # Sort peaks by intensity and select the 10 most intense peaks
    peak_intensities = properties['peak_heights'] = calibratedCOMP1D[peaks]
    top_peaks_idx = np.argsort(peak_intensities)[-68:]  # Indices of the 10 highest peaks
    top_peaks = peaks[top_peaks_idx]  # Pixel positions of the 10 highest peaks
    top_peaks_intensity = calibratedCOMP1D[top_peaks]  # Intensities of the 10 highest peaks

    sorted_peak_pixels = np.sort(top_peaks)
    wavelengths = top_peaks_intensity
    sorted_known_wavelengths = np.sort(wavelengths)

    print(len(sorted_peak_pixels))

    # Plot the 1D comparison spectrum with the 10 most intense peaks marked
    plt.figure(figsize=(10, 6))
    plt.plot(calibratedCOMP1D, label='Comparison Spectrum')
    plt.plot(top_peaks, top_peaks_intensity, 'ro', label='Top 10 Intense Peaks')
    plt.xlabel('Pixel Position')
    plt.ylabel('Intensity')
    plt.title('1D Comparison Lamp Spectrum')
    return sorted_peak_pixels, sorted_known_wavelengths

def polynomialFitWavelengthCalibration(calibratedCOMP1D, sorted_peak_pixels, sorted_known_wavelengths):
    # 3. Fit a polynomial (linear or quadratic depending on distortion)
    fit_poly = Polynomial.fit(sorted_peak_pixels, sorted_known_wavelengths, deg=6)  # Use deg=2 if needed

    # 4. Apply the mapping to convert pixels to wavelengths
    pixel_indices = np.arange(len(calibratedCOMP1D))
    calibrated_wavelengths = fit_poly(pixel_indices)

    plt.figure(figsize=(10, 4))

    # Scatter of known data points
    plt.scatter(sorted_peak_pixels, sorted_known_wavelengths, color='blue', label='Matched Peaks')

    # Polynomial fit curve
    x_fit = np.linspace(min(sorted_peak_pixels), max(sorted_peak_pixels), 1000)
    y_fit = fit_poly(x_fit)
    plt.plot(x_fit, y_fit, color='red', label=f'Degree {fit_poly.degree()} Fit')

    plt.xlabel("Pixel Position")
    plt.ylabel("Wavelength (Angstrom)")
    plt.title("Polynomial Fit for Wavelength Calibration")
    plt.grid(True)
    plt.legend()
    plt.show()

    return calibrated_wavelengths

sorted_peak_pixels, sorted_known_wavelengths = oneDcomparisonLampSpectrum(calibratedCOMP1D)
calibrated_wavelengths = polynomialFitWavelengthCalibration(calibratedCOMP1D, sorted_peak_pixels, sorted_known_wavelengths)


# 1D SPectrum with known line positions 
def plot_wavelength_pixel_solution(
    spectrum, pixel_positions, wavelength_values, wavelength_labels=None, annotate=True, 
    figsize=(12, 5), invert_x=False, xlim=None, ylim=None
):

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(spectrum*150, label='1D Spectrum')
    ax.set_xlabel("Pixel Position")
    ax.set_ylabel("Intensity")
    ax.set_title("1D Spectrum with Known Line Positions")
    ax.grid(True)

    for i, pix in enumerate(pixel_positions):
        wl = wavelength_values[i]
        label = f"{wl:.1f} "
        if wavelength_labels is not None:
            label += f" ({wavelength_labels[i]})"
        ax.axvline(pix, color='red', linestyle='--', alpha=0.5)
        if annotate:
            ax.text(pix, spectrum[int(pix)] -2000000, label, 
                    rotation=90, ha='right', va='bottom', fontsize=12, color='red')

    if invert_x:
        ax.invert_xaxis()
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)

    plt.tight_layout()
    plt.show()

pixels_solution = np.array([212, 308, 436, 495, 514,543, 644,755,779,810,823,929,957,985,1007,1243,1272,1340,1425,
                            1499,1787, 1794,1842,1882,1909,1930,1966,1973,2626,2680])
wavelength_solution = np.array([3719.935, 3859.400, 4045.450,4131.830, 4158.590, 4199.860, 4348.064,4510.733, 
                                4545.080, 4589.900, 4609.560,4764.890, 4806.070,4847.820,4879.900, 5227.300,
                                5269.650,5371.489,5496.210, 5606.732,6074.338,6096.1630,6163.594,6266.495,
                               6304.850,6334.428,6382.991 ,6402.246,6965.4300,7032.4127 ])#7383.980, 7514.651])
wavelength_label = np.array(["FeI", "FeI", "FeI" ,"ArII","ArI", "ArI","FeI", "ArII", "ArII", "AII","AII","AII","AII", 
                             "AII","AII","AII", "AII","FeI","FeI","FeI", "AI","NeI","NeI","NeI","NeI","NeI","NeI",
                            "NeI","NeI","AI"])

plot_wavelength_pixel_solution(
    calibratedCOMP1D, 
    pixels_solution, 
    wavelength_solution, 
    wavelength_labels=wavelength_label,
    annotate=False, 
    invert_x=False
)


# Wavelength solution for FAST Spectograph

# Define a polynomial function to fit the data
def polynomial(x, a, b, c, d):
    return  a * x**3 + b * x**2 + c * x + d


def plotWavelengthFASTSpectograph():
    # Fit the polynomial to the data
    popt, pcov = curve_fit(polynomial, pixels_solution, wavelength_solution,)
    a, b, c, d = popt
    # Print out the polynomial coefficients
    print("Polynomial coefficients:", popt)

    # Generate wavelengths using the polynomial solution across the full pixel range
    pixel_range = np.arange(0, max(pixels_solution) +20)
    calibrated_wavelengths = polynomial(pixel_range, *popt)

    # Plot the original data points and the fitted polynomial solution
    plt.figure(figsize=(10, 6))
    plt.plot(pixels_solution, wavelength_solution, 'ro', label='Data Points')
    plt.plot(pixel_range, calibrated_wavelengths, label='Fitted Wavelength Solution')
    plt.xlabel('Pixel Position')
    plt.ylabel('Wavelength (Angstrom)')
    plt.title('Wavelength Solution for FAST Spectrograph')
    plt.legend()
    plt.show()
    return a,b,c,d

coeffs = plotWavelengthFASTSpectograph()


# 1D Wavelength Calibrated Science Spectrum
def apply_wavelength_solution(pixels, coeffs):
    return np.polyval(coeffs, pixels)

corrected_science_image, raw_science_image = processScienceImage(science[top:bottom, :], master_bias[top:bottom, :], master_flat[top:bottom, :])
science_spectrum_1d = np.sum(corrected_science_image, axis=0)
pixels = np.arange(len(calibratedCOMP1D))

wavelength_solution_adjusted = apply_wavelength_solution(pixels, coeffs)
wavelength_solution_adjusted_science = apply_wavelength_solution(pixels, coeffs)

# Define only the H-alpha spectral line
lines = {
    'H-alpha': 6563,
    'H-beta': 4861,
    'H-gamma': 4340,
    'H-delta': 4101,
    # 'He-I': 5876,
    'Na': 5890,
    'O-III': 5007
}

def plot_1d_spectrum_with_lines(wavelengths, spectrum_1d):
    plt.figure(figsize=(10, 6))

    # Plot the 1D spectrum
    plt.plot(wavelengths, spectrum_1d, label='Corrected Lamp Spectrum')
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('Intensity')
    plt.title(f'1D Wavelength-Calibrated Science Spectrum {top}-{bottom}')
    plt.grid(True)

    # Plot the vertical line for the H-alpha line
    for line_name, line_wavelength in lines.items():
        plt.axvline(x=line_wavelength, color='blue', linestyle='--', alpha=0.7)
        plt.text(line_wavelength +5, np.max(spectrum_1d) * 0.5, line_name, rotation=90, verticalalignment='bottom', color='red')

    plt.legend()
    plt.show()

# Plot the 1D spectrum with only the H-alpha line using the wavelength solution
plot_1d_spectrum_with_lines(wavelength_solution_adjusted, calibratedCOMP1D)

def plot_1d_spectrum_with_lines(wavelengths, spectrum_1d):
    plt.figure(figsize=(10, 6))

    # Plot the 1D spectrum
    plt.plot(wavelengths, spectrum_1d, label='Corrected Science Spectrum')
    plt.xlabel(' Observed Wavelength (Angstroms)')
    plt.ylabel('Intensity')
    # plt.xlim(5800, 6000)
    # plt.ylim(-0.2, 1.2)
    plt.title(f'1D Wavelength-Calibrated Science Spectrum M 104')
    plt.grid(True)

    # Plot the vertical line for the H-alpha line
    for line_name, line_wavelength in lines.items():
        plt.axvline(x=line_wavelength, color='blue', linestyle='--', alpha=0.7)
        plt.text(line_wavelength +5, np.max(spectrum_1d) * 0.5, line_name, rotation=90, verticalalignment='bottom', color='red')

    plt.legend()
    plt.show()

plot_1d_spectrum_with_lines(wavelength_solution_adjusted_science, science_spectrum_1d)

def plot_1d_spectrum_with_lines_redshift(wavelengths, spectrum_1d):
    # z = 0.0028597821118390996  # NCG4157
    # z = 0.011333828522920208    # UGC4972
    z = 0.002412514602062056  # M104
    wavelength_rest = wavelengths / (1 + z)
    plt.figure(figsize=(10, 6))

    # Plot the 1D spectrum
    plt.plot(wavelength_rest, spectrum_1d, label='Corrected Science Spectrum')
    plt.xlabel(' Observed Wavelength (Angstroms)')
    plt.ylabel('Intensity')
    # plt.xlim(6500, 6700)
    # plt.ylim(-0.2, 1.5)
    plt.title(f'1D Wavelength and Redshift-Calibrated Science Spectrum M 104')
    plt.grid(True)

    # Plot the vertical line for the H-alpha line
    for line_name, line_wavelength in lines.items():
        plt.axvline(x=line_wavelength, color='blue', linestyle='--', alpha=0.7)
        plt.text(line_wavelength +5, np.max(spectrum_1d) * 0.5, line_name, rotation=90, verticalalignment='bottom', color='red')

    plt.legend()
    plt.show()

plot_1d_spectrum_with_lines_redshift(wavelength_solution_adjusted_science, science_spectrum_1d)


poly = [-1.68082138e-07,  6.07623032e-04, 8.74859881e-01, 3.55690966e+03]
g = np.poly1d(poly)
pixels = np.arange(len(science[0]))
wavelengths = g(pixels)

ind_halpha = np.where(np.abs(wavelengths - 6563*1.0025) < 25)[0]

halpha_intensity = science[:, ind_halpha].mean(-1)

fig = plt.figure()
plt.plot(halpha_intensity)
plt.show()

# Collapse the 2D image into a 1D spectrum (sum over spatial axis, usually vertical)
spectrum_1d = np.sum(science[0:10], axis=0)  # or axis=1 depending on orientation
spectrum_75 = np.sum(science[70:90], axis=0) 
spectrum_140 = np.sum(science[150:161], axis=0)

fig = plt.figure()
# Plot intensity vs. pixel
plt.plot(spectrum_1d, label='Row 0–10')
plt.plot(spectrum_75, label='Row 75–85')
plt.plot(spectrum_140, label='Row 150–161')
plt.xlabel('Pixel (or Wavelength)')
plt.ylabel('Intensity (Counts)')
plt.title('1D Extracted Spectrum Zoomed in H-alpha')
# plt.xlim(1970,2040)
plt.grid(True)
plt.legend()
plt.show()

fig = plt.figure(figsize=(10,4))
plt.imshow(science, aspect='auto', vmax=np.percentile(science, 95),
       vmin=np.percentile(science, 10), extent=[wavelengths[0], wavelengths[-1], 0, 1])

plt.show()