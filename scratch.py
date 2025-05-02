# def master_flat_fits(flat_files):
#     print("Reading in the flat frames now")
#     master_flat = []
#     master_bias = fits.getdata('masterBiasSpectroscopy.fits') 

#     for fn in flat_files:
#         data = fits.getdata(fn)
#         data = data - master_bias
#         master_flat.append(data[None])

#     master_flat = np.concatenate(master_flat, axis=0)
#     master_flat = np.median(master_flat, axis=0)
#     master_flat = master_flat / np.median(master_flat)
#     print("Created the master flat")

#     hdu = fits.PrimaryHDU(master_flat)
#     hdu.writeto('masterFlatSpectroscopy.fits', overwrite=True)

#     test_hdu = fits.open('masterFlatSpectroscopy.fits')
#     test_hdu.info()

# def plotting(bias, flats, science):
#     fig = plt.figure(figsize=(12,5))

#     plt.subplot(411)
#     plt.imshow(bias, aspect='auto', vmax=np.percentile(bias, 90))
#     plt.colorbar()
#     plt.text(2300, 30, 'Master bias', color='white')
#     plt.ylabel('Spatial pixel')

#     plt.subplot(412)
#     plt.imshow(flats, aspect='auto', vmax=np.percentile(flats, 90), vmin=np.percentile(flats, 10))
#     plt.colorbar()
#     plt.text(2300, 30, 'Master flat', color='white')
#     plt.ylabel('Spatial pixel')

#     plt.subplot(413)
#     data = fits.getdata(science)
#     plt.imshow(data, aspect='auto', vmax=np.percentile(data, 90))
#     plt.text(2200, 30, 'Science frame (raw)', color='white')
#     plt.colorbar()
#     plt.ylabel('Spatial pixel')

#     plt.subplot(414)
#     data = data - bias
#     data = data / flats

#     plt.imshow(np.log(data), aspect='auto', )
#     plt.xlabel('Wavelength pixel')
#     plt.ylabel('Spatial pixel')
#     plt.text(2300, 30, 'Science frame', color='white')
#     plt.colorbar()

#     # fig = plt.figure()
#     # plt.hist(flats.flatten(), bins=100, log=True)
#     # plt.show()

#     return data

# def subplots(science):
#     fig = plt.figure(figsize=(12,5))
    
#     plt.imshow(science, aspect='auto', vmax=np.percentile(science, 99),
#                vmin=np.percentile(science, 50))
#     plt.xlabel('Wavelength pixel')
#     plt.ylabel('Spatial pixel')
#     plt.colorbar()

#     plt.subplot(2, 1, 2)
#     plt.plot(science[75:85].mean(0))
#     plt.xlabel('Wavelength pixel')
#     plt.ylabel('Flux (arbitrary)')
#     plt.semilogy()
#     plt.show()

# def calibrated_compFile():
#     # Define the directory and specific file
#     directory = 'FAST/'
#     comp_file = '0047.COMP.fits'  # Specify the COMP file you want to process

#     # Load the COMP file (before calibration)
#     comp_file_path = os.path.join(directory, comp_file)
#     with fits.open(comp_file_path) as hdul:
#         comp_data_raw = hdul[0].data

#         print(f"\nViewing FITS file: {comp_file_path}")
#         print("Header:")
#         print(repr(hdul[0].header))
#         print("\nData shape:")
#         print(hdul[0].data.shape)

#     # Trim the raw data (same trimming as for bias and flat)
#     comp_data_trimmed = comp_data_raw[:, :]

#     # Normalize and subtract bias (using combined bias data)
#     bias_subtracted_comp = comp_data_trimmed - combined_bias

#     # Trim both the bias-subtracted COMP data and the normalized flat to match the smaller size
#     min_width = min(bias_subtracted_comp.shape[1], normalized_flat.shape[1])

#     bias_subtracted_comp_trimmed = bias_subtracted_comp[:, :min_width]
#     normalized_flat_trimmed = normalized_flat[:, :min_width]

#     # Apply flat correction (normalize by combined flat)
#     calibrated_comp = bias_subtracted_comp_trimmed / normalized_flat_trimmed

#     # Plot the raw, trimmed, and calibrated COMP image
#     fig, axs = plt.subplots(1, 1, figsize=(12, 4))

#     # Display the calibrated COMP file
#     axs.imshow(calibrated_comp, cmap='gray', origin='lower', aspect='auto')
#     axs.set_title('Calibrated COMP File (After Calibration)')
#     axs.set_xlabel('X Pixel')
#     axs.set_ylabel('Y Pixel')

#     plt.tight_layout()
#     plt.show()

#     # Create a 1D spectrum by summing over the spatial dimension (Y-axis)
#     comp_spectrum_1d = np.sum(calibrated_comp, axis=0)