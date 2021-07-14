import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from photutils import aperture_photometry, CircularAperture
from astropy.visualization import simple_norm
from photutils import CircularAperture
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import sigma_clipped_stats
from photutils import detect_threshold

#Load the FITS file

#NEW_FILES
file_name = 'COADD-FAST-2623409.fits'
#file_name = 'CANARICAM-CN-REDUCED-EXTRACTION.fits'
#file_name = 'CANARICAM-OLD-EXTRACTION.fits'

#Open the FITS file and get the data in the hdu
input_file = get_pkg_data_filename(file_name)
hdu_input = fits.open(input_file)
hdu1_input = fits.open(input_file)
data = hdu_input[0].data
threshold = data.copy()

threshold_value = -3000

write_out = False

if file_name.__contains__('COADD'):
    centroid_x = 40
    centroid_y = 40
    fwhm = 2.80
    radius = 40
    end = 39
else:
    centroid_x = 35
    centroid_y = 35
    fwhm = 3.32
    radius = 35
    end = 34

#Find the centroid
positions = centroid_x, centroid_y
apertures = CircularAperture(positions, r=4.)
norm = ImageNormalize(stretch=SqrtStretch())
##plt.imshow(data, origin='lower', norm=norm, interpolation='nearest')
##apertures.plot(color='red', lw=1.5, alpha=0.5)
##plt.show()

def plot_curve_of_growth(centroid_x, centroid_y):
    aperture_sum = np.array([])
    radii = np.array([])
    
    #Threshold because of cross talk (Threshold to the average (1/3 of that?) of the whole image?) SD, MEAN IN SKY. MEAN - 3*SD Threshold to zero there. (SEND FITS TO CHRIS)
    #Look at the standard deviation (3 times the SD - the average?) threshold below that.
    for j in range(len(threshold)):
        for i in range(len(threshold)):
            if threshold[j][i]<threshold_value:
                threshold[j][i] = 0

    positions = centroid_x, centroid_y
    apertures = CircularAperture(positions, r=4.)
    norm = ImageNormalize(stretch=SqrtStretch())
##    plt.imshow(threshold, origin='lower', norm=norm, interpolation='nearest')
##    apertures.plot(color='red', lw=1.5, alpha=0.5)
##    plt.show()

    for i in range(1, end+2):

        #Determine Aperture Size
        radius = i

        if radius >= end+2:
            break

        #Create the Aperture for the Standard Star
        position = [(centroid_x), (centroid_y)]
        aperture = CircularAperture(position, r=radius)

        phot_table = aperture_photometry(threshold, aperture)
        aperture_sum = np.append(aperture_sum, phot_table['aperture_sum'].data)
        radii = np.append(radii, i)

    plt.title("Curve of Growth")
    plt.xlabel("Radius (pixels)")
    plt.ylabel("Aperture Sum (counts)")
    plt.plot(radii, aperture_sum, color="blue")
    plt.show()

    for i in range(len(aperture_sum)):
        if aperture_sum[i] >= round(aperture_sum[-1]*0.95):
            print('Total Flux: ' + str(round(aperture_sum[-1])))
            print('95% of Total Flux: ' + str(round(aperture_sum[-1]*0.95)))
            break

    print("Aperture Radius = " + str(i))
    #print(aperture_sum)

    return aperture_sum

def aperture_sum(centroid_x, centroid_y, radius):
    for j in range(len(threshold)):
        for i in range(len(threshold)):
            if threshold[j][i]<threshold_value:
                threshold[j][i] = 0
    position = [(centroid_x), (centroid_y)]
    aperture = CircularAperture(positions,r=radius)
    phot_table = aperture_photometry(threshold, aperture)
    aperture_sum = phot_table['aperture_sum'].data

    test = aperture_sum*0.95
    print(aperture_sum,test)

    r = 1.00
    
    while True:
        position = [(centroid_x), (centroid_y)]
        aperture = CircularAperture(position, r)
        phot_table = aperture_photometry(threshold, aperture)
        aperture_sum = phot_table['aperture_sum'].data
        #print(aperture_sum)
        if aperture_sum <= test:
            r += 0.01
            #print(round(r, 2))
        else:
            break
    print(round(r, 2))
    
    return

#plot_curve_of_growth(centroid_x, centroid_y)
aperture_sum(centroid_x, centroid_y, radius)
