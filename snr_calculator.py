import math
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import simple_norm
from photutils import CircularAperture, CircularAnnulus
from photutils import DAOStarFinder
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
#from photutils.aperture import aperture_photometry

#Load the FITS file

#GTC PUPIL
#file_name = 'gtc_CC-circ1_psf_08.70_1.1_0240-COADD.fits'
#file_name = 'gtc_CC-open_psf_08.70_1.1_0240-COADD.fits'

#NEW DRIFT SCAN
#file_name = 'COADD-FAST-2623409.fits'
#file_name = 'COADD-MEDIUM-2624338.fits'
#file_name = 'COADD-SLOW-2624795.fits'
#file_name = 'COADD-SLOW-2625110.fits'

#NEW CHOP NOD
file_name = 'CANARICAM-CN-REDUCED-EXTRACTION.fits'
#file_name = 'CANARICAM-CN-REDUCED-NEW.fits'

#OLD CHOP NOD
#file_name = 'CANARICAM-OLD.fits'

##NOTES: full_box_noise and annulus_standard_deviation
##are actually noise per pixel.

#Open the FITS file and get the data in the hdu
input_file = get_pkg_data_filename(file_name)
hdu = fits.open(input_file)
data = hdu[0].data

#Define the aperture radius and annulus radii
if file_name.__contains__('COADD-FAST'):
    r = 15.39
    #r = 8.66    
else:
    r = 15.39
    #r = 8.66
    
radius = r
    
#Define the subframe size used in coaddition
extraction_size = 40

#Define the subframe size used for photometric extraction of chop-nod
cn_extraction_size = 35

#Define the box size used for statistics
box_length = 10

##Handles the Drift Scan file for SNR calc
if file_name.__contains__('COADD'):
    
    #Define the origins of the boxes
    first_box_x, first_box_y = 0, 0
    second_box_x, second_box_y = extraction_size*2, 0
    third_box_x, third_box_y = 0, extraction_size*2
    fourth_box_x, fourth_box_y = extraction_size*2, extraction_size*2

    first_box = data[first_box_y:first_box_y+box_length, first_box_x:first_box_x+box_length]
    second_box = data[second_box_y:second_box_y+box_length, second_box_x-box_length:second_box_x]
    third_box = data[third_box_y-box_length:third_box_y, third_box_x:third_box_x+box_length]
    fourth_box = data[fourth_box_y-box_length:fourth_box_y, fourth_box_x-box_length:fourth_box_x]

    first_box_mean, first_box_median, first_box_std = round(np.mean(first_box)), round(np.median(first_box)), \
                                                      round(np.std(first_box))

    print('-------------------------------------')
    print('BOX STATISTICS: DRIFT-SCAN')
    print('-------------------------------------')
    print('First Box Mean: ' + str(round(first_box_mean)),'\nFirst Box Median: '+str(round(first_box_median)),
          '\nFirst Box Noise: '+str(round(first_box_std)))
    print('-------------------------------------')

    second_box_mean, second_box_median, second_box_std = round(np.mean(second_box)), round(np.median(second_box)), \
                                                         round(np.std(second_box))
    #print('-------------------------------------')
    print('Second Box Mean: ' + str(round(second_box_mean)),'\nSecond Box Median: '+str(round(second_box_median)),
          '\nSecond Box Noise: '+str(round(second_box_std)))
    print('-------------------------------------')

    third_box_mean, third_box_median, third_box_std = round(np.mean(third_box)), round(np.median(third_box)), \
                                                      round(np.std(third_box))
    #print('-------------------------------------')
    print('Third Box Mean: ' + str(round(third_box_mean)),'\nThird Box Median: '+str(round(third_box_median)),
          '\nThird Box Noise: '+str(round(third_box_std)))
    print('-------------------------------------')

    fourth_box_mean, fourth_box_median, fourth_box_std = round(np.mean(fourth_box)), round(np.median(fourth_box)), \
                                                         round(np.std(fourth_box))
    #print('-------------------------------------')
    print('Fourth Box Mean: ' + str(round(fourth_box_mean)),'\nFourth Box Median: '+str(round(fourth_box_median)),
          '\nFourth Box Noise: '+str(round(fourth_box_std)))
    #print('-------------------------------------')

    full_box = np.array([])

    for j in range(len(first_box)):
        for i in range(len(first_box)):
            full_box = np.append(full_box, first_box[j][i])

    for j in range(len(second_box)):
        for i in range(len(second_box)):
            full_box = np.append(full_box, second_box[j][i])

    for j in range(len(third_box)):
        for i in range(len(third_box)):
            full_box = np.append(full_box, third_box[j][i])

    for j in range(len(fourth_box)):
        for i in range(len(fourth_box)):
            full_box = np.append(full_box, fourth_box[j][i])

    #Define the Centroid for the star/object
    centroid_x = extraction_size
    centroid_y = extraction_size

    #Create the Aperture and Annulus for the Standard Star
    position = [(centroid_x),(centroid_y)]
    aperture = CircularAperture(position, r=radius)

    #Define the values for the Aperture Mask/Data
    #The data portion being the array that corresponds to aperture pixels
    aperture_masks = aperture.to_mask(method='center')
    aperture_data = aperture_masks.multiply(data)

##Handles the Chop Nod file for SNR calc or the sub_frame test file
elif file_name.__contains__('CANARICAM'):

    #Find the centroid
    mean, median, std = sigma_clipped_stats(data, sigma = 50.0)

    for i in range(4,50):
        daofind = DAOStarFinder(fwhm=i, threshold=20.*std)
        sources = daofind(data-median)

        if(sources is not None):
            #print('fwhm: ' + str(i))
            break

    for col in sources.colnames:
        sources[col].info.format = '%.8g'

    #print(sources)

    for i in range(len(sources)):
        if (sources['peak'][i]==np.amax(sources['peak'])):
            centroid_x = round(sources['xcentroid'][0])
            centroid_y = round(sources['ycentroid'][0])

    #DEBUG: Graph to check DAOStarFinder
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=4.)
    norm = ImageNormalize(stretch=SqrtStretch())
    plt.imshow(data, cmap='Greys', origin='lower', norm=norm, interpolation='nearest')
    apertures.plot(color='blue', lw=1.5, alpha=0.5)
    plt.show()
    
    data = data[abs(round(centroid_y)-cn_extraction_size):round(centroid_y)+cn_extraction_size,
           abs(round(centroid_x)-cn_extraction_size):round(centroid_x)+cn_extraction_size]
    
    #Define the origins of the boxes
    first_box_x, first_box_y = 0, 0
    second_box_x, second_box_y = cn_extraction_size*2, 0
    third_box_x, third_box_y = 0, cn_extraction_size*2
    fourth_box_x, fourth_box_y = cn_extraction_size*2, cn_extraction_size*2

    first_box = data[first_box_y:first_box_y+box_length, first_box_x:first_box_x+box_length]
    second_box = data[second_box_y:second_box_y+box_length, second_box_x-box_length:second_box_x]
    third_box = data[third_box_y-box_length:third_box_y, third_box_x:third_box_x+box_length]
    fourth_box = data[fourth_box_y-box_length:fourth_box_y, fourth_box_x-box_length:fourth_box_x]

    first_box_mean, first_box_median, first_box_std = round(np.mean(first_box)), round(np.median(first_box)), \
                                                      round(np.std(first_box))
    
    print('-------------------------------------')
    print('BOX STATISTICS: CHOP-NOD')
    print('-------------------------------------')
    print('First Box Mean: ' + str(round(first_box_mean)),'\nFirst Box Median: '+str(round(first_box_median)),
          '\nFirst Box Noise: '+str(round(first_box_std)))
    print('-------------------------------------')

    second_box_mean, second_box_median, second_box_std = round(np.mean(second_box)), round(np.median(second_box)), \
                                                         round(np.std(second_box))
    #print('-------------------------------------')
    print('Second Box Mean: ' + str(round(second_box_mean)),'\nSecond Box Median: '+str(round(second_box_median)),
          '\nSecond Box Noise: '+str(round(second_box_std)))
    print('-------------------------------------')

    third_box_mean, third_box_median, third_box_std = round(np.mean(third_box)), round(np.median(third_box)), \
                                                      round(np.std(third_box))
    #print('-------------------------------------')
    print('Third Box Mean: ' + str(round(third_box_mean)),'\nThird Box Median: '+str(round(third_box_median)),
          '\nThird Box Noise: '+str(round(third_box_std)))
    print('-------------------------------------')

    fourth_box_mean, fourth_box_median, fourth_box_std = round(np.mean(fourth_box)), round(np.median(fourth_box)), \
                                                         round(np.std(fourth_box))
    #print('-------------------------------------')
    print('Fourth Box Mean: ' + str(round(fourth_box_mean)),'\nFourth Box Median: '+str(round(fourth_box_median)),
          '\nFourth Box Noise: '+str(round(fourth_box_std)))
    #print('-------------------------------------')

    full_box = np.array([])

    for j in range(len(first_box)):
        for i in range(len(first_box)):
            full_box = np.append(full_box, first_box[j][i])

    for j in range(len(second_box)):
        for i in range(len(second_box)):
            full_box = np.append(full_box, second_box[j][i])

    for j in range(len(third_box)):
        for i in range(len(third_box)):
            full_box = np.append(full_box, third_box[j][i])

    for j in range(len(fourth_box)):
        for i in range(len(fourth_box)):
            full_box = np.append(full_box, fourth_box[j][i])

    #Define the centroid
    centroid_x = cn_extraction_size
    centroid_y = cn_extraction_size

    #Create the Aperture and Annulus for the Standard Star
    position = [(centroid_x),(centroid_y)]
    aperture = CircularAperture(position, r=radius)

    #Define the values for the Aperture Mask/Data
    #The data portion being the array that corresponds to aperture pixels
    aperture_masks = aperture.to_mask(method='center')
    aperture_data = aperture_masks.multiply(data)

##    hdu = fits.PrimaryHDU(data)
##    hdu1 = fits.HDUList([hdu])
##    hdu.writeto('CANARICAM-CN-REDUCED-EXTRACTION.fits', overwrite='True')

else:
    print('ERROR: No valid file found!')

def snr_ds_box(aperture_data):    
    #Store the values for the counts of the aperture and annulus
    aperture_counts = 0

    #Store the number of pixels for the aperture and annulus
    aperture_pixels = 0

    #Store the values within the aperture and annulus
    aperture_values = np.array([])

    #Statistics for the box
    full_box_mean, full_box_median, full_box_noise = round(np.mean(full_box)), round(np.median(full_box)), \
                                                     round(np.std(full_box))
##    print(full_box_mean, full_box_median, full_box_noise)
    
    #Add the median to the whole image
    for j in range(len(data)):
        for i in range(len(data[j])):
            #data[j][i] = data[j][i] - annulus_mean
            data[j][i] = data[j][i] - full_box_median            

    #Aperture of whole image with mean
    aperture_data = aperture_masks.multiply(data)
        
    #Get the aperture counts
    for j in range(len(aperture_data)):
        for i in range(len(aperture_data)):
            if aperture_data[j][i]!=0:
                aperture_pixels += 1
                aperture_counts += aperture_data[j][i]
                aperture_values = np.append(aperture_values, aperture_data[j][i])
                
    aperture_counts = round(aperture_counts)
    #Statistics for the Aperture
    aperture_mean = round(np.mean(aperture_values))
    aperture_median = round(np.median(aperture_values))
    aperture_std = round(np.std(aperture_values))
    #signal_to_noise = round(aperture_counts/full_box_noise)
    signal_to_noise = round(aperture_counts/(full_box_noise*math.sqrt(aperture_pixels)))

        #DEBUG Value print outs. FIX
    print('-------------------------------------')
    print('SIGNAL TO NOISE: DRIFT-SCAN')
    print('-------------------------------------')
    print('Centroid X: ' + str(centroid_x))
    print('Centroid Y: ' + str(centroid_y))
    print('-------------------------------------')
    print('Aperture Radius: ' + str(radius))
    print('-------------------------------------')
    print('Box Mean: ' + str(full_box_mean))
    print('Box Median: ' + str(full_box_median))
    print('Noise per pixel: ' + str(full_box_noise))
    print('-------------------------------------')
    print('Aperture Counts: ' + str(aperture_counts))
    print('Aperture Pixels: ' + str(aperture_pixels))
    print('Aperture Mean: ' + str(aperture_mean))
    print('Aperture Median: ' + str(aperture_median))
    print('Aperture Standard Deviation: ' + str(aperture_std))
    print('-------------------------------------')
    print('Signal to Noise: ' + str(signal_to_noise))
    print('-------------------------------------')

    return

def snr_cn_box(aperture_data):
    #Store the values for the counts of the aperture
    aperture_counts = 0

    #Store the number of pixels for the aperture
    aperture_pixels = 0

    #Store the values within the aperture
    aperture_values = np.array([])

    #Statistics for the box
    full_box_mean = round(np.mean(full_box))
    full_box_median = round(np.median(full_box))
    full_box_noise = round(np.std(full_box))
            
    #Add the mean to the whole image with mean
    for j in range(len(data)):
        for i in range(len(data[j])):
            #data[j][i] = data[j][i] - full_box_mean
            data[j][i] = data[j][i] - full_box_median

    #Aperture of whole image with mean    
    aperture_data = aperture_masks.multiply(data)
            
    #Get the aperture counts
    for j in range(len(aperture_data)):
        for i in range(len(aperture_data)):
            if aperture_data[j][i]!=0:
                aperture_pixels += 1
                aperture_counts += aperture_data[j][i]
                aperture_values = np.append(aperture_values, aperture_data[j][i])

    aperture_counts = round(aperture_counts)
    aperture_mean = round(np.mean(aperture_values))
    aperture_median = round(np.median(aperture_values))
    aperture_std = round(np.std(aperture_values))
    #signal_to_noise = round(aperture_counts/full_box_noise)
    signal_to_noise = round(aperture_counts/(full_box_noise*math.sqrt(aperture_pixels)))

    #DEBUG Value print outs. FIX
    print('-------------------------------------')
    print('SIGNAL TO NOISE: CHOP-NOD')
    print('-------------------------------------')
    print('Centroid X: ' + str(centroid_x))
    print('Centroid Y: ' + str(centroid_y))
    print('-------------------------------------')
    print('Aperture Radius: ' + str(radius))
    print('-------------------------------------')
    print('Box Mean: ' + str(full_box_mean))
    print('Box Median: ' + str(full_box_median))
    print('Noise per pixel: ' + str(full_box_noise))
    print('-------------------------------------')
    print('Aperture Counts: ' + str(aperture_counts))
    print('Aperture Pixels: ' + str(aperture_pixels))
    print('Aperture Mean: ' + str(aperture_mean))
    print('Aperture Median: ' + str(aperture_median))
    print('Aperture Standard Deviation: ' + str(aperture_std))
    print('-------------------------------------')
    print('Signal to Noise: ' + str(signal_to_noise))
    print('-------------------------------------')

    return

if file_name.__contains__('COADD'):
    snr_ds_box(aperture_data)
else:
    snr_cn_box(aperture_data)
    
