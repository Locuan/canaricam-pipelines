import os
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from matplotlib.colors import LogNorm, NoNorm, Normalize, PowerNorm
from scipy import fft, ndimage
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
import statistics
import numpy as np

#Input Directory
input_path = 'C:\\Users\\Locuan\\Documents\\GitHub\\gtc-drift-scan\\Drift Scan - Raw Data'
output_path = 'C:\\Users\\Locuan\\Documents\\GitHub\\gtc-drift-scan\\'

#Input File
# Fast Drift
file_name = input_path+'\\0002623409-20200727-CANARICAM-Imaging.fits'
drift_type = 'FAST-2623409'

input_file = get_pkg_data_filename(file_name)

#Open the file
hdu_input = fits.open(input_file)

#Get the exentions and the number of images in the MEF FITS file
extension = hdu_input.pop()
frame = extension.shape[0]

dim_y = extension.shape[2]
dim_x = extension.shape[3]

#Reattach the extension used to get dimensions
hdu_input.append(extension)

###DEBUG: Base HDU information
##hdu_input.info()

#List to contain each individual fits images
individual_image = []

for images in range(frame):
    individual_image.append(extension.data[images][0])

# DEFINE IF FINAL COADDED FRAME IS TO BE OUTPUTTED AS FITS FILE
output = True

# DEFINE IF TO DEBUG THE PROGRAM
debug = False

# DEFINE INITIAL INDEXES FOR ROLLING MEDIAN
# FOR THE DIFFERENT DRIFT FILES, THE ROLLING MEDIAN VALUES WILL CHANGE
global i, j, k, l, m

# ROLLING MEDIAN VALUES FAST DRIFT
i,j,k,l,m = 1, 27, 53, 79, 105

#DEFINE THE FIVE DITHER FRAMES
global frame_one, frame_two, frame_three, frame_four, frame_five
frame_one = individual_image[i]
frame_two = individual_image[j]
frame_three = individual_image[k]
frame_four = individual_image[l]
frame_five = individual_image[m]

#DEFINE THE ARRAYS TO STORE ROLLING MEDIAN STATS
global noise, mean, median
noise = np.array([])
mean = np.array([])
median = np.array([])

# SET THE EXTRACTION SIZE FOR THE PHOTOMETRIC EXTRACTION
extraction_size = 40

# SET THE SOURCE REMOVAL BOX SIZE (size = value*2)
removal_size = 14

# DEFINE TRAMLINE LOCATION (TRAMLINES APPEAR AT EVERY STEP OF THE VALUE)
tramline_location = 16

# DEFINE THE FULL ITERATION OVER THE FRAMES
#If I start this at a value other than prior to 53 it won't update the rolling background.
#This is due to the indexes only updating if it meets 53.
start = 1
end = 2901

#DEFINE THE COADDITION ARRAY
coadd = np.zeros((extraction_size*2,extraction_size*2),dtype=int)

#Create the Median Frame
def create_median_frame(frame_one, frame_two, frame_three, frame_four):
    median_frame = np.zeros((dim_y,dim_x),dtype=int)
    median_test = []
    for j in range(dim_y):
        for i in range(dim_x):
            median_test.append(frame_one[j][i])
            median_test.append(frame_two[j][i])
            median_test.append(frame_three[j][i])
            median_test.append(frame_four[j][i])
            median_test.append(frame_five[j][i])
            median_frame[j][i]=statistics.median(median_test)
            median_test = []
    return median_frame

#PLOTTING FUNCTIONS
def plot_image_space(data, title):
    plt.imshow(data, origin='lower')
    if title.__contains__('-'):
        plt.clim(-500,500)
    plt.colorbar()
    plt.xlabel('X Position (px)')
    plt.ylabel('Y Position (px)')
    plt.title(title)
    plt.show()
    return

def plot_fourier_space(data, title):
    plt.imshow(np.abs(data), norm=LogNorm(vmin=5), origin='lower')
    plt.colorbar()
    plt.xlabel('X Position (px)')
    plt.ylabel('Y Position (px)')
    plt.title(title)
    plt.show()
    return

#CREATE THE ROLLING MEDIAN
for iterate in range(start, end):
    print('FRAME: ' + str(iterate))
    
    median_frame = create_median_frame(frame_one, frame_two, frame_three, frame_four)
    sub_frame = individual_image[iterate]-median_frame
    print('Dither Frames: ' + str(i) + ' ' + str(j) + ' ' + str(k) + ' ' +str(l) + ' ' + str(m))
##    print('Median Frame Statistics: \nStandard Deviation: ' + str(round(np.std(median_frame))) +
##          '\nMean: ' + str(round(np.mean(median_frame))) + '\nMedian: ' + str(round(np.median(median_frame))))

##    #UPDATE THE MEDIAN FRAME STATISTICS ARRAYS
##    noise = np.append(noise, np.std(median_frame))
##    mean = np.append(mean, np.mean(median_frame))
##    median = np.append(median, np.median(median_frame))
    
    # FIND CENTROID OF SUBFRAME
    cen_mean, cen_median, cen_std = sigma_clipped_stats(sub_frame, sigma = 3.0)
    daofind = DAOStarFinder(fwhm=2.0, threshold=3.*cen_std)
    sources = daofind(sub_frame - cen_median)

    # DEBUG: PLOT ORIGINAL FRAME, MEDIAN FRAME, AND SUBFRAME
    if debug == True:
        plot_image_space(individual_image[iterate], 'FRAME: ' + str(iterate))
        plot_image_space(median_frame, 'MEDIAN FRAME')
        plot_image_space(sub_frame, 'FRAME ' + str(iterate) + ' - MEDIAN FRAME')

    if (sources is None):
        print('Unable to calculate Centroid...')     
        
    for col in sources.colnames:
        sources[col].info.format = '%.8g'
    #print(sources)
    for num in range(len(sources)):
        if (sources['peak'][num]==np.amax(sources['peak'])):
            centroid_x = round(sources['xcentroid'][num])
            centroid_y = round(sources['ycentroid'][num])

    # If Centroid Calculation is innacurate discard data.
    # Or if Centroid position does not allow extraction size discard data.
    if (centroid_x >= dim_x or centroid_x <= 0) or (centroid_y >= dim_y or centroid_y <= 0):
        print('Unable to calculate Centroid...')

    # CREATE SOURCELESS FRAME
    sourceless_frame = sub_frame.copy()
    sourceless_frame[centroid_y-removal_size:centroid_y+removal_size,
                     centroid_x-removal_size:centroid_x+removal_size] = 0

    # DEBUG: PLOT THE SOURCELESS SUBFRAME
    if debug == True:
        plot_image_space(sourceless_frame, 'SOURCELESS FRAME ' + str(iterate) + ' - MEDIAN FRAME')

    # CREATE FFT FROM SOURCELESS FRAME
    im_fft = fft.fft2(sourceless_frame.copy())

    # DEBUG: PLOT FOURIER OF SOURCELESS SUBFRAME
    if debug == True:
        plot_fourier_space(im_fft, 'FFT SOURCELESS FRAME ' + str(iterate) + ' - MEDIAN FRAME')

    # SELECT THE TRAM LINES
    tramlines = np.zeros((dim_y,dim_x),dtype=complex)
    for y in range(dim_y):
        for x in range(tramline_location, dim_x, tramline_location):
            tramlines[y,x]=im_fft[y,x]

    # DEBUG: PLOT THE TRAMLINES
    if debug == True:
        plot_fourier_space(tramlines, 'TRAM LINES FRAME ' + str(iterate))

    # CREATE THE NOISE FRAME
    noise_frame = fft.ifft2(tramlines.copy()).real
    noise_frame = np.round(noise_frame.copy())
    noise_frame = noise_frame.astype(int)

    # DEBUG: PLOT THE NOISE FRAME
    if debug == True:
        plot_image_space(noise_frame, 'NOISE FRAME ' + str(iterate))

    # SUBTRACT NOISE FRAME FROM SUB_FRAME
    clean_frame = sub_frame - noise_frame

    # DEBUG: PLOT THE CLEAN FRAME
    if debug == True:
        plot_image_space(clean_frame, 'FRAME ' + str(iterate) + ' - MEDIAN FRAME - NOISE FRAME')

    # CREATE THE PHOTOMETRIC EXTRACTION
    phot_extraction = clean_frame[(centroid_y-extraction_size):(centroid_y+extraction_size),
                                  (centroid_x-extraction_size):(centroid_x+extraction_size)]

    # COADD THE FRAMES
    coadd+=phot_extraction

    #UPDATE THE ROLLING MEDIAN
    if iterate == k and m != 2900:
        i+=1
        j+=1
        k+=1
        l+=1
        m+=1
        frame_one = individual_image[i]
        frame_two = individual_image[j]
        frame_three = individual_image[k]
        frame_four = individual_image[l]
        frame_five = individual_image[m]

# DEBUG: PLOT THE COADDED FRAME
if debug == True:
    plot_image_space(coadd, 'COADDED FRAME') 

if output == True:
    # OUTPUT THE COADDED FRAME
    hdu = fits.PrimaryHDU(coadd)
    hdu1 = fits.HDUList([hdu])
    hdu.writeto(os.path.join(output_path,'COADD-'+drift_type+'.fits'),overwrite=True)
    print('Outputted .fits to default folder location...\nFinished Program Run...')
else:
    print('Finished Program Run...')

###CREATE THE FILES FOR THE MEDIAN FRAME STATISTICS
##with open(output_path+'median_frame_noise.txt','w') as filehandle:
##    for listitem in noise:
##        filehandle.write('%s\n' % listitem)
##
##with open(output_path+'median_frame_mean.txt','w') as filehandle:
##    for listitem in mean:
##        filehandle.write('%s\n' % listitem)
##
##with open(output_path+'median_frame_median.txt','w') as filehandle:
##    for listitem in median:
##        filehandle.write('%s\n' % listitem)
