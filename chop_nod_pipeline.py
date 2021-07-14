import os
import numpy as np
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from matplotlib.colors import LogNorm
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from scipy import fft, ndimage
import matplotlib.pyplot as plt

#CanariCam Image Reduction Pipeline
#Works for ABBA Chop/Nod Cycles
#Author: Amílcar R. Torres-Quijano
#Revision Date: 2021/06/24

#Filename will be the only input for our pipeline
file_name = '0002622599-20200727-CANARICAM-Imaging.fits'
input_file = get_pkg_data_filename(file_name)

#Open the file
hdu_input = fits.open(input_file)

#What is the length of the RAW data?
#This is used to iterate
#print(len(hdu_input))
extension_count = len(hdu_input)-1
#print('extension Count = ' + str(extension_count))

#Get the Number of Frames per Chop within the data
#This is used to iterate.
extension = hdu_input.pop()
frame = extension.shape[0]
chop = extension.shape[1]
#print('Frame = ' + str(frame) + ' \nChop = ' +str(chop))

#Reattach the extension used to get dimensions
hdu_input.append(extension)
#hdu_input.info()

#List to contain each individual fits images
individual_image = []

#Define centroid pos
global pos_x, pos_y, neg_x, neg_y
pos_x = 0
pos_y = 0
neg_x = 0
neg_y = 0

#Define sourceless frame
global sourceless_frame, tramline_frame, noise_frame
sourceless_frame = np.zeros((240, 320),dtype=int)
tramline_frame = np.zeros((240, 320),dtype=complex)
noise_frame = np.zeros((240, 320),dtype=int)

###DEBUG: Number of extensions, Frames per chop, Chops
##print('Extension Count: ' + str(extension_count))
##print('Frames per chop: ' + str(frame))
##print('Chops: ' + str(chop))

##Store each individual image from the MEF file in
##the individual_image array.
for images in range(extension_count):
    extension = hdu_input.pop()
    for j in range(frame):
        for i in range(chop):
            individual_image.append(extension.data[j][i])
            #print(individual_image)
            #print('arr['+str(j)+']['+str(i)+']')            

###DEBUG: Number of individual fits images
##print(len(individual_image))

#Store the sum of each chop-nod pair.
chop_nod_pair = []

#Store each Chop on a Different List
nod1 = 0
nod2 = 1

#PLOTTING FUNCTIONS
def plot_image_space(data, title):
    plt.imshow(data, origin='lower')
    #plt.clim(-500,500)
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

def find_centroid(data, beam):
    mean, median, std = sigma_clipped_stats(data, sigma = 3.0)
    daofind = DAOStarFinder(fwhm=3.0, threshold=3.*std)
    sources = daofind(data-median)

    if (sources is None):
        print('Unable to calculate Centroid...')

    for col in sources.colnames:
        sources[col].info.format = '%.8g'

    for num in range(len(sources)):
        if (sources['peak'][num]==np.amax(sources['peak'])):
            global pos_x, pos_y
            pos_x = round(sources['xcentroid'][num])
            pos_y = round(sources['ycentroid'][num])

    if beam == 'A1':
        global neg_x, neg_y
        neg_x = pos_x
        neg_y = pos_y - 72
    if beam == 'A2':
        #global neg_x, neg_y
        neg_x = pos_x
        neg_y = pos_y - 72
    if beam == 'B1':
        neg_x = pos_x
        neg_y = pos_y - 72
    if beam == 'B2':
        neg_x = pos_x
        neg_y = pos_y - 72

##    print(pos_x, pos_y)
##    print(neg_x, neg_y)

    return

def create_sourceless(data):
    global sourceless_frame
    sourceless_frame = data.copy()
    sourceless_frame[pos_y-20:pos_y+20,
                     pos_x-20:pos_x+20] = 0
    sourceless_frame[neg_y-20:neg_y+20,
                     neg_x-20:neg_x+20] = 0
    #plot_image_space(sourceless_frame, 'SOURCELESS')
    return

def create_noise_frame(data):
    for y in range(240):
        for x in range(16, 320, 16):
            global tramline_frame
            tramline_frame[y,x]=data[y,x]
    #plot_fourier_space(tramline_frame, 'TRAMLINES')
    global noise_frame
    noise_frame = fft.ifft2(tramline_frame.copy()).real
    noise_frame = np.round(noise_frame.copy())
    noise_frame = noise_frame.astype(int)
    #plot_image_space(noise_frame,'NOISE FRAME')
    return

#Run the process for as many extensions available per MEF file
#Each set of ABBA cycle is 4.
iterate = int(extension_count/4)
for j in range(iterate):
    print('FRAME: ' + str(j))
    #Store the frames of each Chop Cycle per Nod.
    beamC1A2 = []
    beamC2A2 = []
    beamC1B2 = []
    beamC2B2 = []
    beamC1B1 = []
    beamC2B1 = []
    beamC1A1 = []
    beamC2A1 = []

    #Store the subtraction of chops from each beam.
    sub_beamA1 = []
    sub_beamA2 = []
    sub_beamB1 = []
    sub_beamB2 = []

    #Store the value of the sum of each the respective beams
    sum_beamA1 = 0
    sum_beamA2 = 0
    sum_beamB1 = 0
    sum_beamB2 = 0

    ##DEBUG: The values of the counting variables used to iterate over the individual fits images
    #print('nod1: ' + str(nod1))
    #print('nod2: ' + str(nod2))
    for i in range(frame):
        beamC1A2.append(individual_image[i+nod1])
        beamC2A2.append(individual_image[i+nod2])
        #print('Nod 1 individual_image['+str(i+nod1)+'] ' + 'Nod 2 individual_image['+str(i+nod2)+']')
        nod1 += 1
        nod2 += 1

    nod1 += frame
    nod2 += frame
    ##print(nod1)
    ##print(nod2)
        
    for i in range(frame):
        beamC1B2.append(individual_image[i+nod1])
        beamC2B2.append(individual_image[i+nod2])
        #print('Nod 1 individual_image['+str(i+nod1)+'] ' + 'Nod 2 individual_image['+str(i+nod2)+']')
        nod1 += 1
        nod2 += 1

    nod1 += frame
    nod2 += frame
    ##print(nod1)
    ##print(nod2)

    for i in range(frame):
        beamC1B1.append(individual_image[i+nod1])
        beamC2B1.append(individual_image[i+nod2])
        #print('Nod 1 individual_image['+str(i+nod1)+'] ' + 'Nod 2 individual_image['+str(i+nod2)+']')
        nod1 += 1
        nod2 += 1

    nod1 += frame
    nod2 += frame
    ##print(nod1)
    ##print(nod2)

    for i in range(frame):
        beamC1A1.append(individual_image[i+nod1])
        beamC2A1.append(individual_image[i+nod2])
        #print('Nod 1 individual_image['+str(i+nod1)+'] ' + 'Nod 2 individual_image['+str(i+nod2)+']')
        nod1 += 1
        nod2 += 1    

    nod1+=frame
    nod2+=frame

    #Subtract the Chops from each beam.
    sd = 1
    for i in range(frame):
        #DATA A1
        dataA1 = beamC1A1[i]-beamC2A1[i]
        find_centroid(dataA1, 'A1')
        create_sourceless(dataA1)
        im_fft = fft.fft2(sourceless_frame)
        create_noise_frame(im_fft)
        clean_frame_A1 = dataA1-noise_frame
##        plot_image_space(dataA1, 'BEAM SUB')
##        plot_image_space(sourceless_frame, 'SOURCELESS')
##        plot_fourier_space(im_fft, 'FFT BEAM SUB')
##        plot_image_space(clean_frame_A1, 'CLEAN FRAME A1')
        
        #DATA A2
        dataA2 = beamC1A2[i]-beamC2A2[i]
        find_centroid(dataA2, 'A2')
        create_sourceless(dataA2)
        im_fft = fft.fft2(sourceless_frame)
        create_noise_frame(im_fft)
        clean_frame_A2 = dataA2-noise_frame
##        plot_image_space(dataA2, 'BEAM SUB')
##        plot_image_space(sourceless_frame, 'SOURCELESS')
##        plot_fourier_space(im_fft, 'FFT BEAM SUB')
##        plot_fourier_space(tramline_frame, 'TRAMLINES')
##        plot_image_space(noise_frame,'NOISE FRAME')
##        plot_image_space(clean_frame_A2, 'CLEAN FRAME A2')


        #DATA B1
        dataB1 = beamC1B1[i]-beamC2B1[i]
        find_centroid(dataB1, 'B1')
        create_sourceless(dataB1)
        im_fft = fft.fft2(sourceless_frame)
        create_noise_frame(im_fft)
        clean_frame_B1 = dataB1-noise_frame
##        plot_image_space(dataB1, 'BEAM SUB')
##        plot_image_space(sourceless_frame, 'SOURCELESS')
##        plot_fourier_space(im_fft, 'FFT BEAM SUB')
##        plot_fourier_space(tramline_frame, 'TRAMLINES')
##        plot_image_space(noise_frame,'NOISE FRAME')
##        plot_image_space(clean_frame_B1, 'CLEAN FRAME B1')
        
        #DATA B2
        dataB2 = beamC1B2[i]-beamC2B2[i]
        find_centroid(dataB2, 'B2')
        create_sourceless(dataB2)
        im_fft = fft.fft2(sourceless_frame)
        create_noise_frame(im_fft)
        clean_frame_B2 = dataB2-noise_frame
##        plot_image_space(dataB2, 'BEAM SUB')
##        plot_image_space(sourceless_frame, 'SOURCELESS')
##        plot_fourier_space(im_fft, 'FFT BEAM SUB')
##        plot_fourier_space(tramline_frame, 'TRAMLINES')
##        plot_image_space(noise_frame,'NOISE FRAME')
##        plot_image_space(clean_frame_B2, 'CLEAN FRAME B2')       

        
        sub_beamA1.append(clean_frame_A1)
        sub_beamA2.append(clean_frame_A2)
        sub_beamB1.append(clean_frame_B1)
        sub_beamB2.append(clean_frame_B2)

    ###DEBUG: The length of arrays regarding the sum of each individual beam.
    ##print(len(sub_beamA1))
    ##print(len(sub_beamA2))
    ##print(len(sub_beamB1))
    ##print(len(sub_beamB2))

    #Sum each of the respective beams
    for i in range(frame):
        sum_beamA1 += sub_beamA1[i]
        sum_beamA2 += sub_beamA2[i]
        sum_beamB1 += sub_beamB1[i]
        sum_beamB2 += sub_beamB2[i]

    reduction = (sum_beamA1+sum_beamA2)-(sum_beamB1+sum_beamB2)
    chop_nod_pair.append(reduction)

#Finalize the image reduction process by summing each chop/nod pair.
final_reduce = 0
for j in range(iterate):
    final_reduce += chop_nod_pair[j]

#print(final_reduce)

#Write out the output file.
hdu = fits.PrimaryHDU(final_reduce)
hdu1 = fits.HDUList([hdu])
hdu.writeto('CANARICAM-CN-REDUCED-NEW.fits', overwrite=True)
