import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import simple_norm
from matplotlib.colors import LogNorm, SymLogNorm

#NEW DRIFT SCAN
#file_name = 'COADD-FAST-2623409.fits'
#file_name = 'ORIG-2623409-FAST-COADD-FULL-26.fits'

#NEW CHOP NOD
file_name = 'CANARICAM-CN-REDUCED-EXTRACTION.fits'
#file_name = 'CANARICAM-CN-REDUCED-NEW.fits'

#OLD CHOP NOD
#file_name = 'CANARICAM-OLD.fits'

output=False

#Open the FITS file and get the data in the hdu
input_file = get_pkg_data_filename(file_name)
hdu = fits.open(input_file)
data = hdu[0].data

def plot(data, title):
    plt.imshow(data,origin='lower')
    plt.colorbar()
    plt.xlabel('X Position (px)')
    plt.ylabel('Y Position (px)')
    plt.title(title)
    plt.show()
    return

def plot_log(data,title):
    plt.imshow(data,norm=SymLogNorm(linthresh=100000, linscale=0.5, vmin=None, vmax=None, clip=False, base=10),origin='lower')
    plt.colorbar()
    plt.xlabel('X Position (px)')
    plt.ylabel('Y Position (px)')
    plt.title(title)
    plt.show()

def threshold(data):
    for j in range(data.shape[0]):
        for i in range(data.shape[1]):
            if data[j][i]<0:
                data[j][i] = 1
    return data

#plot(threshold(data), 'FAST DRIFT - FINAL')
plot_log(data, 'Chop Nod')

if output == True:
    hdu1 = fits.PrimaryHDU(data)
    hdu2 = fits.HDUList([hdu1])
    hdu.writeto('THRESHOLD-SPIE-26F-COADD.fits',overwrite=True)
