import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import simple_norm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import SymLogNorm, LogNorm

#NEW INPUT FILES
file_drift = 'COADD-FAST-2623409.fits'
file_chop = 'CANARICAM-CN-REDUCED-EXTRACTION.fits'
file_chop_old = 'CANARICAM-OLD-EXTRACTION.fits'

def plot_3d_log(file_name, title):
    input_file = get_pkg_data_filename(file_name)
    hdu = fits.open(input_file)
    hdu1 = fits.open(input_file)

    data = hdu[0].data

    ##FOR DRIFT-SCAN
    if file_name.__contains__('COADD'):
        data = data[25:55,25:55]
    
    ##FOR CHOP-NOD
    elif file_name.__contains__('CANARICAM'):
        data = data[20:50,20:50]
    
    ##THE STUFF FOR 3D PLOTS
    x = range(data.shape[1])
    y = range(data.shape[0])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(x,y)
    test = ax.plot_surface(X,Y,data,cmap='magma',norm=SymLogNorm(linthresh=100000, linscale=0.5, vmin=None, vmax=None, clip=False, base=10))

##    #For NOISE
##    test = ax.plot_surface(X,Y,data,cmap='magma',norm=SymLogNorm(linthresh=100000, linscale=0.5, vmin=-5000, vmax=15000, clip=False, base=10))
    
    print(np.amax(data))
    ax.set_title(title)
    
    #zlim max
    ax.set_zlim(0, 45000000)

    #Default Viewing Angle
    ax.view_init(8,-130)

##    #zlim Noise
##    ax.set_zlim(-15000, 9000)
##    
##    #NEW viewing angle
##    ax.view_init(2,-140)
    
    ax.set_xlabel('X Position (px)')
    ax.set_ylabel('Y Position (px)')
    ax.set_zlabel('Counts')
    plt.show()
    return

def plot_3d_log_wire(file_name, title):
    input_file = get_pkg_data_filename(file_name)
    hdu = fits.open(input_file)
    hdu1 = fits.open(input_file)

    data = hdu[0].data

    ##FOR DRIFT-SCAN
    if file_name.__contains__('COADD'):
        data = data[25:55,25:55]
    
    ##FOR CHOP-NOD
    elif file_name.__contains__('CANARICAM'):
        data = data[20:50,20:50]
    
    ##THE STUFF FOR 3D PLOTS
    x = range(data.shape[1])
    y = range(data.shape[0])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(x,y)
    test = ax.plot_wireframe(X,Y,data,cmap='magma',norm=SymLogNorm(linthresh=100000, linscale=0.5, vmin=None, vmax=None, clip=False, base=10))
    print(np.amax(data))
    ax.set_title(title)

    #zlim max
    ax.set_zlim(0, 45000000)

    #Default Viewing Angle
    ax.view_init(8,-130)

##    #zlim Noise
##    ax.set_zlim(-15000, 9000)
##    
##    #NEW viewing angle
##    ax.view_init(2,-140)
    
    ax.set_xlabel('X Position (px)')
    ax.set_ylabel('Y Position (px)')
    ax.set_zlabel('Counts')
    plt.show()
    return

##plot_3d_log(file_drift, 'Drift Scan')
##plot_3d_log(file_chop, 'Chop Nod (With Noise Subtraction)')
##plot_3d_log(file_chop_old, 'Chop Nod (Standard Reduction)')

plot_3d_log_wire(file_drift, 'Drift Scan')
plot_3d_log_wire(file_chop, 'Chop Nod (With Noise Subtraction)')
plot_3d_log_wire(file_chop_old, 'Chop Nod (Standard Reduction)')
