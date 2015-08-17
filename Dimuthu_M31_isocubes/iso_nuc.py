"""
Purpose : Reads ISOCVF cube from the nucleus and extrac a spectrum using the same aperture that we used for that region
          in IRS cubes. This also returns the spectral parameters for both ISO and IRS spectra.
          
Region  : Nucleus of M31

Inputs  : This program reads two fits files.
          cvf76100601.fits - This is a IDL structure based cube of the nucleus.I used it just to get the wavelengths.
          test.fits  - Same fits file but have converted so that it can be read in python.


Outputs : Flux density (MJy / sr) vs Wavelength (mu) plot. # Optional
          A data file which has the final flux and wavelength values.

Notes  : Apertures are rectangular apertures. Therefore the function "squre" takes the pixel coordinates of the four corners
         of the aperture (Which were found in DS9) and generates the parameters of the equations of 4 sides of the aperture.

         The wavelength range overlaps in ISOCVF cubes around 9 mu. This was fixed manually by removing the data that overlaps.


"""

import numpy as np
import matplotlib.pyplot as plt        
import sys
import string, math, glob, os.path                         
import astropy.io.fits as fits
import pylab as py 


def squre (x1,y1,x2,y2,x3,y3,x4,y4):
    """
    This function calculates the slope (m) and the interception (c) of the 4 sides of the aperture.

    Inputs : Pixel coordinates of (x,y) fore each corner of the aperture. These are found using DS9.

    Outputs : Slopes and interceptions of 4 sides of the aperture.

    Note : Used y = mx + c linear equation to get m and c.
    """

    m1 = float(y1-y2) / (x1-x2)
    m2 = float(y3-y2) / (x3-x2)
    m3 = float(y4-y3) / (x4-x3)
    m4 = float(y4-y1) / (x4-x1)

    c1 = y1 - (m1*x1)
    c2 = y2 - (m2*x2)
    c3 = y3 - (m3*x3)
    c4 = y4 - (m4*x4)

    return m1,m2,m3,m4,c1,c2,c3,c4


def getplot(image):
    """
    Takes a 3D array of images and calculates the flux density within a specified aperture.

    Input : 3D array.

    Output : Returns the Flux density within the aperture (Units : MJy/sr) and prints out the fits image showing the
             aperture area in balck.

    """
   
    imageSize = np.shape(image)
    sum_array = [] # Defines an array to collect flux values within the aperture at each wavelength.
    sum1 = 0
    magSize = imageSize
    m = 0
    
    m1,m2,m3,m4,c1,c2,c3,c4 = squre(15,12,11,16,17,22,20,18) #Pixel coordinates of (x,y) fore each corner of the aperture.

    # Looping through all the pixels and adds the pixels into sum_array[] if they are within the aperture.
    for k in range( 0, magSize[0]):

        for i in range( 0, magSize[1]):
               
            for j in range( 0, magSize[2]):
                   # Sets NAN numbers to 0.
                   if np.isnan(image[k,i,j]):
                       image[k,i,j] = 0
                
                   if (j >= (i*m1) + c1) and ( j <= (i*m2) + c2) and (j <= (i*m3) + c3) and ( j >= (i*m4) + c4):
                       value = image[k,i,j]  # If u want to convert units, use *3.046118*(CDELT**2)*(10**2) #unit conversion jy/pix
                       
                       sum1 = sum1 + (value)
                       m = m+1
                     
        sum_array.append(sum1/m)
        m = 0
        sum1 = 0             
    return sum_array

def getwavelength(image) :
    """Gets the wavelengths out from the fits file which has an IDL based structure.
    Input : Fits file that has IDL structure.
    Output : Wavelength range in a 1D array."""

    wl = []
    zaxis = (np.shape(image))[0]
    for i in range(0,zaxis):
        wave = (image[i])[2]
        #print wave,i
        wl.append(wave)

    return wl


def plotvaluesNuc():
    isofits = fits.getdata('cvf76100601.fits') # Reading the IDL based fits file
    wl = getwavelength(isofits) # Gets the wavelength range

    image = fits.getdata('test.fits')  # Reads the python friendly fits file

    flux = getplot(image)

    # Saving the wavelengthand the flux into a file
    with open('isonuctable', 'w') as f:
        for x,y in zip(wl,flux):
            f.write(str(x))
            f.write('\t'+str(y)+'\n')
    f.close()

    # After writing the wavlenght and flux into two columns I manually made two tables out of that removing the overlap at 9 mu.
    iso1 = np.loadtxt("isonuctable1")
    iso2 = np.loadtxt("isonuctable2")

    iso1[0:,1] = iso1[0:,1] + 1.93  # 1st part is shifted up by 1.93

    #Puts Flux into Y and wavelength into X and plots Y vs X and the IRS spectrum in the same figure. 
    Y = []
    X = []

    for (i,aval) in enumerate(iso1[0:,0]):
        Y.append((iso1[0:,1])[i])
        X.append(aval)

    for (i,aval) in enumerate(iso2[0:,0]):
        Y.append((iso2[0:,1])[i])
        X.append(aval)

    #plt.plot(X,Y,'-',linewidth=4.0,label = 'ISOCAM')

    nuc = np.loadtxt("nucFLUX.dat") #IRS data
    irs_flux = nuc[:,1]
    irs_wl = nuc[:,0]

    # Returns ISO and IRS wavelengths and Flux densities.
    return X,Y,irs_wl,irs_flux


x,y,irs_wl,irs_flux = plotvaluesNuc()






