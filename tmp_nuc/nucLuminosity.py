import numpy as np                     
#import matplotlib.pyplot as plt
import math

""" This code takes the the flux from the nucleus of M31 within the IRAS filter
and subtracts the starlight and calculates the 12 micron luminosity and from that it calculates
the bolometric luminosity of the nucleus"""


nuc= np.loadtxt("nucFLUX")  # this file has flux only between 8.5 to 15 microns.
CDELT = 0.000513888895512  # Deg/pix for SL1 #PB: 1.84"/pix
pix_area = 17*27 # Number of pixels
C = 2.999*10**8  # Spped of light 


y = nuc[:,1]  # Flux in MJy/sr/pix
x = nuc[:,0]  # wavelength in microns


# Calculating starlight using the balckbody function given in PAHFIT
lam = x   
B = 1.4387752*10**4/ (lam*5000) #PB: h*c/(lambda*kB*T) [unitless] with lambda in microns
A = 3.97289*10**13/((lam**3)*(np.exp(B) - 1)) # Planck B_nu [MJy/sr] with lambda in microns
starlight = A*1.05*10**(-10) # PB: what is this number? Must be some kind of scaling factor for starlight (PAHFIT output?)

#y = y - starlight  # subtracting starlight
x = x*10**(-6)  # converting microns into meeters 
y = y*C/(x**2)  # converting f_nu [MJy/sr] into f_lambda [MJy/sr*m*s] #PB flam = (nu/lambda)*fnu = (c/lambda^2)*fnu


newY = y*pix_area*3.046118*(CDELT**2)*(10**2)  # Converting flambda into Jy #PB: mult by solid angle & factor 1e6 
#                                              # solid angle in sr: pix_area*3.046118*(CDELT**2)*(1e-4)
#                                              # conversion MJy->Jy multiply by 1e6
                                               # so this number is in Jy/m*s = 1e-26 W/m^3 
y = newY*10**(-26)  # Converting y into Jy continues #PB: no, this converts y into W/m^3


I = np.trapz(y, x)  # Integrate the spectrum #PB inetgral of F_lam [W/m^3]*dlam[m] is in W/m^2, units of flux


flux = I

D = 780*10**3*(3.08567758*10**16)  # 780 Kpc into meters

L = flux*4*math.pi*D**2  # This is the 12 micron luminosity

# The bolometric luminosity is just 5 times this # PB: factor 5 is from Spinoglio & Malkan 1989, ApJ 342, 83. 
# PB: better value is 1/0.07 = 14.3 (Spinoglio et al 1995)

BL = L*5*10**7   # from W to ergs by multiplying 10*7



print BL




