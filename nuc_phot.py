from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, Column
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area # need astropy 1.0+ for this
from photutils import SkyRectangularAperture, aperture_photometry
import numpy as np
import glob
import matplotlib.pyplot as plt
from astropy.analytic_functions.blackbody import blackbody_nu

# center of SL nucleus map
ap_ctr = SkyCoord(10.684029,41.269439,frame='icrs', unit='deg')
# SL extraction aperture is 50" wide, 30" high, with PA of 45 deg
phot_ap = SkyRectangularAperture(ap_ctr,w=50*u.arcsec,h=30*u.arcsec,theta=45*u.degree)

def dophot(img_list, aperture = phot_ap):
    ''' Do aperture photometry on a list of images
    input: img_list: list of images
           phot_ap: aperture to use (same for all)

    output: astropy table with photomety for all
    '''
    # make an empty table to hold the output
    phot_tab = Table(names=('img_name','xcenter', 'ycenter', 'aperture_sum','MJy_counts'), dtype=('S30','f8', 'f8', 'f8','f8'))

    # loop over the images
    for (img_num,img_name) in enumerate(img_list):
        phot_tab.add_row(('',0.0,0.0,0.0,0.0))
        img = fits.open(img_name)
        out_tab = aperture_photometry(img, aperture)
#        print out_tab

        # store name of image
        phot_tab['img_name'][img_num] = img_name

        # copy output into final table
        for col in ['xcenter', 'ycenter', 'aperture_sum']:
            phot_tab[col][img_num] = out_tab[col][0]

        # calibrate photmetry 
        obs_val = out_tab['aperture_sum'][0] * out_tab['aperture_sum'].unit
        phot_tab['MJy_counts'][img_num] = calib_phot(obs_val, img, output_units='MJy')

    # done loop over images
    return(phot_tab)

# TODO: include surfce-brightness-to-flux conversion
def calib_phot(input_value, img, output_units='MJy'):
    '''
    Convert the aperture_sum value to output_units
    input: input_value: value to be converted, with units
           img: image HDU, for ancillary/calibration data
           output_units: units to convert output value to

    output: calibrated value
    '''
    #  do we already have unit information?
    if input_value.unit.is_unity(): # means unit isn't given in table

        # so try to figure out what to do from image header
        hdr = img[0].header
        if 'BUNIT' in hdr:
            # shouldn't get here if coming from photutils, but anyway..
            obs_val = input_value * u.Unit(hdr['BUNIT']) # might fail if BUNIT badly formatted..
        elif 'MAGZP' and 'VEGAFLUX' in hdr: # convert to Jansky
            mag = hdr['MAGZP'] - 2.5*np.log10(input_value)
            obs_val = hdr['VEGAFLUX']*10.0**(-0.4*mag) * u.Jy
        else:
            print 'Not enough info to calibrate'
    # do we need a surface-brightness to flux conversion?
    elif input_value.unit == u.MJy/u.sr: #        (this is not perfectly general but oh well)
        hdr = img[0].header
        wcs = WCS(hdr)
        # proj_plane_pixel_area returns values in same units as CDELT,etc: likely deg^2
        pxarea = (proj_plane_pixel_area(wcs) * (u.degree**2)).to(u.arcsec**2) 
        intermed = input_value.to(u.Jy/u.arcsec**2) # not strictly necessary but easier to follow
        obs_val = intermed * pxarea
    else:
        obs_val = input_value

    #  now do the conversion
    try:
        calib_val = obs_val.to(output_units).value
    except UnitsError:
        print 'Problem with unit conversion'
        return(None)

    return(calib_val)

def makeplot(photdat=None):

    if photdat == None:
        # dophotometry photometry
        imglist = glob.glob('m31nuc_f1*.fits')+glob.glob('m31_2mass_*.fits')+glob.glob('m31_?_bgsub_bc_nuc.fits')    
    #    imglist[0], imglist[1] = imglist[1],imglist[0] # need this swap because H is before J in alpha but not wavelength order!
        photdat = dophot(imglist)
    photwaves = np.array(([1.1, 1.6, 1.6,1.2,2.2,3.6,4.5,5.8,8]))
    photcorr = np.array([1.0,1.0,1.0,1.0,1.0,0.91,0.94,0.68,0.74]) # IRAC extd src correction
    photvals = photdat['MJy_counts']* photcorr

#    load the IRS spectrum and convert to MJy
    nuc_wave,nuc_irs = np.loadtxt('../pb_m31_spectra/nucFLUX',unpack=True)
    nuc_irs = nuc_irs*((1500*u.arcsec**2).to(u.sr).value)

#   read the nu Pav spectrum
    nupav_wave, nupav_flux = np.loadtxt('nu_pav_spect.txt',unpack=True,usecols=[0,1])
    # normalize
    find_8micron = np.searchsorted(nupav_wave,8)
    nupav_flux = nupav_flux*(photvals[-1]/nupav_flux[find_8micron])

#   create a RJ tail to compare to
    bb_wl = np.arange(1.0,22,0.4)
    bb = blackbody_nu(bb_wl*u.micron,5000)
    # normalize it to the IRAC flux at 8um
    find_8micron = np.searchsorted(bb_wl,8)
    bb = bb*(photvals[-1]/bb[find_8micron].value)

    # plot
    f,ax=plt.subplots()
    ax.plot(nuc_wave,nuc_irs*1e6,ls='solid',marker=None, lw=2,label= 'M31 IRS')
    ax.plot(bb_wl, bb.value*1e6, ls='dashed', color='k',marker=None, lw=2, label = '5000K BB' )
    ax.plot(nupav_wave, nupav_flux*1e6, ls = 'solid', color='k',marker=None, lw=2,label='nu Pav')
    ax.plot(photwaves[0:2],photvals[0:2]*1e6,ms=10,label='HST') # factor 1e6 makes plot in Jy, gives nice scale
    ax.plot(photwaves[2:5],photvals[2:5]*1e6,ms=10,label='2MASS') # factor 1e6 makes plot in Jy, gives nice scale
    ax.plot(photwaves[5:],photvals[5:]*1e6,ms=10,label='IRAC') # factor 1e6 makes plot in Jy, gives nice scale
    ax.set_xlabel('Wavelength [micron]')
    ax.set_ylabel('Flux density [Jy]')
    ax.legend(loc='best')
    ax.set_xlim(0,22)
    ax.set_ylim(0,10)
    return(photvals)
