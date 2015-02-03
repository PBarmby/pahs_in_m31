from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, Column
from astropy.wcs import WCS
from astropy.wcs.utils import celestial_pixel_scale
from photutils import SkyRectangularAperture, aperture_photometry
import numpy as np

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
        print out_tab

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
        pxarea = (celestial_pixel_scale(wcs).to(u.arcsec).value)**2
        intermed = input_value.to(u.Jy/u.arcsec**2) # not strictly necessary but easier to follow
        obs_val = intermed*(pxarea*u.arcsec**2)
    else:
        obs_val = input_value

    #  now do the conversion
    try:
        calib_val = obs_val.to(output_units).value
    except UnitsError:
        print 'Problem with unit conversion'
        return(None)

    return(calib_val)
