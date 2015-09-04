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

# usage:
# photometry on full SL nucleus region:
#   cube_phot = nuc_phot.dophot()  # use defaults
# photometry on smaller 9arcsec regions:
#   north_phot = nuc_phot.dophot(aperture=nuc_phot.phot_ap_north,apcorr= nuc_phot.photcorr_small_ap)  
#   cen_phot = nuc_phot.dophot(aperture=nuc_phot.phot_ap_cen ,apcorr= nuc_phot.photcorr_small_ap)  

# ctr of galaxy
glx_ctr = SkyCoord(10.68475, 41.269028, frame='icrs', unit='deg') 

# photometric apertures
phot_ap_cen = SkyCircularAperture(glx_ctr,r=4*u.arcsec)

imglist = ['m31nuc_part_sil_map_se.fits']+glob.glob('m31_?_bgsub_bc_nuc.fits')    
photwaves = np.array(([10,3.6,4.5,5.8,8]))
photcorr_extd = np.array([1.0,0.91,0.94,0.68,0.74]) # IRAC extd src correction, from handbook
photcorr_small_ap = np.array([1.0,1.07,1.08,1.076,1.087]) # IRAC pt src aperture correction, for 4-pix radius ap

def do_surf_phot(img_list=imglist, aperture_list = phot_ap, band_waves = photwaves, apcorr = photcorr_extd):
    ''' Do aperture photometry on a list of images
    input: img_list: list of images
           phot_ap: aperture to use (same for all)
           band_waves: labels for wavelengths of the images in img_list
           apcorr: multiplicative photometric correction for each band

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

    # apply IRAC aperture correction
    phot_tab['MJy_counts'] = phot_tab['MJy_counts']* apcorr
    # add wavelength info
    phot_tab.add_column(Column(name='Wavelength', data= band_waves, unit='micron'))    
    phot_tab.sort('Wavelength')

    return(phot_tab)

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
        elif 'MAGZP' and 'VEGAFLUX' in hdr: # convert from mag to Jansky
            print 'magnitude to flux'
            mag = hdr['MAGZP'] - 2.5*np.log10(input_value)
            obs_val = hdr['VEGAFLUX']*10.0**(-0.4*mag) * u.Jy
        else:
            print 'Not enough info to calibrate'
    # surface-brightness to flux conversion
    elif input_value.unit == u.MJy/u.sr: #        (this is not perfectly general but oh well)
        print 'surface brightness to flux'
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


def makeplot_v2(photdat, norm_wave=8, normval = None, specfile='../pb_m31_spectra/nucFLUX',spec_area =1500):
    '''plots photometry and some spectra on same plot'''
    # difference btw this one and makeplot() is how the normalization is done

    photwaves = photdat['Wavelength'] # known wavelengths
    photvals = photdat['MJy_counts']

    if normval == None:
        normval = photvals[np.searchsorted(photwaves, norm_wave)]

#    load the IRS spectrum and convert to MJy
    nuc_wave,nuc_irs = np.loadtxt(specfile,unpack=True, usecols=[0,1])
    nuc_irs = nuc_irs*((spec_area*u.arcsec**2).to(u.sr).value)
    # normalize
    nuc_irs = spect_norm(nuc_wave,nuc_irs, norm_wave, normval)

#   read the nu Pav spectrum
    nupav_wave, nupav_flux = np.loadtxt('nu_pav_spect.txt',unpack=True,usecols=[0,1])
    # normalize
    nupav_flux = spect_norm(nupav_wave,nupav_flux, norm_wave, normval)

#   create a RJ tail to compare to
    bb_wl = np.arange(1.0,22,0.4)
    bb = blackbody_nu(bb_wl*u.micron,5000)
    # normalize
    bb = spect_norm(bb_wl,bb, norm_wave, normval)

    # plot
    f,ax=plt.subplots()
    ax.plot(bb_wl, bb.value*1e6, ls='dashed', color='k',marker=None, lw=2, label = '5000K BB' )
    ax.plot(nupav_wave, nupav_flux*1e6, ls = 'solid', color='k',marker=None, lw=2,label='nu Pav')
    ax.plot(nuc_wave,nuc_irs*1e6,ls='solid',marker=None, lw=2,label= 'M31 IRS')
    ax.plot(photwaves[5:],photvals[5:]*1e6,ms=10,label='IRAC') # factor 1e6 makes plot in Jy, gives nice scale
    ax.set_xlabel('Wavelength [micron]')
    ax.set_ylabel('Flux density [Jy]')
    ax.legend(loc='best')
    ax.set_xlim(3,22)
    ax.set_ylim(0,4)
    return

