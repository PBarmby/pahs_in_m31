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

# center of SL nucleus map
ap_ctr = SkyCoord(10.684029,41.269439,frame='icrs', unit='deg')
# SL extraction aperture is 50" wide, 30" high, with PA of 45 deg
phot_ap = SkyRectangularAperture(ap_ctr,w=50*u.arcsec,h=30*u.arcsec,theta=45*u.degree)

# smaller apertures also used 
ap_ctr_north = SkyCoord(10.683541,41.272717,frame='icrs', unit='deg') # see 11_3downNuc.reg
phot_ap_north = SkyRectangularAperture(ap_ctr_north,w=9*u.arcsec,h=9*u.arcsec,theta=45*u.degree)
ap_ctr_cen = SkyCoord(10.68445,41.269065,frame='icrs', unit='deg') # see m31nuc_sl1_nucCentre.tbl
phot_ap_cen = SkyRectangularAperture(ap_ctr_cen,w=9*u.arcsec,h=9*u.arcsec,theta=45*u.degree)

imglist = glob.glob('m31nuc_f1*.fits')+glob.glob('m31_2mass_*.fits')+glob.glob('m31_?_bgsub_bc_nuc.fits')    
photwaves = np.array(([1.11, 1.61, 1.6,1.2,2.2,3.6,4.5,5.8,8]))
photcorr_extd = np.array([1.0,1.0,1.0,1.0,1.0,0.91,0.94,0.68,0.74]) # IRAC extd src correction, from handbook
photcorr_small_ap = np.array([1.0,1.0,1.0,1.0,1.0,1.07,1.08,1.076,1.087]) # IRAC pt src aperture correction, for 4-pix radius ap

def dophot(img_list=imglist, aperture = phot_ap, band_waves = photwaves, apcorr = photcorr_extd):
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

def makeplot(photdat=None, specfile='../pb_m31_spectra/nucFLUX',spec_area =1500):
    '''plots photometry and some spectra on same plot'''
    if photdat == None:
        # do photometry 
        photdat = dophot(imglist)
    photwaves = photdat['Wavelength'] # known wavelengths
    photvals = photdat['MJy_counts']

#    load the IRS spectrum and convert to MJy
    nuc_wave,nuc_irs = np.loadtxt(specfile,unpack=True,usecols=[0,1])
    nuc_irs = nuc_irs*((spec_area*u.arcsec**2).to(u.sr).value)

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
#    ax.set_ylim(0,10)
    return(photvals)

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

def spect_norm(spect_wave, spect_flux, norm_wave, norm_val):
    find_norm_pt = np.searchsorted(spect_wave,norm_wave)
    norm_spect = spect_flux *(norm_val/spect_flux[find_norm_pt])
    return(norm_spect)
