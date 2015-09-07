from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, Column, hstack, vstack
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area # need astropy 1.0+ for this
from photutils import SkyCircularAperture, aperture_photometry
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



imglist = ['m31nuc_part_sil_map_se.fits','m31_1_bgsub_bc_nuc.fits']
#photwaves = np.array(([10,3.6,4.5,5.8,8]))
#photcorr_extd = np.array([1.0,0.91,0.94,0.68,0.74]) # IRAC extd src correction, from handbook
#photcorr_small_ap = np.array([1.0,1.07,1.08,1.076,1.087]) # IRAC pt src aperture correction, for 4-pix radius ap

def do_surf_phot(img_list=imglist, aperture_radii=np.arange(1,20,1)*u.arcsec):
    ''' Do aperture photometry on a list of images in a bunch of apertures
    input: img_list: list of images
           aperture_list: apertures to use (same for all)

    output: astropy table with photomety for all
    '''

    # make an empty table to hold the output
    rcol = Column(aperture_radii.value, name=('r_arcsec'))
    phot_tab = Table([rcol])

    # loop over the images
    for (img_num,img_name) in enumerate(img_list):

        img = fits.open(img_name)
        photlist = multi_ap_phot(img, glx_ctr, aperture_radii, sb=True)

        # calibrate photmetry: TODO
#        obs_val = out_tab['aperture_sum'][0] * out_tab['aperture_sum'].unit
#        cal_phot = calib_phot(obs_val, img, output_units='MJy')

        phot_tab.add_column(Column(photlist, name=img_name))
    # done loop over images

    return(phot_tab)

def multi_ap_phot(img, ap_ctr, aperture_radii, sb=False):

    flux = []
    for radius in aperture_radii:
        flux.append(aperture_photometry(img, SkyCircularAperture(ap_ctr, radius)))
    phot_table = vstack(flux)
    if sb:
        areas = np.zeros(len(aperture_radii))
        surf_br = np.zeros(len(aperture_radii))
        r = aperture_radii.value
        areas[0] = r[0]**2
        areas[1:] = r[1:]**2-r[:-1]**2
        areas = areas * np.pi
        surf_br[0] = phot_table['aperture_sum'][0]/areas[0]
        surf_br[1:] = (phot_table['aperture_sum'][1:]- phot_table['aperture_sum'][:-1])/areas[1:]
        for i in range(0,len(r)):
            print r[i], phot_table['aperture_sum'][i], areas[i],surf_br[i]
        return(surf_br)
    else:
        return(phot_table['aperture_sum'])



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


def makeplot(photdat):

    # plot
    f,ax=plt.subplots()
    for col in photdat.colnames[1:]:
        ax.plot(photdat['r_arcsec'],photdat[col]/photdat[col][0], ls='solid', lw=2, label=col)
    ax.set_xlabel('Radius [arcsec]')
    ax.set_ylabel('Surface brightness [relative]')
    ax.legend(loc='best')
    f.show()
    return

