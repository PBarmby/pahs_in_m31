#! /usr/bin/python
# Sophia Lianou - Plotting script for Dimuthu's spectral maps
import os
import sys
import os.path
import numpy
import scipy
import math
import matplotlib
import pylab
import pyfits
import astropy
import montage_wrapper as montage
import aplpy
import pyregion
#
from astropy.io import fits
import CommonAxes
#
#
# SL1
gc = aplpy.FITSFigure('fits/9micronIntFit.fits', subplot=(1, 1, 1), north=True) 
CommonAxes.stax(gc)
#
gc.show_colorscale(vmin=0,vmax=12, cmap='YlGnBu', smooth=None, stretch='linear') #
gc.add_colorbar()
gc.colorbar.show()
gc.colorbar.set_location('top')
gc.colorbar.set_width(0.2)  # arbitrary units, default is 0.2
gc.colorbar.set_pad(0.1)  # arbitrary units, default is 0.05
gc.colorbar.set_axis_label_font(size=20, weight='normal')
gc.colorbar.set_font(size='large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
gc.colorbar.set_axis_label_text('Flux density (MJy/sr)')
gc.colorbar.set_axis_label_pad(10)
gc.colorbar.set_axis_label_rotation(0)
#
gc.add_label(0.15, 0.95, 'IRS SL1', color='b', size='x-large', family='sans-serif', relative=True)
#
gc.save('IRS_SL1.eps')
#
###
# SL1
gc = aplpy.FITSFigure('fits/m31nuc_contu.fits', subplot=(1, 1, 1), north=True) 
gc.show_regions('fits/11_3downNuc.reg')
gc.show_regions('fits/CentreNuc.reg')
CommonAxes.stax(gc)
#
gc.show_colorscale(vmin=0, cmap='YlGnBu', smooth=1, kernel='gauss',stretch='linear') #
gc.add_colorbar()
gc.colorbar.show()
gc.colorbar.set_location('top')
gc.colorbar.set_width(0.2)  # arbitrary units, default is 0.2
gc.colorbar.set_pad(0.1)  # arbitrary units, default is 0.05
gc.colorbar.set_axis_label_font(size=20, weight='normal')
gc.colorbar.set_font(size='large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
gc.colorbar.set_axis_label_text('Flux density (W/m$^2$/sr)')
gc.colorbar.set_axis_label_pad(10)
gc.colorbar.set_axis_label_rotation(0)
#
gc.add_label(0.15, 0.95, 'IRS SL1', color='b', size='x-large', family='sans-serif', relative=True)
#
gc.save('IRS_SL1_2.eps')
#
###
# ISOCAM
gc = aplpy.FITSFigure('fits/103.fits', subplot=(1, 1, 1), north=True) 
gc.show_regions('fits/region4cubes.reg')
CommonAxes.stax(gc)
#
gc.show_colorscale(vmin=0, cmap='YlGnBu', smooth=1, kernel='gauss',stretch='linear') #
gc.add_colorbar()
gc.colorbar.show()
gc.colorbar.set_location('top')
gc.colorbar.set_width(0.2)  # arbitrary units, default is 0.2
gc.colorbar.set_pad(0.1)  # arbitrary units, default is 0.05
gc.colorbar.set_axis_label_font(size=20, weight='normal')
gc.colorbar.set_font(size='large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
gc.colorbar.set_axis_label_text('Flux density (MJy/sr)')
gc.colorbar.set_axis_label_pad(10)
gc.colorbar.set_axis_label_rotation(0)
#
gc.add_label(0.15, 0.95, 'ISOCAM', color='b', size='x-large', family='sans-serif', relative=True)
#
gc.save('Isocam.eps')
 
