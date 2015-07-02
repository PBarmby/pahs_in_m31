import astropy
from astropy.io import fits
import montage_wrapper as montage
import aplpy
import pyregion
import matplotlib.pyplot as plt
# Plots 11.3 micron fits image of M31 with two regions on it.


gc = aplpy.FITSFigure('m31nuc.fits', subplot=(1, 1, 1), north=True) #
gc.show_regions('CentreNuc.reg') # upload some regions in the ds9 format
gc.show_regions('11_3downNuc.reg') # upload some regions in the ds9 format

# the following are just some format options
gc.axis_labels.set_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
gc.tick_labels.set_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
gc.tick_labels.set_yformat('dd:mm:ss')
gc.tick_labels.set_xformat('hh:mm:ss')
gc.ticks.set_length(10)  # points
gc.ticks.set_linewidth(2)  # points
gc.ticks.set_color('black')
gc.ticks.set_minor_frequency(5)
plt.rcParams['lines.marker']= 'None' # colorbar needs solid borders and no markers 
#
gc.show_colorscale(cmap='YlGnBu', smooth=None, stretch='linear', pmin=50) # the color of the map
#
gc.add_colorbar()
gc.colorbar.show()
gc.colorbar.set_location('top')
gc.colorbar.set_width(0.2)  # arbitrary units, default is 0.2
gc.colorbar.set_pad(0.22)  # arbitrary units, default is 0.05
gc.colorbar.set_axis_label_font(size=20, weight='normal')
gc.colorbar.set_font(size='large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
gc.colorbar.set_axis_label_text('Flux (W/m$^2$sr)')
gc.colorbar.set_axis_label_pad(10)
gc.colorbar.set_axis_label_rotation(0)
#
gc.add_scalebar(0.00198414) # in degrees
gc.scalebar.show(0.00198414)
gc.scalebar.set_label('100pc')
gc.scalebar.set_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
gc.scalebar.set(linestyle='solid', color='black', linewidth=3)
#
gc.add_label(0.15, 0.9, 'IRS SL1', size='x-large', family='sans-serif', relative=True)
gc.add_label(0.6, 0.4, 'Nucleus', size='x-large', family='sans-serif', relative=True)
#