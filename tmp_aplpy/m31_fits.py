#! /usr/bin/python
# Sophia Lianou - Plotting script for Dimuthu's m31 maps
import astropy
import montage_wrapper as montage
import aplpy

gc = aplpy.FITSFigure('m31_4_fixed.fits', subplot=(1, 1, 1), north=True) #
gc.axis_labels.set_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
gc.tick_labels.set_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
gc.tick_labels.set_yformat('dd.dd')
gc.tick_labels.set_xformat('dd:mm:ss')
gc.ticks.set_length(10)  # points
gc.ticks.set_linewidth(2)  # points
gc.ticks.set_color('black')
gc.ticks.set_minor_frequency(5)
#
gc.show_colorscale(cmap='YlGnBu', smooth=None, stretch='linear') # the color of the map
#
gc.add_colorbar()
gc.colorbar.show()
gc.colorbar.set_location('top')
gc.colorbar.set_width(0.2)  # arbitrary units, default is 0.2
gc.colorbar.set_pad(0.1)  # arbitrary units, default is 0.05
gc.colorbar.set_axis_label_font(size=20, weight='normal')
gc.colorbar.set_font(size='large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
gc.colorbar.set_axis_label_text('Flux density (MJy) per sr')
gc.colorbar.set_axis_label_pad(10)
gc.colorbar.set_axis_label_rotation(0)
#
gc.add_scalebar(0.00198414) # in degrees
gc.scalebar.show(0.00198414)
gc.scalebar.set_label('100pc')
gc.scalebar.set_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
gc.scalebar.set(linestyle='solid', color='black', linewidth=3)
#
gc.save('m31_4_fixed.eps')
