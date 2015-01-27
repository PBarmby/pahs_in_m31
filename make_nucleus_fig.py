import aplpy
import matplotlib.pyplot as plt

def main(img_file, region_file=None, colorbar=True, scalebar=True):
    '''Makes Figure 13 of paper
       use either m31nuc_irs10.fits or m31nuc_irs11.3.fits  as img_file
    '''
    gc = aplpy.FITSFigure(img_file, subplot=(1, 1, 1), north=True) 
    
    # format options
    gc.axis_labels.set_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
    gc.tick_labels.set_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
    gc.tick_labels.set_yformat('dd:mm:ss')
    gc.tick_labels.set_xformat('hh:mm:ss')
    gc.ticks.set_length(10)  # points
    gc.ticks.set_linewidth(2)  # points
    gc.ticks.set_color('black')
    gc.ticks.set_minor_frequency(5)
    plt.rcParams['lines.marker']= 'None' # colorbar needs solid borders and no markers 
    
    gc.show_colorscale(cmap='YlGnBu', smooth=None, stretch='linear', pmin=50) # the color of the map
    if region_file!= None:
        gc.show_regions(region_file) # upload some regions in the ds9 format

    if colorbar:
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
    
    if scalebar:
        gc.add_scalebar(0.00735) 
        gc.scalebar.show(0.00735) # in degrees: 0.00735deg = 26.45arcsec = 100 pc
        gc.scalebar.set_label('100 pc')
        gc.scalebar.set_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
        gc.scalebar.set(linestyle='solid', color='black', linewidth=3)
    
    return
