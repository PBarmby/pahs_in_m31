import aplpy
import matplotlib.pyplot as plt

# use:
# fig10a = make_nucleus_fig.main('m31nuc_irs11.3.fits',region_file='extract_apertures.reg')
# fig10b = make_nucleus_fig.main('m31nuc_irs10.fits',region_file='nuclear_region.reg',contour_img='m31_1_bgsub_nuc.fits')
# fig10a.save('fig10a.eps',dpi=150)
# fig10b.save('fig10b.eps',dpi=150)
def main(img_file, region_file=None, scalebar=True, colorbar=False, contour_img=None):
    '''Makes Figure 10 of paper
       use either  m31nuc_irs11.3.fits or m31nuc_irs10.fits as img_file
       for regions: use extract_apertures.reg  on irs11.3
                    and nuclear_region.reg on irs10
       for contours: use m31_1_bgsub_nuc.fits
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

    if scalebar:
        gc.add_scalebar(0.00735) 
        gc.scalebar.show(0.00735) # in degrees: 0.00735deg = 26.45arcsec = 100 pc
        gc.scalebar.set_label('100 pc')
        gc.scalebar.set_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
        gc.scalebar.set(linestyle='solid', color='black', linewidth=3)

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

    if contour_img!= None:
        gc.show_contour(contour_img, colors='k', overlap=True, smooth=5)
    
    return(gc)
