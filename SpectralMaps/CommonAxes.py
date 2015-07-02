#! /usr/bin/python
# Common Axes for Sophia's plots
def stax(common):
  common.set_frame_color('black')
  common.set_axis_labels_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
  common.set_tick_labels_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
  common.set_tick_labels_format(xformat='dd:mm:ss',yformat='dd:mm:ss')
  common.ticks.set_length(10)  # points
  common.ticks.set_linewidth(2)  # points
  common.ticks.set_color('black')
  common.ticks.set_minor_frequency(5)
  #
  common.add_scalebar(0.00198414)
  common.scalebar.show(0.00198414)
  common.scalebar.set_label('100pc')
  common.scalebar.set_font(size='x-large', weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
  common.scalebar.set(linestyle='solid', color='black', linewidth=4)
 
