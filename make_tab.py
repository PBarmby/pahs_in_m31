from astropy.table import Table, Column, join, vstack
import astropy.units as u
import numpy as np
import string, math, os


# list of lines to be extracted; would be better not to hardcode
pah_wave_lab=['PAH5.7', 'PAH6.2','PAH7.4','PAH7.6','PAH7.9', 'PAH8.3', 'PAH8.6', 'PAH10.7','PAH11.23','PAH11.33',\
'PAH12.0','PAH12.62','PAH12.69','PAH14.0','PAH16.45','PAH17.04','PAH17.39','PAH17.87'] 

pah_complex_list = {'PAH7.7': ['PAH7.4','PAH7.6','PAH7.9'], 'PAH11.3': ['PAH11.23','PAH11.33'],\
'PAH12.7': ['PAH12.62','PAH12.69'], 'PAH17.0': ['PAH16.45','PAH17.04','PAH17.39','PAH17.87']}
# Table 3 of Smith et al (2007) lists all the lines and gives complexes as:
# 7.7 complex: 7.4,7.6,7.9
# 11.3 complex: 11.2,11.3
# 12.7 complex: 12.6,12.7
# 17.0 complex: 16.4, 17.0,17.4,13.9
# from comments in Dimuthu code looks like 6.7,13.5, 14.2, 15.9 are not included in PAHFIT output

# atomic line list
atomic_wave_lab=['ArII', 'ArIII', 'SIV', 'NeII', 'NeIII', 'SIII'] 

startcol = 2 # first table column that contains numeric values
upperlim_tol=1e-20

# lists of columns to output in paper tables, with formatting info
atm_cols = [('Pub_ID','%12s'), ('ArII','%.1f'), ('ArIII', '%.1f'), ('SIV','%.1f'), ('NeII','%.1f'),('NeIII','%.1f'),('SIII','%.1f')]
pah_cols=[('Pub_ID','%12s'),('PAH5.7','%.1f'), ('PAH6.2','%.1f'),('PAH7.7','%.1f'),('PAH8.3','%.1f'),('PAH8.6','%.1f'), ('PAH10.7','%.1f'),\
('PAH11.3','%.1f'),('PAH12.0','%.1f'),('PAH12.7','%.1f'),('PAH17.0','%.1f')]

# make all the tables
# INCOMPLETE
def doall(writefits=False, write_latex=False):
    tab_eqw = get_linelists('eqw_filenames.dat',suffix='EQW2.dat',skipr=0, wave_lab=pah_wave_lab)
    tab_pah = get_linelists('PAHfilenames.dat', suffix='PAH.dat',skipr=1,wave_lab=pah_wave_lab)
    tab_atm = get_linelists('Atomiclines_fnames',suffix='_ato_Line.dat', skipr=1,wave_lab=atomic_wave_lab)
    # conversion factor 35.26 goes from W/m^2/sr to nW/m^2, assuming 1500arcsec^2 extraction area.
#    cf = 35.26
    cf= 1e9
    nwms = 1.0e-9* u.W/(u.m*u.m)
    tab_atm_new = convert_linelist(tab_atm, conv_factor = cf, complex_list={}, fix_upper_lim=True, colunit=nwms)
    tab_pah_new = convert_linelist(tab_pah, conv_factor = cf, complex_list=pah_complex_list, fix_upper_lim=False, colunit=nwms)
    tab_eqw_new = convert_linelist(tab_eqw, conv_factor = 1.0, complex_list=pah_complex_list, fix_upper_lim=False, colunit=u.micron)
#    tab_eqw_norm = norm_pah(tab_eqw_new) # TODO: figure out why this is failing
    tab_atm_new = add_pub_id(tab_atm_new, "id_map")
    tab_eqw_new = add_pub_id(tab_eqw_new, "id_map")
    tab_pah_new = add_pub_id(tab_pah_new, "id_map")
#    add_pub_id(tab_eqw_norm, "id_map")
    if writefits: # write to FITS tables
        tab_atm_new.write('m31_atomic.fits', format='fits')
        tab_eqw_new.write('m31_pah_eqw.fits', format='fits')
        tab_pah_new.write('m31_pah_str.fits', format='fits')
#        tab_eqw_norm.write('m31_pah_eqw_norm.fits', format='fits')
    # write to Latex tables
    if write_latex:
        make_latex_table_rows(tab_atm_new, col_list = atm_cols, outfile = 'm31_atomic.tex')
        make_latex_table_rows(tab_eqw_new, col_list = pah_cols, outfile = 'm31_pah_eqw.tex')
        make_latex_table_rows(tab_pah_new, col_list = pah_cols, outfile = 'm31_pah_str.tex')
    # join into two big tables
#    big_tab1 = join(tab_pah_new, tab_atm_new, table_names = ['PAH_str','Atom'])
#    big_tab2 = join(tab_eqw_new, tab_eqw_norm)
    return



# mod from m31gemini/analysis/table_proc.py
def get_linelists(filelist_file, suffix='PAH.dat', wave_lab=pah_wave_lab, skipr=1):
    """generate a table combining linelists for different objects"""

    # first set up the empty table
    collist = ['ID', 'Filename'] 
    for i,wavelength in enumerate(wave_lab): # make column names
        collist += [wave_lab[i],wave_lab[i]+'_unc']
    dt_list = ['a20'] * 2 + ['f4'] * (len(collist) -2)  # list of data types: one string plus bunch of floats
    linetab = Table(names=collist, dtype = dt_list) # generate a new table

    filelist = np.loadtxt(filelist_file, dtype='string')

    # now put stuff in the table
    for f in filelist: # loop through each file
        vals, uncerts = np.loadtxt(f, unpack=True,usecols=[0,1],skiprows=skipr)
        obj_dict={}
        obj_dict['Filename'] = f 
        obj_dict['ID'] = f[:string.find(f,suffix)]
        for i,line_lab in enumerate(wave_lab): # put the columns from the file into one row of a table
            if np.isnan(uncerts[i]) or np.abs(vals[i] < upperlim_tol) :  # non-detection or upper limit
                obj_dict[line_lab] = 0.0           
                obj_dict[line_lab+'_unc'] = np.nan    
            else:
                obj_dict[line_lab] = vals[i]
                obj_dict[line_lab+'_unc'] = uncerts[i]
        linetab.add_row(obj_dict)
    return(linetab)

def convert_linelist(in_tab, conv_factor = 1.0e9, complex_list = pah_complex_list, fix_upper_lim=False, colunit='None'):
    """ process raw line list:
        multiply by conv_factor
        add lines in complexes
        add units to table header
        replace values less than their uncertainties with upper limits
    """

    tab = in_tab.copy() # returns a copy of the table, original is left unaltered

    # fix upper limits first
    if fix_upper_lim:
        for i in range(0, len(tab)): # loop over rows
            thisrow = tab[i]
            specfile = 'spectra/%sFLUX' % thisrow['ID']
            for col in thisrow.colnames[startcol:-1:2]: #  check each detection and see if unc> value
                if thisrow[col+'_unc'] > thisrow[col] or np.isnan(thisrow[col+'_unc']): # applies to atomic lines
                    thisrow[col+'_unc'] = np.nan
                    thisrow[col] = compute_upper_limit(specfile, col)
 
    # then just multiply everything by conversion factor (NaNs are OK here), add units
    for col in tab.colnames[startcol:]:
        tab[col] *= conv_factor
        tab[col].units = colunit
        tab[col].format = '%.3e'
    # TODO: put the conversion factor in the header!

    # now compute complexes
    for complex in complex_list.keys():
        compl_val = np.zeros(len(tab))
        compl_unc = np.zeros(len(tab))
        for feat in complex_list[complex]:
            if feat not in tab.colnames:
                print 'warning: missing feature %s' % feat
                continue
            compl_val += tab[feat]
            compl_unc += (tab[feat+'_unc'])**2 # NB: problems here with non-detections
        tab[complex] = compl_val
        tab[complex+'_unc'] = np.sqrt(compl_unc)
    # end of loop over complexes

    return(tab)


# originally from Dimuthu_M31_Atomic_lines/limit.py
#Purpose : To calculate the upperlimits for missing atomic lines in M31 IRS spectra
#          S/N = Signal to nose ratio (We used 3)
#          RMS = Root mean square of the noise of the spectrum at position of the atomic line.
#          N = Number of data points within 3 sigma width of a Gaussian profile.
#          del_lam = wavelength difference between two data points. (This was found manually)
#Region  : Works on all regions of M31.
#
#Inputs  : Data file which has wavelength (in mu) and flux density values (in MJy/sr) of a mid_IR spectrum.
#
#Outputs : Upper limit value of the line intensity of the specified atomic line.
#
#Notes   : This program works only to find upper limits for ArII, ArIII, NeII, NeIII,  SIV 
#

def compute_upper_limit(specname, feature, SNR=3):

    # Selecting wavelength range.
    if feature== "NeIII" :
	    line = 15.5 # Position of the line
	    FWHM = 0.126  # Full Width at Half Maxima (Found using PAHFIT) in microns
	    low_lim = 15.0  # Beginning of the wavelength range of the atomic line.
	    up_lim = 16.0   # End of the wavelength range of the atomic line.
	    del_lam = 0.0917  # in microns
    elif feature=="SIV" :
	    line = 10.5
	    FWHM = 0.09
	    low_lim =10.0
	    up_lim = 11.0
	    del_lam = 0.062
    elif feature=="NeII" :
	    line = 12.8
	    FWHM = 0.09
	    low_lim =12.5
	    up_lim = 13.5
	    del_lam = 0.062
    elif feature=="ArII" :
	    line = 7.0
	    FWHM = 0.0476
	    low_lim =6.5
	    up_lim = 7.5
	    del_lam = 0.031
    elif feature=="ArIII" :
	    line = 9.0
	    FWHM = 0.09
	    low_lim =8.5
	    up_lim = 9.5
	    del_lam = 0.062
    else:
        print 'Feature must be one of NeII, NeIII, ArII, ArIII, SIV'
        return
	
    spec = np.loadtxt(specname)  # Load the Flux values
    RMS = slice_spec(spec,low_lim, up_lim)
    F,Limit = get_N(SNR,FWHM,RMS,del_lam,line)
    return(Limit)


#I got this code (boxSmooth1D) from Dr. Pauline Barmby.
def boxSmooth1D(yin, nsmooth):
    '''Smooth an input array yin using a boxcar average of width nsmooth. 
    Boxcar averaging replaces each data value yin[i] by the average of
    (adjacent) data values in the range i-nsmooth/2 to i+nsmooth/2.'''
  
    if not isinstance(nsmooth, int) or nsmooth<1:
        print 'nsmooth must be a positive integer!'
        return
    if not isinstance(yin, np.ndarray):
        print 'yin must be an array!'
        return
    if nsmooth>len(yin):
        print 'boxSmooth1D called with nsmooth > length of input array'
        return

    nlow = int(math.floor(nsmooth/2.0))
    nhigh = int(math.ceil(nsmooth/2.0))

    ysmooth = np.zeros(len(yin))
    for n in range(0,len(yin)):
        # these 3 lines correctly deal with the edges of the input array
        nelem_lower = max(0,n-nlow)
        nelem_upper=min(len(yin)-1,n+(nhigh-1))
        num = nelem_upper-nelem_lower+1
	#print "num:",num
      #  print "n= %d: nlower = %d, nupper =%d, num = %d)\n" % (n, nelem_lower, nelem_upper, num)
        for ns in range(nelem_lower,nelem_upper+1):
                   ysmooth[n]=ysmooth[n]+yin[ns]
        ysmooth[n]=ysmooth[n]/num
    return ysmooth


# Reads a file and gets the flux and uncertainty and slice it and return the RMS
def slice_spec(spec,low_lim, up_lim):
    """ This function slices a spectrum at a given wavelength range and
    gets the continum subtracted spectrum and finally calculates the
    Root Mean Square (RMS) of the spectrum. The continuum has been determined by
    smoothing the same spectrum using the boxSmooth1D() function."""
    
    wave = spec[:,0]
    flux = spec[:,1]
    flux = flux[(wave >= low_lim) & (wave <= up_lim)]
    wave = wave[(wave >= low_lim) & (wave <= up_lim)]
    smoothed = boxSmooth1D(flux, 10)
#    plt.plot(wave,flux,'-')  # Plots the sliced spectrum
#    plt.plot(wave,smoothed,'-') # Plots the smoothes spectrum
#    plt.show()

    flux = flux - smoothed  # Continuum subtraction
#    plt.plot(wave,flux,'-')  # Plots the continuum subtracted spectrum,
    RMS = np.std(flux)  # Calculates the RMS

    return RMS

def get_N(SNR,FWHM,RMS,del_lam,line):
    """Calculates the upper limit for a specific atomic line. Here we assume a Gaussian profile
    for the atomic line profile where 'FWHM' is its FWHM."""
    
    sigma = FWHM/float(2.35) # Gets  sigma from the FWHM
#    print sigma
    sqrtN = math.sqrt(6*sigma/ float(del_lam)) # Get the square root of N. Here N is calculated by finding the
    # number of data points within a 3 sigma width (Both sides => 6 sigma).
#    print sqrtN

    # F is the upper limit of the atomic line intensity.
    F = (3*(10**(-6)))*SNR*2*RMS*sqrtN*del_lam/float((line)**2)   # Unit conversion has been done [ C * 10^-26 /lambda^2 ]

    Limit = (3*RMS*FWHM/(2.35*0.3989*line**2))*10**(-6)*2.99  # This is the changed version of finding upperlimits 
    # MJy*mu/sr to W/m^2/sr

    return F, Limit



# first process EQW_ERR files to compute std deviation
# then put these together with measurements to make 2-column files to be read by get_linelists
def process_EQW_unc(filelist_file='eqwUNC_filenames.dat', prefix='EQW_ERR_', suffix='.dat', newsuffix='EQWUNC.dat', eqwsuffix='EQW.dat'):
    filelist = np.loadtxt(filelist_file, dtype='string')
    for filename in filelist: # loop over all files in filelist
        base = filename[string.find(filename,prefix)+len(prefix):string.find(filename,suffix)] #pull prefix & suffix off file name
        newf = base + newsuffix # make new output filename
        dat = np.loadtxt(filename, skiprows=1) # read data
        nlines = dat.shape[1] # find out how many columns
        std_dat = np.zeros(nlines)
        for i in range(0,dat.shape[1]): # loop over all columns in file
            std_dat[i] = dat[:,i].std() # compute std deviation
            if std_dat[i] < upperlim_tol:
                std_dat[i] = np.nan
        np.savetxt(newf, std_dat, fmt = '%.5e') 
        # this is kind of hacky but I'm in a hurry
        eqw_fname = base+eqwsuffix # find the EQW file
        newf2 = base + 'EQW2' + suffix # make a new output file name
        sysstr = 'tail -18 %s | paste - %s > %s' % (eqw_fname, newf, newf2) # paste them together
        os.system(sysstr)
    return

# NB: probably this should weight by uncertainty.
def norm_pah(in_tab):
    """ produce a normalized PAH line list:
        divide each feature + uncertainty by the average of the
        feature strength over all regions
    """

    tab = in_tab.copy() # returns a copy of the table, original is left unaltered

    # loop over all the columns
    for col in tab.colnames[startcol:-1:2]:
        normfact = in_tab[col][in_tab[col]>0].mean() # compute the average over the good values
#        print col, normfact
        tab[col] *= 1.0 / normfact  # divide by the average
        tab[col+'_unc'] *= 1.0 / normfact
        tab.rename_column(col, col+'norm') # rename the columns so we know what we did
        tab.rename_column(col+'_unc', col+'norm_unc')
    return(tab)

# from /Volumes/data/m31gemini/analysis/table_proc.py
# returns correctly-formatted latex string for a single table cell
#   including uncertainties if applicable
def uncert_str(tab_row, col_name, value_fmt):
    if col_name+'_pe' in tab_row.colnames and col_name+'_me' in tab_row.colnames: # two-sided uncertainties
    	formatter_str = '$' + value_fmt + '^{+' + value_fmt +'}_{-' + value_fmt+ '}$'
    	final_str = formatter_str % (tab_row[col_name], tab_row[col_name+'_pe'], tab_row[col_name+'_me'])
    elif col_name+'_unc' in tab_row.colnames: # one-sided uncertainties
        if np.isnan(tab_row[col_name+'_unc']): # either an upper limit or no data
            if np.isnan(tab_row[col_name]) or tab_row[col_name]<upperlim_tol : # no data
                final_str = '\\dots'                
            else: # upper limit
                formatter_str = '$<' + value_fmt +'$'
                final_str = formatter_str % (tab_row[col_name])
        else: # regular uncertainty
            formatter_str = '$' + value_fmt + '\\pm' + value_fmt +'$'
            final_str = formatter_str % (tab_row[col_name], tab_row[col_name+'_unc'])

    else: # no uncerts
        if 's' in value_fmt: # it's a string, do some replacements
                newval = string.replace(tab_row[col_name],'_','\\_')
                if 'Ref' in col_name:
                    formatter_str =  '\\citet{%s}'
                    newval = string.strip(newval)
                else:
                    formatter_str =  '' + value_fmt +''
        else:
            formatter_str =  '$' + value_fmt +'$' 
            newval = tab_row[col_name]
        final_str = formatter_str % (newval)
    return(final_str)

# contruct latex-formatted table rows
# (need this to deal with uncertainty columns)
def make_latex_table_rows(intab, col_list, outfile):

    outf = open(outfile,'w') # overwrites input
    # keep track of what columns are in the output
    formatted_line = '% '
    for j in range (0,len(col_list)):
            formatted_line += ' %s ' % col_list[j][0] # column names
    outf.write(formatted_line + '\n')
    # now write the individual rows
    for i in range(0,len(intab)):
        formatted_line = ''
        for j in range (0,len(col_list)):
            formatted_line += uncert_str(intab[i],col_list[j][0],col_list[j][1]) # one column entry
            if j<len(col_list)-1: formatted_line += ' & ' # column separator
        formatted_line +='\\\\\n' # end-of-line marker
        outf.write(formatted_line)
    outf.close()
    return

# reads ID file and uses this to insert 'publishable ID' column into tab
def add_pub_id(tab, id_file):
    tab_ids = Table.read(id_file, format='ascii.commented_header')
    tab = join(tab_ids, tab, keys='ID')
    # should error-trap the case where not all objects in tab are listed in tab_ids
    return(tab)
