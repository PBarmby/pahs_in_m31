from astropy.table import Table, Column, join, vstack
import astropy.units as u
import numpy as np
import string, math, os

# list of lines to be extracted; would be better not to hardcode
pah_wave_lab=['PAH5.7', 'PAH6.2','PAH7.4','PAH7.6','PAH7.9', 'PAH8.3', 'PAH8.6', 'PAH10.7','PAH11.23','PAH11.33',\
'PAH12.0','PAH12.62','PAH12.69','PAH14.0','PAH16.45','PAH17.04','PAH17.39','PAH17.87'] 

pah_complex_list = {'PAH7.7': ['PAH7.4','PAH7.6','PAH7.9'], 'PAH11.3': ['PAH11.23','PAH11.33'],\
'PAH12.7': ['PAH12.62','PAH12.69'], 'PAH17.0': ['PAH16.45','PAH17.04','PAH17.39','PAH17.87']}

# from table 3 of Smith et al (2007):
# 7.7 complex: 7.4,7.6,7.9
# 11.3 complex: 11.2,11.3
# 12.7 complex: 12.6,12.7
# 17.0 complex: 16.4, 17.0,17.4,13.9

# from comments in Dimuthu code looks like 6.7,13.5, 14.2, 15.9 are not included

atomic_wave_lab=['ArII', 'ArIII', 'SIV', 'NeII', 'NeIII', 'SIII'] 

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
            if np.isnan(uncerts[i]): 
                obj_dict[line_lab] = 0           # this is what Dimuthu did for PAH lines
                obj_dict[line_lab+'_unc'] = 0    # perhaps a little dubious
            else:
                obj_dict[line_lab] = vals[i]
                obj_dict[line_lab+'_unc'] = uncerts[i]
        linetab.add_row(obj_dict)
    return(linetab)

def convert_linelist(in_tab, conv_factor = 1.0e9, complex_list = pah_complex_list, fix_upper_lim=False):
    """ process raw line list:
        multiply by conv_factor
        add lines in complexes
        add units to table header (NOT DONE)
        replace values less than their uncertainties with upper limits
    """

    tab = in_tab.copy() # returns a copy of the table, original is left unaltered
    wms = u.W/(u.m*u.m)
    # first just multiply everything by conversion factor (NaNs are OK here), add units
    for col in tab.colnames[2:]:
        tab[col] *= conv_factor
        tab[col].units = wms
    # TODO: put the factor of 1e9 in the header!

    # now compute complexes
    for complex in complex_list.keys():
        compl_val = np.zeros(len(tab))
        compl_unc = np.zeros(len(tab))
        for feat in complex_list[complex]:
            if feat not in tab.colnames:
                print 'warning: missing feature %s' % feat
                continue
            compl_val += tab[feat]
            compl_unc += (tab[feat+'_unc'])**2
        tab[complex] = compl_val
        tab[complex+'_unc'] = np.sqrt(compl_unc)
    # end of loop over complexes

    if fix_upper_lim:
        # now check each detection and see if unc> value
        for col in tab.colnames[2:-1:2]:
            tab[col][tab[col+'_unc']>tab[col]] = np.nan
            tab[col][tab[col]<1e-20] = np.nan        
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

def compute_upper_limit(spec, feature, SNR=3):

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
	
    spec = np.loadtxt("irc4FLUX")  # Load the Flux values
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
    RMS = std(flux)  # Calculates the RMS

    return RMS

def get_N(SNR,FWHM,RMS,del_lam,line):
    """Calculates the upper limit for a specific atomic line. Here we assume a Gaussian profile
    for the atomic line profile where 'FWHM' is its FWHM."""
    
    sigma = FWHM/float(2.35) # Gets  sigma from the FWHM
    print sigma
    sqrtN = sqrt(6*sigma/ float(del_lam)) # Get the square root of N. Here N is calculated by finding the
    # number of data points within a 3 sigma width (Both sides => 6 sigma).
    print sqrtN

    # F is the upper limit of the atomic line intensity.
    F = (3*(10**(-6)))*SNR*2*RMS*sqrtN*del_lam/float((line)**2)   # Unit conversion has been done [ C * 10^-26 /lambda^2 ]

    Limit = (3*RMS*FWHM/(2.35*0.3989*line**2))*10**(-6)*2.99  # This is the changed version of finding upperlimits 
    # MJy*mu/sr to W/m^2/sr

    return F, Limit





# from Dimuthu_M31_EQWs_PAHFIT/eqw_std.py
"""
Purpose : To combine EQWs of dust features and their uncertainties and normalize them.
          Functionality of each function is described within the code.

Region  : All regions of M31 except the nucleus. 

Inputs  : 1) Table of EQWs which has 18 columns (Dust features) and 500 rows (Data from 500 spectra).
          2) Data files of EQWs of 18 dust features for each region.

Outputs : 1) Data file which has the uncertainties of EQWs of 10 dust features from all regions.
          2) Data file which has the EQWS of 10 dust features from all regions.
          3) Two data files of normalized EQWs and their uncertainties.
          4) A table which has the EQW +/- Uncertainty in a readable version for Latex.

Notes  : 

"""
# Calling all the file names and getting the uncertainty of EQWs.
#for name in eqw_filenames:
#    EQW_values = np.loadtxt(name, skiprows = 1 )
#    eqw_err = get_eqw_std(EQW_values)  # finding the standard deviation of EQWs of dust features from all the regions.
#    eqw_unc.append(eqw_err)
#
#np.savetxt('eqw_unc.dat', eqw_unc,fmt='%.2f' )  # Saving a data file which has uncertainties of EQWs of 10 dust features from all the regions.
#    
#eqw_of_regions = np.loadtxt("eqw_of_regions.dat" , dtype = 'string')  # Reading all the file names of files that has eqw values of regions
#
#EQW_array = []  # This array is defined to put EQWs of all the regions.
## This For loop reads all the files from all the regions and add EQWs to EQW_array.
#for name in eqw_of_regions:
#    EQWs = np.loadtxt(name, skiprows = 1 )
#    EQWs[np.isnan(EQWs)] = 0  # Setting nan to zero
#    EQW_array.append(EQWs)
#
#EQW_array = zip(*EQW_array)  # change rows in to columns in EQW_array
#
#
#np.savetxt('EQW_original.txt', EQW_array,fmt='%.2f')  # Saving EQWs in one data file.
#
#Num_regions = np.shape(eqw_of_regions)[0]  # Number of regions 
#Num_features = 11 # Number of dust features that we are interested in
#Combined_EQW = np.zeros((Num_regions,Num_features))  # Defines an array to fill data
#
##EQW_array = np.array([EQW_array])
##np.savetxt('EQW_original.txt', EQW_array)
#
#
##print EQW_array[2]
#
#Combined_EQW[:,0] = EQW_array[0]  # 5.7 mu dust feature
#Combined_EQW[:,1] = EQW_array[1]  # 6.2 mu dust feature
#Combined_EQW[:,2] = np.array(EQW_array[2]) + np.array(EQW_array[3]) + np.array(EQW_array[4]) # 7.7 mu dust feature
#Combined_EQW[:,3] = EQW_array[5]  # 8.3 mu dust feature
#Combined_EQW[:,4] = EQW_array[6]  # 8.6 mu dust feature
#Combined_EQW[:,5] = EQW_array[7]  # 10.7 mu dust feature
#Combined_EQW[:,6] = np.array(EQW_array[8]) + np.array(EQW_array[9])  # 11.3 mu dust feature
#Combined_EQW[:,7] = EQW_array[10]  # 12.0 mu dust feature
#Combined_EQW[:,8] = np.array(EQW_array[11]) + np.array(EQW_array[12])  # 12.7 mu dust feature
#Combined_EQW[:,9] = EQW_array[13] # 14.0 mu dust feature
#Combined_EQW[:,10] = np.array(EQW_array[14]) + np.array(EQW_array[15]) + np.array(EQW_array[16]) + np.array(EQW_array[17])  # 17.0 mu dust feature
#
#
#
##Combined_EQW = np.array([Combined_EQW])
##Combined_EQW =  [ "%.2f" % i for i in Combined_EQW]
##Combined_EQW = np.around(Combined_EQW, decimals=2)
#np.savetxt('EQW_combined.dat', Combined_EQW,fmt='%.2f')
#
#
#########################################################################################
## Finding the normalization factor for EQWs and normalize both EQW and Uncertainty
## Normalization is done by finding the average EQW of each dust feature and dividing them by the average.
#Norm_Fact_Arr= []
#eqw_unc = np.loadtxt("eqw_unc.dat")
#
#for i in range(Num_features) :
#    PAHline = Combined_EQW[:,i]
#    PAHline = PAHline[np.where(PAHline!=0)]  # Removing the regions where EQW = 0
#    avrg = np.mean(PAHline)
#    if avrg != 0 :
#        Norm_fac = 1./float(avrg)
#        Norm_Fact_Arr.append(Norm_fac)
#    Combined_EQW[0:,i] = (Combined_EQW[0:,i])*(Norm_fac)
#    eqw_unc[0:,i] = (eqw_unc[0:,i])*float(Norm_fac)
#
#np.savetxt('Normslized_EQW.txt', Combined_EQW,fmt='%.2f') # Saving normalized EQWs.
#np.savetxt('Normslized_EQWunc.txt', eqw_unc,fmt='%.2f')  # Saving uncertainties of the normalized EQWs.
#
#
#
#############################################################################################
##Read two tables of values and corresponding error values and put them together with +- mark
#def put_value_err_together(value,error):
#
#    regions = np.shape(value)[0] # number of regions
#    features = np.shape(value)[1] # number of features
#    concated = np.zeros((regions,features),dtype=np.str)
#    #concated = ["               " for x in range(features)]
#    total = []
#    print concated
#    #concated = np.empty([regions,features], dtype=str)
#    for i in range(regions):
#        concated = ["                   " for x in range(features)]
#        for j in range(features):
#            concated[j] = str("%.1f" % value[i,j]) + '$\pm$' + str("%.1f" % error[i,j] + '        ')
#            print concated[j]
#        total.append(concated)
#    print np.shape(total)
#    np.savetxt('EQW_&_Error.txt', total, fmt='%s' )
#    #print total
#
#
##put_value_err_together(Combined_EQW,eqw_unc)
#UnnormalizedCombined_EQW= np.loadtxt("EQW_combined.dat")
#Unnormalizedeqw_unc= np.loadtxt("eqw_unc.dat")
##If you comment the above two lines and call the put_value_err_together(Combined_EQW,eqw_unc) it will give the normalized eqw table
#put_value_err_together(UnnormalizedCombined_EQW,Unnormalizedeqw_unc)



# first process EQW_ERR files to compute std deviation
# then put these together with measurements to make 2-column files to be read by get_linelists
def process_EQW_unc(filelist_file, prefix='EQW_ERR_', suffix='.dat', newsuffix='EQWUNC.dat', eqwsuffix='EQW.dat'):
    filelist = np.loadtxt(filelist_file, dtype='string')
    for filename in filelist: # loop over all files in filelist
        base = filename[string.find(filename,prefix)+len(prefix):string.find(filename,suffix)] #pull prefix & suffix off file name
        newf = base + newsuffix # make new output filename
        dat = np.loadtxt(filename, skiprows=1) # read data
        nlines = dat.shape[1] # find out how many columns
        std_dat = np.zeros(nlines)
        for i in range(0,dat.shape[1]): # loop over all columns in file
            std_dat[i] = dat[:,i].std() # compute std deviation
        np.savetxt(newf, std_dat, fmt = '%.5e') 
        # this is kind of hacky but I'm in a hurry
        eqw_fname = base+eqwsuffix # find the EQW file
        newf2 = base + 'EQW2' + suffix # make a new output file name
        sysstr = 'tail -18 %s | paste - %s > %s' % (eqw_fname, newf, newf2) # paste them together
        os.system(sysstr)
    return

# make all the tables
# INCOMPLETE
def doall():
    tab_eqw = get_linelists('eqw_filenames.dat',suffix='EQW2.dat',skipr=0)
    tab_pah = get_linelists('PAHfilenames.dat')
    tab_atm = get_linelists('Atomiclines_fnames',suffix='_ato_Line.dat', wave_lab=atomic_wave_lab,)
    tab_atm_new = convert_linelist(tab_atm, conv_factor = XX, complex_line_list={}, fix_upper_lim=True)
    tab_eqw_new = convert_linelist(tab_eqw, conv_factor = XX, complex_line_list=pah_complex_list, fix_upper_lim=False)
    tab_pah_new = convert_linelist(tab_pah, conv_factor = XX, complex_line_list=pah_complex_list, fix_upper_lim=False)
    tab_eqw_norm = norm_pah(tab_eqw_new)
    # join into one big table?
    # add publishable ID
    # write to FITS table
    # write to Latex table
    return

def norm_pah():
