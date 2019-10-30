""" Module for dealing with spectral line data
classes: Line, LineDataHeader, LineData, Accum, NetCDFLineHeader
methods: read_obsnum_ps, read_obsnum_otf, count_otf_spectra, 
uses: 
author: FPS                                                                                            date: September 2019 (documented)
changes:
python 3
"""

import numpy as np
import math
from scipy.interpolate import interp1d
from lmtslr.spec.spec import SpecBankCal, SpecBankData
from lmtslr.ifproc.ifproc import IFProcData, IFProcCal
from lmtslr.utils.roach_file_utils import create_roach_list, \
    lookup_roach_files
from lmtslr.utils.ifproc_file_utils import lookup_ifproc_file
import netCDF4

def read_obsnum_ps(obsnum,list_of_pixels,bank,use_calibration,tsys=150.,path='/data_lmt/'):
    """ Reads the spectral line data from WARES spectrometer for a particular obsnum.

        This is a useful utility function for use in line processing. It combines a number
        of tasks that are commonly done together to read a set of spectral line data.

    Args:
        obsnum (int): sequential id (ObsNum) of the observation
        list_of_pixels (list): identifies which elements of the spectrometer are to be read
        use_calibration (bool): set True if we are to use calibration scan for cal. 
                                    False just multiplies by system temperature
        tsys (float): system temperature to use of use_calibration is False
        path (string): path to the top of the data_lmt directory (usually '/data_lmt/')
    Returns:
        I: ifproc_data object with IF Processor parameters
        S: spec_bank_data object with the actual spectral data
    """
    # look up files to match pixel list
    roach_list = create_roach_list(list_of_pixels)
    files,nfiles = lookup_roach_files(obsnum,roach_list,path=path+'spectrometer/')
    ifproc_file = lookup_ifproc_file(obsnum,path=path+'ifproc/')

    # create the spec_bank object.  This reads all the roaches in the list "files"
    I = IFProcData(ifproc_file)
    S = SpecBankData(files,I,pixel_list=list_of_pixels,bank=bank)

    # check whether to use calibration and open necessary file
    if use_calibration == True:
        S.cal_flag = False
        calobsnum = S.calobsnum
        cal_files,ncalfiles = lookup_roach_files(calobsnum,roach_list,path=path+'spectrometer/')
        ifproc_cal_file = lookup_ifproc_file(calobsnum,path=path+'ifproc/')
        ICal = IFProcCal(ifproc_cal_file)
        SCal = SpecBankCal(cal_files,ICal,pixel_list=list_of_pixels)
        check_cal = SCal.test_cal(S)
        if check_cal > 0:
            print('WARNING: CAL MAY NOT BE CORRECT')

        # reduce all spectra - calibrated
        for ipix in range(S.npix):
            S.roach[ipix].reduce_ps_spectrum(type=2,normal_ps=True,calibrate=True,tsys_spectrum=SCal.roach[ipix].tsys_spectrum)
    else:
        # reduce all spectra - uncalibrated
        for ipix in range(S.npix):
            S.roach[ipix].reduce_ps_spectrum(type=2,normal_ps=True,calibrate=False,tsys_no_cal=tsys)

    return I,S

def read_obsnum_bs(obsnum,list_of_pixels,bank,use_calibration,tsys=150.,path='/data_lmt/'):
    """ Reads the spectral line data from WARES spectrometer for a particular obsnum.

        This is a useful utility function for use in line processing. It combines a number
        of tasks that are commonly done together to read a set of spectral line data.

    Args:
        obsnum (int): sequential id (ObsNum) of the observation
        list_of_pixels (list): identifies which elements of the spectrometer are to be read - for BS observation this is just the two pixels used for the switch
        use_calibration (bool): set True if we are to use calibration scan for cal. 
                                    False just multiplies by system temperature
        tsys (float): system temperature to use of use_calibration is False
        path (string): path to the top of the data_lmt directory (usually '/data_lmt/')
    Returns:
        I: ifproc_data object with IF Processor parameters
        S: spec_bank_data object with the actual spectral data
    """
    # look up files to match pixel list
    roach_list = create_roach_list(list_of_pixels)
    files,nfiles = lookup_roach_files(obsnum,roach_list,path=path+'spectrometer/')
    ifproc_file = lookup_ifproc_file(obsnum,path=path+'ifproc/')

    # create the spec_bank object.  This reads all the roaches in the list "files"
    I = IFProcData(ifproc_file)
    S = SpecBankData(files,I,pixel_list=list_of_pixels,bank=bank)

    # check whether to use calibration and open necessary file
    if use_calibration == True:
        S.cal_flag = False
        calobsnum = S.calobsnum
        cal_files,ncalfiles = lookup_roach_files(calobsnum,roach_list,path=path+'spectrometer/')
        ifproc_cal_file = lookup_ifproc_file(calobsnum,path=path+'ifproc/')
        ICal = IFProcCal(ifproc_cal_file)
        SCal = SpecBankCal(cal_files,ICal,pixel_list=list_of_pixels)
        check_cal = SCal.test_cal(S)
        if check_cal > 0:
            print('WARNING: CAL MAY NOT BE CORRECT')

        # reduce the two spectra - calibrated 
        S.roach[0].reduce_ps_spectrum(type=2,normal_ps=False,calibrate=True,tsys_spectrum=SCal.roach[0].tsys_spectrum)
        S.roach[1].reduce_ps_spectrum(type=2,normal_ps=True,calibrate=True,tsys_spectrum=SCal.roach[1].tsys_spectrum)

    else:
        # reduce the two spectra - uncalibrated
        S.roach[0].reduce_ps_spectrum(type=2,normal_ps=False,calibrate=False, tsys_no_cal=tsys)
        S.roach[1].reduce_ps_spectrum(type=2,normal_ps=True,calibrate=False,tsys_no_cal=tsys)

    return I,S

def read_obsnum_otf(obsnum,list_of_pixels,bank,use_calibration,tsys=150.,path='/data_lmt/'):
    """ Reads the spectral line data from WARES spectrometer for a particular obsnum.

        This is a useful utility function for use in line processing. It combines a number
        of tasks that are commonly done together to read a set of spectral line data.

    Args:
        obsnum (int): sequential id (ObsNum) of the observation
        list_of_pixels (list): identifies which elements of the spectrometer are to be read
        use_calibration (bool): set True if we are to use calibration scan for cal. 
                                    False just multiplies by system temperature
        tsys (float): system temperature to use of use_calibration is False
        path (string): path to the top of the data_lmt directory (usually '/data_lmt/')
    Returns:
        I: ifproc_data object with IF Processor parameters
        S: spec_bank_data object with the actual spectral data
    """
    # look up files to match pixel list
    roach_list = create_roach_list(list_of_pixels)
    files,nfiles = lookup_roach_files(obsnum,roach_list,path=path+'spectrometer/')
    ifproc_file = lookup_ifproc_file(obsnum,path=path+'ifproc/')

    # create the spec_bank object.  This reads all the roaches in the list "files"
    I = IFProcData(ifproc_file)
    S = SpecBankData(files,I,pixel_list=list_of_pixels,bank=bank)

    # check whether to use calibration and open necessary file
    if use_calibration == True:
        S.cal_flag = False
        calobsnum = S.calobsnum
        cal_files,ncalfiles = lookup_roach_files(calobsnum,roach_list,path=path+'spectrometer/')
        ifproc_cal_file = lookup_ifproc_file(calobsnum,path=path+'ifproc/')
        ICal = IFProcCal(ifproc_cal_file)
        SCal = SpecBankCal(cal_files,ICal,pixel_list=list_of_pixels)
        check_cal = SCal.test_cal(S)
        if check_cal > 0:
            print('WARNING: CAL MAY NOT BE CORRECT')

        # reduce all spectra - calibrated
        for ipix in range(S.npix):
            S.roach[ipix].reduce_spectra(type=1,calibrate=True,tsys_spectrum=SCal.roach[ipix].tsys_spectrum)
    else:
        # reduce all spectra - uncalibrated
        for ipix in range(S.npix):
            S.roach[ipix].reduce_spectra(type=1,calibrate=False,tsys_no_cal=tsys)

    return I,S

def count_otf_spectra(S, list_of_pixels):
    """ count the total number of spectra in a specbank 
    inputs:
    S - a spec object with the roach data.
    list_of_pixels - list of pixel ID's to be processed.
    """
    total_spectra = 0
    for ipix in list_of_pixels:
        i = S.find_pixel_index(ipix)
        n_spectra = len(S.roach[i].xmap[S.roach[i].ons])
        total_spectra = total_spectra + n_spectra
        print(ipix,n_spectra,total_spectra)
    print('Total Number of Spectra = %d'%(total_spectra))
    return total_spectra



class Line(object):
    """ Class with methods that operate on spectral line data.
        
        The methods here are common to any spectral line reduction.  They
        assume a common format with arrays of data for the x axis and y axis of
        the spectrum.
    """
    def __init__(self,iarray,xarray,yarray,xlabel):
        """ initializes an instance of the Line class

        Args:
            iarray (np array): array of channel numbers
            xarray (np array): array of "x axis" values.  Could be velocity, frequency, etc.
            yarray (np array): spectral line data values corresponding to xarray/iarray.
            xlabel (string): label for the x axis; used in plots
        """
        self.nchan = len(xarray)
        self.iarray = iarray
        self.xarray = xarray
        self.yarray = yarray
        self.dxdc = self.xarray[1]-self.xarray[0]
        self.xname = xlabel

    def eliminate(self,list):
        """ sets yarray values to nan so that channels can be ignored in functions
        
        Args:
            list (list): list of channels that are to be eliminated
        Result:
            yarray (np.array): yarray is modified by setting eliminated channels to np.nan
        """ 
        for i in list:
            self.yarray[np.where(self.iarray==i)] = np.nan

    def baseline(self,list,n,baseline_order=0):
        """ removes a baseline from the working spectrum

        Args:
            list (list): list of lists specifies indices in the spectrum where baseline is to be fit
            n (int): not used
            baseline_order (int): order of polynomial to be fit to baseline (used in np.polyfit)
        Result:
            baseline (np array): polynomial from baseline fit
            rms (float): rms of data about polynomial fit in regions
            yarray (np.array): yarray is modified by subtracting the baseline polynomial
        """
        input_baseline_list = np.array(list)
        baseline_list = input_baseline_list[np.where(np.isfinite(self.yarray[input_baseline_list])==True)[0]]
        n_baseline_list = len(baseline_list)
        params = np.polyfit(baseline_list,self.yarray[baseline_list],baseline_order)
        self.baseline = np.polyval(params,np.arange(self.nchan))
        residuals = self.yarray[baseline_list] - np.polyval(params,baseline_list)
        self.rms = np.sqrt(np.dot(residuals.transpose(),residuals)/n_baseline_list)
        self.yarray = self.yarray - self.baseline
    
    # special definition of generic line_statistics method
    def line_stats(self,list,n):
        """ computes line statistics for a spectrum                                
        Args:
            list (list): indices in spectrum for line integration
            n (int): not used
        Line Statistics are:                                                                               
            ymax = maximum value in channel_list                                                               
            cmax = channel of maximum value                                                                    
            xmax = velocity of maximum value                                                                   
            yint = integral of line over channel list in units of K km/s                                       
            yerr = error on YINT                                                                               
            xmean = first moment of line integrated over channel list    
            xwidth = estimate of line width based on Peak value and integrated intensity                       
        """
        input_channel_list = np.array(list)
        channel_list = input_channel_list[np.where(np.isfinite(self.yarray[input_channel_list])==True)[0]]
        n_channel_list = len(channel_list)
        self.ymax = np.max(self.yarray[channel_list])
        self.cmax = channel_list[np.where(self.yarray[channel_list] == self.ymax)[0]]
        self.xmax = self.xarray[self.cmax]
        self.yint = np.sum(self.yarray[channel_list])*self.dxdc
        self.yerr = self.dxdc*self.rms*np.sqrt(n_channel_list)
        self.xmean = np.sum(self.yarray[channel_list]*self.xarray[channel_list])/np.sum(self.xarray[channel_list])
        self.xwidth = self.yint/self.ymax

    def smo(self,nw,type='box'):
        """ smooths a spectrum using a smoothing function

        Args:
            nw (int): number of channels in smoothing function
            type (string): identifies smoothing function
                           'box' (default) - boxcar smoothing function
                           'hanning' - hanning window
                           'hamming' - hamming window
        Result:
            yarray replaced with smoothed version of spectrum sampled at original spacing
        """
        if type == 'box':
            window = np.ones(nw)
        elif type == 'hanning':
            window = np.hanning(nw)
        elif type == 'hamming':
            window = np.hamming(nw)
        else:
            print('smo: window %s not supported'%(type))
        denom = np.sum(window)
        self.yarray = np.convolve(self.yarray,window,mode='same')/denom

    def xlist(self,regions):
        """ generates a list of lists with indices corresponding to x axis values
        Args:
            regions (list): a list of lists of regions in spectrum identified according to x axis values
        Returns:
            channel_list (list): list of lists with channel numbers corresponding to input regions
            nchannels (int): number of channels in the list
        """
        nregions = len(regions)
        if self.xarray[1] > self.xarray[0]:
            channel_list = []
            for i in range(nregions):
                idx = np.where( (self.xarray>regions[i][0]) ==  (self.xarray<regions[i][1]) )[0]
                channel_list = channel_list + idx.tolist()
        else:
            channel_list = []
            for i in range(nregions):
                idx = np.where( (self.xarray>regions[i][1]) ==  (self.xarray<regions[i][0]) )[0]
                channel_list = channel_list + idx.tolist()

        nchannels = len(channel_list)
        return channel_list, nchannels

    def clist(self,regions):
        """ generates a list of lists with indices corresponding to channel number
        Args:
            regions (list): a list of lists of channels in spectrum
        Returns:
            channel_list (list): list of lists with indices corresponding to input regions
            nchannels (int): number of channels in the list
        """
        nregions = len(regions)
        channel_list = []
        for i in range(nregions):
            idx = np.where( (self.iarray>regions[i][0]) == (self.iarray<regions[i][1]) )[0]
            channel_list = channel_list + idx.tolist()
        nchannels = len(channel_list)
        return channel_list, nchannels


class LineDataHeader(object):
    ''' creates header information for line observation and provides methods to transform "x" axis.
    units of velocity: km/s
    units of frequency: Hertz
    '''
    def __init__(self,ifproc,bank,nchan,bandwidth):
        ''' constructor for LineDataHeader
        inputs:
        ifproc - an ifproc object with the ifproc information
        bank - selects which bank of the spectrometer is relevant 
        nchan - number of channels in the spectrum
        bandwidth - bandwidth of the spectrum (MHz)
        '''
        self.c = 299792.458 # speed of light in km/s
        self.bank = bank
        self.v_source = ifproc.velocity # km/s
        self.v_system = ifproc.velocity_system
        self.z_source = ifproc.line_redshift[bank] # z
        self.frest = ifproc.line_rest_frequency[bank] * np.float64(1.0e9) # Hz
        self.dop_track  = ifproc.doppler_track
        self.vbary = ifproc.barycenter_velocity # km/s
        self.vobs = ifproc.observatory_velocity # km/s
        self.fsky = ifproc.sky_frequency[bank] * np.float64(1.0e9) # Hz  
        self.flo1  = ifproc.lo_1_frequency * np.float64(1.0e9) # Hz  
        self.flo2  = ifproc.lo_2_frequency[bank] * np.float64(1.0e9) # Hz  
        self.fif1 = ifproc.if_1_frequency[bank] * np.float64(1.0e9) # Hz  
        self.fif2 = ifproc.if_2_frequency[bank] * np.float64(1.0e9) # Hz  
        self.fsyn = ifproc.synthesizer_frequency[bank] * np.float64(1.0e9) # Hz  
        self.foff = ifproc.frequency_offset[bank] * np.float64(1.0e9) # Hz  
        self.voff = -self.foff/self.frest * self.c # km/z
        self.sideband = ifproc.sideband[bank]

        self.nchan = nchan 
        self.ch0 = float(nchan-1)/2.
        self.bw    = bandwidth * np.float64(1.0e6) # Hz  
        if self.sideband == 1:
            self.dfdc = -self.bw/self.nchan
        else:
            self.dfdc = self.bw/self.nchan
        self.dvdc = -self.dfdc/self.frest * self.c # radio definition

        self.iarray = np.arange(self.nchan)

        # default is system velocity for the x array
        if self.v_system == 0:
            self.v0_lsr = self.v_source
            self.v0_bary = self.v_source + self.vbary
            self.v0_sky = self.vobs + self.v_source
            self.x_vlsr()
        elif self.v_system == 1:
            self.v0_lsr = self.v_source + self.vbary
            self.v0_bary = self.v_source
            self.v0_sky = self.vobs + self.v_source - self.vbary
            self.x_vbary()

    def x_vlsr(self):
        """ sets x axis to be in velocity units wrt the Local Standard of Rest
        """
        self.xv0 = self.v0_lsr + self.voff
        self.xarray = (self.iarray-self.ch0)*self.dvdc + self.xv0
        self.dxdc = self.xarray[1]-self.xarray[0]
        self.xtype = 0
        self.xname = 'VLSR'

    def x_vbary(self):
        """ sets x axis to be in velocity units wrt the Solar System Bary Center
        """
        self.xv0 = self.v0_bary + self.voff
        self.xarray = (self.iarray-self.ch0)*self.dvdc + self.xv0
        self.dxdc = self.xarray[1]-self.xarray[0]
        self.xtype = 0
        self.xname = 'VBARY'

    def x_vsky(self):
        """ sets x axis to be in velocity units wrt topocentric system
        """
        self.xv0 = self.v0_sky + self.voff
        self.xarray = (self.iarray-self.ch0)*self.dvdc + self.xv0
        self.dxdc = self.xarray[1]-self.xarray[0]
        self.xtype = 0  
        self.xname = 'VSKY'

    def x_vsrc(self):
        """ sets x axis to be in velocity units wrt the source velocity
        """
        self.xv0 = self.voff
        self.xarray = (self.iarray-self.ch0)*self.dvdc + self.xv0
        self.dxdc = self.xarray[1]-self.xarray[0]
        self.xtype = 0  
        self.xname = 'VSRCE'

    def x_fsky(self):
        """ sets x axis to be in frequency units accounting for topocentric doppler shift
        """
        self.xarray = (self.iarray-self.ch0)*self.dfdc + self.fsky
        self.dxdc = self.xarray[1]-self.xarray[0]
        self.xtype = 1  
        self.xname = 'FSKY'

    def x_flsr(self):
        """ sets x axis to be in frequency units accounting for Local Standard of Rest  doppler shift
        """
        self.xf0 = self.fsky + self.frest/self.c * (self.v0_sky-self.v0_lsr)
        self.xarray = (self.iarray-self.ch0)*self.dfdc + self.xf0
        self.dxdc = self.xarray[1]-self.xarray[0]
        self.xtype = 1  
        self.xname = 'FLSR'

    def x_fbary(self):
        """ sets x axis to be in frequency units accounting for S.S. Bary Center doppler shift
        """
        self.xf0 = self.fsky + self.frest/self.c * (self.v0_sky-self.v0_bary)
        self.xarray = (self.iarray-self.ch0)*self.dfdc + self.xf0
        self.dxdc = self.xarray[1]-self.xarray[0]
        self.xtype = 1  
        self.xname = 'FLSR'

    def x_fsrc(self):
        """ sets x axis to be in frequency units accounting for source doppler shift
        """
        self.xf0 = self.fsky + self.frest/self.c * self.v0_sky 
        self.xarray = (self.iarray-self.ch0)*self.dfdc + self.xf0
        self.dxdc = self.xarray[1]-self.xarray[0]
        self.xtype = 1
        self.xname = 'FSRCE'

    def v2c(self,v):
        """ converts velocity into channel number
        Args:
            v (float/np array): velocity for conversion
        Returns:
            c (int/np array): channels corresponding to velocity input
        """
        if self.xtype == 0:
            c = np.array(np.round((v-self.xv0)/self.dvdc + self.ch0),dtype=int)
        else:
            print('v2c: xarray is in frequency')
        return(c)

    def f2c(self,f):
        """ converts frequency into channel number
        Args:
            f (float/np array): frequency for conversion
        Returns:
            c (int/np array): channel numbers corresponding to frequency input
        """
        if self.xtype == 1:
            c = np.array(np.round((f-self.xf0)/self.dfdc + self.ch0),dtype=int)
        else:
            print('f2c: xarray is in velocity')
        return(c)

    def c2v(self,c):
        """ converts a channel number into a velocity                         
        Args:
            c (int): input channel number
        Returns:
            v (float): velocity of channel
        """
        if self.xtype == 0:
            v = (c-self.ch0)*self.dvdc + self.xv0
        else:
            print('c2v: xarray is in frequency')
        return(v)

    def c2f(self,c):
        """ converts a channel number into a frequency 
        Args:
            c (int): input channel number
        Returns:
            f (float): frequency of channel
        """
        if self.xtype == 1:
            f = (c-self.ch0)*self.dfdc + self.xf0
        else:
            print('c2f: xarray is in velocity')
        return(v)

    def make_velocity_list(self,regions):
        if self.xtype == 0:
            nregions = len(regions)
            channel_list = []
            for i in range(nregions):
                c0 = self.v2c(regions[i][0])
                c1 = self.v2c(regions[i][1])+1
                if c1>c0:
                    channel_list = channel_list + range(c0,c1)
                else:
                    channel_list = channel_list + range(c1,c0)
            nchannels = len(channel_list)
        else:
            print('make_velocity_list: xarray not in velocity units')
        return channel_list, nchannels

    def make_frequency_list(self,regions):
        if self.xtype == 1:
            nregions = len(regions)
            channel_list = []
            for i in range(nregions):
                c0 = self.f2c(regions[i][0])
                c1 = self.f2c(regions[i][1])+1
                if c1>c0:
                    channel_list = channel_list + range(c0,c1)
                else:
                    channel_list = channel_list + range(c1,c0)
            nchannels = len(channel_list)
        else:
            print('make_frequency_list: xarray not in frequency units')
        return channel_list, nchannels


    def make_channel_list(self,regions):
        """ creates a list of channels given a list of lists of intervals to be processed                         
            for example: regions = [[1,3],[5,8]] creates [1,2,3,5,6,7,8]                                          
            returns the list and the number of elements in the list                                               
        """
        nregions = len(regions)
        channel_list = []
        for i in range(nregions):
            channel_list = channel_list + range(regions[i][0],regions[i][1]+1)
        nchannels = len(channel_list)
        return channel_list, nchannels



class LineData(LineDataHeader):
    ''' provides methods to manipulate line data
    units of velocity: km/s
    units of frequency: Hertz
    '''
    def __init__(self,ifproc,bank,nchan,bandwidth,data):
        ''' constructor for LineData
        inputs:
        ifproc - an ifproc object with the ifproc information
        bank - selects which bank of the spectrometer is relevant 
        nchan - number of channels in the spectrum
        bandwidth - bandwidth of the spectrum (MHz)
        data - the spectral line data
        '''
        LineDataHeader.__init__(self,ifproc,bank,nchan,bandwidth)
        self.data  = data
        # yarray is the working spectrum - initialized as copy of data
        self.yarray = data.copy()

    def line(self):
        return Line(self.iarray, self.xarray, self.yarray, self.xname)

    def revert(self):
        """
        resets the spectrum to the original values
        """
        self.yarray = self.data.copy()

    def cgen(self, c0, c1, dc=1):
        """ generate a Line object with data covering range from channel c0 to channel c1
        Args:
            c0 (int): beginning channel number for new spectrum
            c1 (int): ending channel number for new spectrum
            dc (int): spacing of channels for new specrum
        Returns:
            result (Line object): Line object with resampled spectrum
        """
        new_iarray = np.arange(c0,c1,dc)
        fx = interp1d(self.iarray,self.xarray)
        fy = interp1d(self.iarray,self.yarray)
        new_xarray = fx(new_iarray)
        new_yarray = fy(new_iarray)
        result = Line(new_iarray,new_xarray,new_yarray,self.xname)
        return result

    def xgen(self, x0, x1, dx):
        """ generate a Line object with data covering range of x axis values from  x0 to x1
        Args:
            x0 (int): beginning x value for new spectrum
            x1 (int): ending x value for new spectrum
            dx (int): spacing of samples for new specrum
        Returns:
            result (Line object): Line object with resampled spectrum
        """
        new_xarray = np.arange(x0,x1,dx)
        fi = interp1d(self.xarray,self.iarray)
        fy = interp1d(self.xarray,self.yarray)
        new_iarray = fi(new_xarray)
        new_yarray = fy(new_xarray)
        result = Line(new_iarray,new_xarray,new_yarray,self.xname)
        return result

    def cslice(self, c0, c1):
        """ generate a Line object by taking a slice from original array
        Args:
            c0 (int): starting channel number of slice
            c1 (int): ending channel number of slice
        Returns:
            result (Line object): Line object with the slice
        """
        if c0<c1:
            new_iarray = np.arange(c0,c1)
            new_xarray = self.xarray[c0:c1]
            new_yarray = self.yarray[c0:c1]
        else:
            new_iarray = np.arange(c1,c0,-1)
            new_xarray = self.xarray[c1:c0]
            new_yarray = self.yarray[c1:c0]
        result = Line(new_iarray,new_xarray,new_yarray,self.xname)
        return result

    def vslice(self,v0,v1):
        """ generate a Line object by taking a slice from original array; slice obtained by giving velocity limits
        Args:
            v0 (int): starting velocity of slice
            v1 (int): ending velocity of slice
        Returns:
            result (Line object): Line object with the slice
        """
        c0 = self.v2c(v0)
        c1 = self.v2c(v1)
        if c0<c1:
            new_iarray = self.iarray[c0:c1]
            new_xarray = self.xarray[c0:c1]
            new_yarray = self.yarray[c0:c1]
        else:
            new_iarray = self.iarray[c1:c0]
            new_xarray = self.xarray[c1:c0]
            new_yarray = self.yarray[c1:c0]

        result = Line(new_iarray,new_xarray,new_yarray,self.xname)
        return result

    def fslice(self,f0,f1):
        """ generate a Line object by taking a slice from original array; slice obtained by giving frequency limits
        Args:
            f0 (int): starting freqeuncy of slice
            f1 (int): ending frequency of slice
        Returns:
            result (Line object): Line object with the slice
        """
        c0 = self.f2c(f0)
        c1 = self.f2c(f1)
        if c0<c1:
            new_iarray = self.iarray[c0:c1]
            new_xarray = self.xarray[c0:c1]
            new_yarray = self.yarray[c0:c1]
        else:
            new_iarray = self.iarray[c1:c0]
            new_xarray = self.xarray[c1:c0]
            new_yarray = self.yarray[c1:c0]

        result = Line(new_iarray,new_xarray,new_yarray,self.xname)
        return result

    
class Accum():
    """ manages spectral averaging
    
        All spectra to be averaged must be aligned in channel/index space
        Weighted averages are possible using weight provided when spectrum is added to the stack
    """
    def __init__(self):
        self.count = 0
        self.alist = []
        self.wlist = []

    def load(self,data,wt=1.):
        ndata = len(data)
        if self.count == 0:
            self.nchannels = ndata
            self.count += 1 
            self.alist.append(data)
            self.wlist.append(wt)
        elif self.count>0 and ndata!=self.nchannels:
            print('accum:load - WARNING - number of points in spectrum doesnt match')
        else:
            self.count += 1 
            self.alist.append(data)
            self.wlist.append(wt)

    def ave(self,type=0):
        self.average = np.zeros(self.nchannels)
        a = np.array(self.alist)
        w = np.array(self.wlist)
        for i in range(self.nchannels):
            slist = np.where(np.isfinite(a[:,i]))[0]
            self.average[i] = np.sum(a[slist,i]*w[slist])/np.sum(w[slist])


class NetCDFLineHeader():
    """ manages the header information to write to netcdf file
        assumes an open netcdf file : nc
    """
    def __init__(self,nc):
        self.nc = nc

    def write_line_data_header_variables(self,L):
        """
        writes creates line header variables for a net cdf file
        """
        t = self.nc.createVariable('Header.LineData.Bank','i4')
        self.nc.variables['Header.LineData.Bank'][0] = L.bank

        t = self.nc.createVariable('Header.LineData.VSource','f8')
        self.nc.variables['Header.LineData.VSource'][0] = L.v_source
        t.units = 'km/s'

        t = self.nc.createVariable('Header.LineData.VSystem','f8')
        self.nc.variables['Header.LineData.VSystem'][0] = L.v_system
        t.units = 'km/s'
        
        t = self.nc.createVariable('Header.LineData.ZSource','f8')
        self.nc.variables['Header.LineData.ZSource'][0] = L.z_source
        t.units = 'z'

        t = self.nc.createVariable('Header.LineData.LineRestFrequency','f8')
        self.nc.variables['Header.LineData.LineRestFrequency'][0] = L.frest
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.DopplerTrack','i4')
        self.nc.variables['Header.LineData.DopplerTrack'][0] = L.dop_track
        t.units = 'T/F'

        t = self.nc.createVariable('Header.LineData.VBarycenter','f8')
        self.nc.variables['Header.LineData.VBarycenter'][0] = L.vbary
        t.units = 'km/s'
 
        t = self.nc.createVariable('Header.LineData.VObservatory','f8')
        self.nc.variables['Header.LineData.VObservatory'][0] = L.vobs
        t.units = 'km/s'

        t = self.nc.createVariable('Header.LineData.SkyFrequency','f8')
        self.nc.variables['Header.LineData.SkyFrequency'][0] = L.fsky
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.LO1Frequency','f8')
        self.nc.variables['Header.LineData.LO1Frequency'][0] = L.flo1
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.LO2Frequency','f8')
        self.nc.variables['Header.LineData.LO2Frequency'][0] = L.flo2
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.IF1Frequency','f8')
        self.nc.variables['Header.LineData.IF1Frequency'][0] = L.fif1
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.IF2Frequency','f8')
        self.nc.variables['Header.LineData.IF2Frequency'][0] = L.fif2
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.SynthesizerFrequency','f8')
        self.nc.variables['Header.LineData.SynthesizerFrequency'][0] = L.fsyn
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.FrequencyOffset','f8')
        self.nc.variables['Header.LineData.FrequencyOffset'][0] = L.foff
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.VelocityOffset','f8')
        self.nc.variables['Header.LineData.VelocityOffset'][0] = L.voff
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.Sideband','i4')
        self.nc.variables['Header.LineData.Sideband'][0] = L.sideband
        t.units = 'USB/LSB'

        t = self.nc.createVariable('Header.LineData.NChannels','i4')
        self.nc.variables['Header.LineData.NChannels'][0] = L.nchan
        
        t = self.nc.createVariable('Header.LineData.Bandwidth','f8')
        self.nc.variables['Header.LineData.Bandwidth'][0] = L.bw
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.DeltaFrequency','f8')
        self.nc.variables['Header.LineData.DeltaFrequency'][0] = L.dfdc
        t.units = 'Hz'

        t = self.nc.createVariable('Header.LineData.DeltaVelocity','f8')
        self.nc.variables['Header.LineData.DeltaVelocity'][0] = L.dvdc
        t.units = 'km/s'

        t = self.nc.createVariable('Header.LineData.LSRFrameV0','f8')
        self.nc.variables['Header.LineData.LSRFrameV0'][0] = L.v0_lsr
        t.units = 'km/s'

        t = self.nc.createVariable('Header.LineData.BarycenterFrameV0','f8')
        self.nc.variables['Header.LineData.BarycenterFrameV0'][0] = L.v0_bary
        t.units = 'km/s'

        t = self.nc.createVariable('Header.LineData.SkyFrameV0','f8')
        self.nc.variables['Header.LineData.SkyFrameV0'][0] = L.v0_sky
        t.units = 'km/s'

        t = self.nc.createVariable('Header.LineData.XType','i4')
        self.nc.variables['Header.LineData.XType'][0] = L.xtype
        t.long_name = 'Axis Type'

        iarray = self.nc.createVariable('Header.LineData.ChannelNumber','i4',('nchan',))
        iarray[0:L.nchan] = L.iarray[0:L.nchan]
        iarray.long_name = 'Channel Number'

        # the following are the generic spectrum axis parameters

        t = self.nc.createVariable('Header.SpectrumAxis.CRPIX','f8')
        self.nc.variables['Header.SpectrumAxis.CRPIX'][0] = L.ch0
        t.long_name = 'CRPIX'

        t = self.nc.createVariable('Header.SpectrumAxis.CRVAL','f8')
        self.nc.variables['Header.SpectrumAxis.CRVAL'][0] = L.xv0
        t.long_name = 'CRVAL'
        if L.xtype == 0:
            t.units = 'km/s'
        else:
            t.units = 'Hz'

        t = self.nc.createVariable('Header.SpectrumAxis.CDELT','f8')
        self.nc.variables['Header.SpectrumAxis.CDELT'][0] = L.dxdc
        t.long_name = 'CDELT'
        if L.xtype == 0:
            t.units = 'km/s'
        else:
            t.units = 'Hz'

        ctype = self.nc.createVariable('Header.SpectrumAxis.CTYPE','c',('nlabel',))
        ctype[0:len(L.xname)] = L.xname[0:len(L.xname)]
        ctype.long_name = 'CTYPE'

        xarray = self.nc.createVariable('Header.SpectrumAxis.CAXIS','f8',('nchan',))
        xarray[0:L.nchan] = L.xarray[0:L.nchan]
        if L.xtype == 0:
            xarray.units = 'km/s'
        else:
            xarray.units = 'Hz'
        xarray.long_name = L.xname

    def read_line_data_header_variables(self,L):
        """
        reads the header data from a net cdf file
        """
        L.bank = self.nc.variables['Header.LineData.Bank'][0]
        L.v_source = self.nc.variables['Header.LineData.VSource'][0]
        L.v_system = self.nc.variables['Header.LineData.VSystem'][0]
        L.z_source = self.nc.variables['Header.LineData.ZSource'][0]
        L.frest = self.nc.variables['Header.LineData.LineRestFrequency'][0]
        L.dop_track = self.nc.variables['Header.LineData.DopplerTrack'][0]
        L.vbary = self.nc.variables['Header.LineData.VBarycenter'][0]
        L.vobs = self.nc.variables['Header.LineData.VObservatory'][0]
        L.fsky = self.nc.variables['Header.LineData.SkyFrequency'][0]
        L.flo1 = self.nc.variables['Header.LineData.LO1Frequency'][0]
        L.flo2 = self.nc.variables['Header.LineData.LO2Frequency'][0]
        L.fif1 = self.nc.variables['Header.LineData.IF1Frequency'][0]
        L.fif2 = self.nc.variables['Header.LineData.IF2Frequency'][0]
        L.fsyn = self.nc.variables['Header.LineData.SynthesizerFrequency'][0]
        L.foff = self.nc.variables['Header.LineData.FrequencyOffset'][0]
        L.voff = self.nc.variables['Header.LineData.VelocityOffset'][0]
        L.sideband = self.nc.variables['Header.LineData.Sideband'][0]
        L.nchan = self.nc.variables['Header.LineData.NChannels'][0]
        L.bw = self.nc.variables['Header.LineData.Bandwidth'][0]
        L.dfdc = self.nc.variables['Header.LineData.DeltaFrequency'][0]
        L.dvdc = self.nc.variables['Header.LineData.DeltaVelocity'][0]
        L.v0_lsr = self.nc.variables['Header.LineData.LSRFrameV0'][0]
        L.v0_bary = self.nc.variables['Header.LineData.BarycenterFrameV0'][0]
        L.v0_sky = self.nc.variables['Header.LineData.SkyFrameV0'][0]
        L.xtype = self.nc.variables['Header.LineData.XType'][0]
        L.iarray = self.nc.variables['Header.LineData.ChannelNumber'][:]
        # the following are the generic spectrum axis parameters
        L.ch0 = self.nc.variables['Header.SpectrumAxis.CRPIX'][0]
        L.xv0 = self.nc.variables['Header.SpectrumAxis.CRVAL'][0]
        L.dxdc = self.nc.variables['Header.SpectrumAxis.CDELT'][0]
        L.xname = netCDF4.chartostring(self.nc.variables['Header.SpectrumAxis.CTYPE'][:])
        L.xarray = self.nc.variables['Header.SpectrumAxis.CAXIS'][:]


    def write_line_header_variables(self,L):
        t = self.nc.createVariable('Header.Line.NChannels','i4')
        self.nc.variables['Header.Line.NChannels'][0] = L.nchan

        iarray = self.nc.createVariable('Header.Line.ChannelNumber','i4',('nchan',))
        iarray[0:L.nchan] = L.iarray[0:L.nchan]
        iarray.long_name = 'Channel Number'
        
        # the following are the generic spectrum axis parameters
        t = self.nc.createVariable('Header.SpectrumAxis.CDELT','f8')
        self.nc.variables['Header.SpectrumAxis.CDELT'][0] = L.dxdc
        
        t = self.nc.createVariable('Header.SpectrumAxis.CRPIX','f8')
        self.nc.variables['Header.SpectrumAxis.CRPIX'][0] = 0.0

        t = self.nc.createVariable('Header.SpectrumAxis.CRVAL','f8')
        self.nc.variables['Header.SpectrumAxis.CRVAL'][0] = L.xarray[0]

        ctype = self.nc.createVariable('Header.SpectrumAxis.CTYPE','c',('nlabel',))
        ctype[0:len(L.xname)] = L.xname[0:len(L.xname)]

        xarray = self.nc.createVariable('Header.SpectrumAxis.CAXIS','f8',('nchan',))
        xarray[0:L.nchan] = L.xarray[0:L.nchan]

    def read_line_header_variables(self,L):
        L.nchan = self.nc.variables['Header.Line.NChannels'][0]
        L.iarray = self.nc.variables['Header.Line.ChannelNumber'][:]
        # the following are the generic spectrum axis parameters
        L.dxdc = self.nc.variables['Header.SpectrumAxis.CDELT'][0]
        L.xname = netCDF4.chartostring(self.nc.variables['Header.SpectrumAxis.CTYPE'][:])
        L.xarray = self.nc.variables['Header.SpectrumAxis.CAXIS'][:]


