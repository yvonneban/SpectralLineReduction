""" Module to read a set of spectra from the spectrometer

classes: spec_bank, roach_spec, spec_bank_data, spec_bank_cal
methods: lookup_roach_files, find_roach_from_pixel, create_roach_list
uses: ifproc, grid
author: FPS
date: May 2018
changes:
python 3
"""
import numpy as np
import math
import os
import sys
import fnmatch
import netCDF4
from scipy import interpolate
from operator import itemgetter
from itertools import groupby
from lmtslr.grid.grid import Grid

from lmtslr.ifproc.ifproc import IFProc

# define all the pixels in the roach boards they appear in
roach_pixels_all = [[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]]



class RoachSpec():
    """ this BASE class deals with a single time series of spectra """
    def __init__(self, obsnum, roach_id, roach_input,
                 nchan, bandwidth, nspec, raw_spec, spec_time,
                 xmap, ymap, pmap, bufpos):
        self.obsnum = obsnum
        self.roach_id = roach_id
        self.roach_input = roach_input
        self.nchan = nchan
        self.bandwidth = bandwidth
        self.nspec = nspec
        self.raw_spec = raw_spec
        self.spec_time = spec_time
        self.xmap = xmap
        self.ymap = ymap
        self.pmap = pmap
        self.bufpos = bufpos

        self.refs,self.ref_ranges,self.nrefs = self.get_ranges(1)
        self.ons,self.on_ranges,self.nons = self.get_ranges(0)
        self.hots,self.hot_ranges,self.nhots = self.get_ranges(3)
        self.skys,self.sky_ranges,self.nskys = self.get_ranges(2)

    def get_ranges(self, value=1):
        """ searches bufpos and finds indices with values corresponding to value
            value is type of observation for search: 0 is ons; 1 is refs; 2 is sky; 3 is hot
        """
        if value == 0:
            idx = np.where((self.bufpos == value) | (self.bufpos >= 100))[0]
        else:
            idx = np.where(self.bufpos == value)[0]
        ranges = []
        for k, g in groupby(enumerate(idx), lambda x:x[0]-x[1]):  
            group = list(map(itemgetter(1), g))
            ranges.append((group[0], group[-1]))
        return(idx,ranges,len(ranges))


    def compute_reference_spectra(self):
        """ computes the average reference spectrum for each block of refs
            these are returned as a list
        """
        self.reference_spectra = np.zeros((self.nrefs,self.nchan))
        for iref in range(self.nrefs):
            the_refs = range(self.ref_ranges[iref][0],self.ref_ranges[iref][1]+1)
            self.reference_spectra[iref,:] = np.mean(self.raw_spec[the_refs,:],axis=0)

    def compute_main_spectra(self):
        """ computes the average main (on) spectrum for each block of on-source observations
            these are returned as a list
        """
        self.main_spectra = np.zeros((self.nons,self.nchan))
        for ion in range(self.nons):
            the_ons = range(self.on_ranges[ion][0],self.on_ranges[ion][1]+1)
            self.main_spectra[ion,:] = np.mean(self.raw_spec[the_ons,:],axis=0)

    def compute_median_spectrum(self):
        """ computes the median spectrum from ALL spectra in the file
        """
        self.median_spectrum = np.median(self.raw_spec[:,:],axis=0)

    def compute_reference_spectrum(self):
        """ computes the average reference spectrum using all ref observations
        """
        self.reference_spectrum = np.mean(self.raw_spec[self.refs,:],axis=0)
        
    def compute_main_spectrum(self):
        """ computes the average main spectrum using all main (on) observations
        """
        self.main_spectrum = np.mean(self.raw_spec[self.ons,:],axis=0)
        
    def reduce_on_spectrum(self, calibrate=False, tsys_spectrum=0,
                           tsys_no_cal=1):
        """ creates a ON spectrum returned as self.on_spectrum. Reduction procedure depends
        on the type parameter:
        calibrate is True when we want to calibrate the ps_spectrum
        tsys_spectrum is a spectrum of tsys values for the calibration used when calibrate=True
        tsys_no_cal is tsys value for the calibration used when calibrate=False
        """
        self.compute_main_spectrum()
        self.on_spectrum = self.main_spectrum
        # calibrate it if requested
        if calibrate == True:
            self.on_spectrum = self.on_spectrum * tsys_spectrum
        else:
            self.on_spectrum = self.on_spectrum * tsys_no_cal

    def reduce_ps_spectrum(self, type=2, normal_ps=True,
                           calibrate=False, tsys_spectrum=0, tsys_no_cal=1):
        """ creates a PS spectrum returned as self.ps_spectrum. Reduction procedure depends
        on the type parameter:
        # type 0 is not defined
        # type 1: compute average of all refs and mains and take the difference
        # type 2: compute individual on's and off's in blocks; difference each pair and average differences
        normal_ps is True when we want (Main-Ref)/Ref.  Set False for (Ref-Main)/Main for beam switch
        calibrate is True when we want to calibrate the ps_spectrum
        tsys_spectrum is a spectrum of tsys values for the calibration used when calibrate=True
        tsys_no_cal is tsys value for the calibration used when calibrate=False
        """
        if type == 1:
            self.compute_main_spectrum()
            self.compute_reference_spectrum()
            if normal_ps == True:
                self.ps_spectrum = (self.main_spectrum-self.reference_spectrum)/self.reference_spectrum
            else:
                self.ps_spectrum = (self.reference_spectrum-self.main_spectrum)/self.main_spectrum
        else: # type == 2:
            self.compute_main_spectra()
            self.compute_reference_spectra()
            if self.nons == self.nrefs:
                ps_list = np.zeros((self.nons,self.nchan))
                if normal_ps == True:
                    for i in range(self.nons):
                        ps_list[i,:] = (self.main_spectra[i,:]-self.reference_spectra[i,:])/self.reference_spectra[i,:]
                else:
                    for i in range(self.nons):
                        ps_list[i,:] = (self.reference_spectra[i,:]-self.main_spectra[i,:])/self.main_spectra[i,:]
                self.ps_spectrum = np.mean(ps_list[:,:],axis=0)
            else:
                print('check number of ons and refs %d %d'%(self.nons,self.nrefs))
            
        # calibrate it if requested
        if calibrate == True:
            self.ps_spectrum = self.ps_spectrum * tsys_spectrum
        else:
            self.ps_spectrum = self.ps_spectrum * tsys_no_cal

    def reduce_spectra(self, type=0, calibrate=False,
                       tsys_spectrum=0, tsys_no_cal=1):
        """ creates a list of all the "on" spectra in the file.  The on spectra are
            reduced according to the value of "type":
            # type 0: use the median spectrum for the whole thing
            # type 1: use a single reference spectra which is average of all refs
            # type 2: use the average of refs which bracket the ons
        """
        spectra = []

        if self.nrefs == 0:
            type = 0

        if type == 0:
            self.compute_median_spectrum()
            for i in self.ons:
                spectra.append( (self.raw_spec[i,:]-self.median_spectrum[:]) / self.median_spectrum[:] )
        elif type == 1:
            if self.nrefs != 0:
                self.compute_reference_spectrum()
                for i in self.ons:
                    spectra.append( (self.raw_spec[i,:]-self.reference_spectrum[:]) / self.reference_spectrum[:] )
            else:
                for i in self.ons:
                    spectra.append( (self.raw_spec[i,:]) )
        else: # type == 2:
            if self.nrefs != 0:
                self.compute_reference_spectra()
                nbins = self.nrefs-1
                for ibin in range(nbins):
                    istart = self.on_ranges[ibin][0]
                    istop = self.on_ranges[ibin][1]
                    for i in range(istart,istop+1):
                        ref = (self.reference_spectra[ibin]+self.reference_spectra[ibin+1])/2.
                        spectra.append( (self.raw_spec[i,:]-ref) / ref)
            else:
                for i in self.ons:
                    spectra.append( (self.raw_spec[i,:]) )
                        
        # save reduced spectra as a 2D numpy array
        # calibrate it if requested 
        if calibrate == False:
            self.reduced_spectra = np.array(spectra) * tsys_no_cal
        else:
            self.reduced_spectra = np.array(spectra) * tsys_spectrum

    def baseline(self, spectrum, baseline_list, n_baseline_list,
                 baseline_order=0):
        """ computes a baseline given a spectrum using baseline_list of channels
        """
        baseline_list = [ i for i in baseline_list if i < len(spectrum) ]
        params = np.polyfit(baseline_list,spectrum[baseline_list],baseline_order)
        x = np.arange(0,self.nchan)
        baseline = np.polyval(params,x)
        residuals = spectrum[baseline_list] - np.polyval(params,baseline_list)
        rms = np.sqrt(np.dot(residuals.transpose(),residuals)/n_baseline_list)
        return(baseline,rms)

    def integrate_spectra(self, channel_list, n_channel_list,
                          baseline_list, n_baseline_list, baseline_order=0, type=0):
        """ computes the integral of the spectrum over some range of channels (channel_list)
            a baseline is removed using the baseline_list
            type =0 fpr YINT  =1 for YMAX
        """
        spectra = []
        for i,isample in enumerate(self.ons):
            baseline,rms = self.baseline(self.reduced_spectra[i],baseline_list,n_baseline_list,baseline_order)
            baselined_spectrum = self.reduced_spectra[i]-baseline
            if type == 0:
                spectra.append(np.sum(baselined_spectrum[channel_list]))
            else:
                spectra.append(np.max(baselined_spectrum[channel_list]))

        return(np.array(spectra))

    def integrate_spectrum(self, on_list, channel_list, n_channel_list,
                           baseline_list, n_baseline_list, baseline_order=0, type=0):
        """ averages reduced spectra over list of indices (on_list) and then
            computes the integral of the spectrum over some range of channels (channel_list)
            a baseline is removed using the baseline_list.
            type =0 fpr YINT  =1 for YMAX
        """
        spectrum = np.mean(self.reduced_spectra[on_list,:],axis=1)[0]
        #print(on_list)
        #print(np.mean(self.xmap[on_list]),np.mean(self.ymap[on_list]))
        #print(spectrum)
        baseline,rms = self.baseline(spectrum,baseline_list,n_baseline_list,baseline_order)
        baselined_spectrum = spectrum - baseline 
        if type == 0:
            result = np.sum(baselined_spectrum[channel_list]) # YINT
        else:
            result = np.max(baselined_spectrum[channel_list]) # YMAX
            
        return(baselined_spectrum,result)

    def compute_tsys_spectrum(self, bdrop=100, edrop=100):
        """ computes a system temperature spectrum and average value of tsys
            set bdrop to be number of channels excluded from beginning of spectrum
            set edrop to be number of channels excluded from end of spectrum
        """
        if ((self.nhots>0) and (self.nskys>0)):
            hot_spectrum = np.median(self.raw_spec[self.hots,:],axis=0)
            sky_spectrum = np.median(self.raw_spec[self.skys,:],axis=0)
            self.tsys_spectrum = 280.*sky_spectrum/(hot_spectrum-sky_spectrum)
            # find the index where tsys_spectrum in finite
            indx_fin = np.where(np.isfinite(self.tsys_spectrum))
            # compute tsys as the mean of finite tsys_spectrum
            self.tsys = np.mean(self.tsys_spectrum[indx_fin][bdrop:self.nchan-edrop])
            # find the index where tsys_spectrum in infinite
            indx_inf = np.where(np.isinf(self.tsys_spectrum))
            # replace infinite tsys_spectrum with the mean
            self.tsys_spectrum[indx_inf] = self.tsys
        else:
            print('ObsNum %d Roach %d does not have calibration data'%(self.obsnum,self.roach_id))
            self.tsys_spectrum = np.zeros(self.nchan)
            self.tsys = 0.

    class LineStatistics():
        def __init__(self, parent, v, spectrum, channel_list, n_channel_list,
                     baseline_list, n_baseline_list, baseline_order=0, pixel_id_0=0,
                     pixel_id_1=1, obspgm=None):
            """ removes a baseline and computes line statistics for a spectrum
            Line Statistics are:
            YMAX = maximum value in channel_list
            CMAX = channel of maximum value
            XMAX = velocity of maximum value
            YINT = integral of line over channel list in units of K km/s
            YERR = error on YINT
            XMEAN = first moment of line integrated over channel list
            XWIDTH = estimate of line width based on Peak value and integrated intensity
            RMS = rms of baseline fit
            """
            self.pixel_id_0 = pixel_id_0
            self.pixel_id_1 = pixel_id_1
            self.obspgm = obspgm
            delta_v = np.abs(v[1]-v[0])
            baseline,rms = parent.baseline(spectrum,baseline_list,n_baseline_list,baseline_order)
            baselined_spectrum = spectrum - baseline
            self.v = v
            self.spectrum = baselined_spectrum
            self.rms = rms
            self.ymax = np.max(baselined_spectrum[channel_list])
            x = np.arange(0,parent.nchan)
            self.cmax = x[channel_list[np.where(baselined_spectrum[channel_list] == self.ymax)[0][0]]]
            self.xmax = v[self.cmax]
            self.yint = np.sum(baselined_spectrum[channel_list])*delta_v
            self.yerr = delta_v*self.rms*np.sqrt(n_channel_list)
            self.xmean = np.sum(baselined_spectrum[channel_list]*v[channel_list])/n_channel_list
            self.xwidth = self.yint/self.ymax

        def to_string(self):
            str = '%s Pix=%02d'%(self.obspgm,self.pixel_id_0)
            if self.pixel_id_0 != self.pixel_id_1:
                str += '/%02d'%self.pixel_id_1
            str += ' ymax=%.3f cmax=%d xmax=%.3f yint=%.3f yerr=%.3f xmean=%.3f xwidth=%.3f rms=%.3f'%(self.ymax,self.cmax,self.xmax,self.yint,self.yerr,self.xmean,self.xwidth,self.rms)
            return str


class SpecBank():
    """ this BASE class deals with a complete "bank" of spectra """
    def __init__(self,roach_files,ifproc_data,pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],time_offset=[-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03],bank=0):
        """ inputs: 
        roach files  is list of files for the roach boards (nominally 4 of them)
        ifproc_data  is data from the corresponding ifproc file
        pixel_list   is a list of pixels to process
        bank         is the bank number for the receiver.
        time_offset  is the time lag for each roach board (s) 
        """
        #self.ifproc = ifproc_data(ifproc_file)
        self.ifproc = ifproc_data
        # most header variables in ifproc object; copy these two since often used
        self.obsnum = self.ifproc.obsnum
        self.calobsnum = self.ifproc.calobsnum
        self.elev = self.ifproc.elev
        self.obspgm = self.ifproc.obspgm
        self.source = self.ifproc.source
        self.map_coord = self.ifproc.map_coord

        # timing offsets for each roach board
        self.time_offset = time_offset

        self.x_interpolation_function = interpolate.interp1d(self.ifproc.time,self.ifproc.azmap,bounds_error=False)
        self.y_interpolation_function = interpolate.interp1d(self.ifproc.time,self.ifproc.elmap,bounds_error=False)
        self.p_interpolation_function = interpolate.interp1d(self.ifproc.time,self.ifproc.parang,bounds_error=False)
        self.b_interpolation_function = interpolate.interp1d(self.ifproc.time,self.ifproc.bufpos,kind='nearest',bounds_error=False)

        self.gaps = np.where(self.ifproc.bufpos[:-1] != self.ifproc.bufpos[1:])[0]
        self.time_gaps_0 = np.append([self.ifproc.time[0]],self.ifproc.time[self.gaps+1])
        self.time_gaps_1 = np.append(self.ifproc.time[self.gaps],[self.ifproc.time[-1]])

        self.nroach = len(roach_files)

        self.roach_pixel_ids = []
        self.roach = []
        for i in range(self.nroach):
            self.roach_pixel_ids = self.roach_pixel_ids + self.read_roach(roach_files[i],i,pixel_list)
        self.roach_pixel_ids = np.array(self.roach_pixel_ids)
        self.npix = len(self.roach)
        self.nchan = self.roach[0].nchan
        self.bandwidth = self.roach[0].bandwidth
        self.channel_0 = (self.nchan-1) / 2.0
        self.velocity_0 = self.ifproc.velocity
        self.line_rest_frequency = self.ifproc.line_rest_frequency[bank]
        self.receiver = self.ifproc.receiver
        self.sideband = self.ifproc.sideband[bank]
        self.tracking_beam = self.ifproc.tracking_beam

        # dfdc sign depends on the sideband (USB:- and LSB+)
        if self.sideband == 1:
            self.dfdc = -self.bandwidth*np.float64(1.0e6)/self.nchan
        else:
            self.dfdc = +self.bandwidth*np.float64(1.0e6)/self.nchan
        self.dvdc = -self.dfdc/(self.line_rest_frequency*np.float64(1.0e9))*np.float64(2.99792458e5)

        self.cal_flag = False
        
        # define the windows for line integrations and baseline fits
        self.nc = 0
        self.clist = []
        self.baseline_order = 0
        self.nb = 0
        self.blist = []

        
    def make_pixel_id_list(self, roach_files, pixel_list):
        """ creates a list of pixel id numbers corresponding to list of roach spectra that is created 
        """
        l = []
        for f in roach_files:
            for i in range(4):
                roach_name = 'roach%d'%i
                if roach_name in f:
                    for p in roach_pixels_all[i]:
                        if p in pixel_list:
                            l = l + [p]
                    break
        return(np.array(l))

    def make_channel_list(self,regions):
        """ creates a list of channels given a list of lists of intervals to be processed
            for example: regions = [[1,3],[5,8]] creates [1,2,3,5,6,7,8]
            returns the list and the number of elements in the list
        """
        nregions = len(regions)
        channel_list = []
        for i in range(nregions):
            # python 3 requires range to be converted to list
            channel_list = channel_list + list(range(regions[i][0],regions[i][1]+1))
        nchannels = len(channel_list)
        return(channel_list,nchannels)

    def make_velocity_list(self,velocity_regions,id='line'):
        """ creates a channel list from a set of velocity intervals
        """
        nregions = len(velocity_regions)
        channel_list = []
        for i in range(nregions):
            c0 = self.v2c(velocity_regions[i][0])
            c1 = self.v2c(velocity_regions[i][1])+1
            # python 3 requires range to be converted to list
            if c1>c0:
                channel_list = channel_list + list(range(c0,c1))
            else:
                channel_list = channel_list + list(range(c1,c0))
        nchannels = len(channel_list)
        if id == 'line':
            self.clist = channel_list
            self.nc = nchannels
        elif id == 'baseline':
            self.blist = channel_list
            self.nb = nchannels
        return(channel_list,nchannels)

    def v2c(self,v):
        """ converts velocity into channel number
        """
        c = np.array(np.round((v-self.velocity_0)/self.dvdc + self.channel_0),dtype=int)
        return(c)
    
    def c2v(self,c):
        """ converts a channel number into a velocity
        """
        v = (c-self.channel_0)*self.dvdc + self.velocity_0
        return(v)

    def create_velocity_scale(self):
        """ computes the velocity scale corresponding to the channels in the spectrum
        """
        return(self.c2v(range(self.nchan)))

    def find_pixel_index(self,ipixel):
        """ looks up the index in our array of spectra corresponding to a specific pixel
        """
        index = np.where(self.roach_pixel_ids == ipixel)
        return (int(index[0])) if len(index[0]) else 0

    def find_map_pixel_index(self,ipixel):
        """ looks up the index in our array of spectra corresponding to a specific pixel
        """
        index = np.where(self.map_pixel_list == ipixel)
        return (int(index[0])) if len(index[0]) else 0

    def read_roach(self,filename,roach_id,pixel_list):
        """ reads a roach file and creates spectrum (roach class) for each input
            filename is the name of the roach file
            roach_id is the number of the roach.  e.g. roach0 would be 0
            pixel_list is the list of pixels we want to process.
        """
        pixel_ids = []
        if os.path.isfile(filename):
            nc = netCDF4.Dataset(filename)
            print('read_roach', filename)
            # header information
            obsnum = nc.variables['Header.Telescope.ObsNum'][0]
            nchan = nc.variables['Header.Mode.numchannels'][0]
            bandwidth = nc.variables['Header.Mode.Bandwidth'][0]
            ninputs = 4

            datatime = nc.variables['Data.Integrate.time'][:]
            rawdata = nc.variables['Data.Integrate.Data'][:]
            input_list = nc.variables['Data.Integrate.Inputs'][:]

            acc_n = nc.variables['Data.Integrate.acc_n'][:]
            sync_time = nc.variables['Data.Integrate.sync_time'][:]
            read_time = nc.variables['Data.Integrate.read_time'][:]
            datatime = datatime - read_time
            
            # get roach index
            roach_index = roach_id
            for i in range(4):
                roach_name = 'roach%d'%i
                if roach_name in filename:
                    roach_index = i
                    break

            for input_chan in range(ninputs):
                # check to see whether this input matches one of the pixels we would like
                # if no match, then don't process it into the list
                if roach_pixels_all[roach_index][input_chan] not in pixel_list:
                    continue
                # if there is a match, then add the pixel id to our list and process
                pixel_ids = pixel_ids + [roach_pixels_all[roach_index][input_chan]]
                ilist  = np.where(input_list == input_chan)
                ilist0  = np.where(input_list == 0)
                raw_spec = rawdata[ilist,:][0]

                print('r:%d inp:%d pix:%d to:%f'%(roach_index,input_chan,roach_pixels_all[roach_index][input_chan],self.time_offset[roach_pixels_all[roach_index][input_chan]]))

                spec_time = datatime[ilist0] + self.time_offset[roach_pixels_all[roach_index][input_chan]]
                nspec = len(spec_time)
                # use the interpolation functions to find map positions,paralactic angle, and bufpos
                # in python 3 we must fix spec_time which was a masked array
                xmap = self.x_interpolation_function(np.ma.getdata(spec_time,subok=False))
                ymap = self.y_interpolation_function(np.ma.getdata(spec_time,subok=False))
                pmap = self.p_interpolation_function(np.ma.getdata(spec_time,subok=False))
                bufpos = self.b_interpolation_function(np.ma.getdata(spec_time,subok=False))

                # correct the interpolated arrays to remove points not in range for interpolation
                l = len(self.time_gaps_1)
                cond = False
                for i in range(l):
                    cond = cond | ((spec_time > self.time_gaps_0[i]) & (spec_time < self.time_gaps_1[i]))
                spec_time = spec_time[cond]
                xmap = xmap[cond]
                ymap = ymap[cond]
                pmap = pmap[cond]
                bufpos = bufpos[cond].astype(int)
                raw_spec = raw_spec[cond]
                
                nspec = len(spec_time)
                indx_small, = np.where(raw_spec[:,1024] < 100)
                if indx_small.any() > 0:
                    print('low_power',obsnum,roach_index,input_chan,indx_small)
                # now we append the individual roach spectrum object to our list
                self.roach.append(RoachSpec(obsnum,roach_id,input_chan,nchan,bandwidth,nspec,raw_spec,spec_time,xmap,ymap,pmap,bufpos))
            nc.close()
        else:
            print('%s does not exist for roach_id=%d'%(self.filename,roach_id))
        return pixel_ids                            

class SpecBankData(SpecBank):
    """ this BASE class deals with a complete "bank" of spectra """
    def __init__(self,roach_files,
                 ifproc_data,
                 pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],
                 time_offset=[-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03,-0.03],bank=0):
        """ inputs: 
        roach files  is list of files for the roach boards (nominally 4 of them)
        ifproc_data  is data from the corresponding ifproc file
        pixel_list   is a list of pixels to process
        time_offset  is a list of the time offsets for the pixels
        """
        SpecBank.__init__(self,roach_files,ifproc_data,pixel_list=pixel_list,time_offset=time_offset,bank=bank)
    
    def create_map_data(self,channel_list,n_channel_list,baseline_list,n_baseline_list,baseline_order=0,pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],type=0):
        """ processes a list of pixel ids to make a set of integrated spectra for mapping (and fitting)
        """
        t_list = []
        data_list = []
        x_list = []
        y_list = []
        p_list = []
        n_list = []
        mp_list = []
        for ipix in pixel_list:
            i = self.find_pixel_index(ipix)
            mp_list.append(ipix)
            t_list.append(self.roach[i].spec_time[self.roach[i].ons])
            x_list.append(self.roach[i].xmap[self.roach[i].ons])
            y_list.append(self.roach[i].ymap[self.roach[i].ons])
            p_list.append(self.roach[i].pmap[self.roach[i].ons])
            n_list.append(len(self.roach[i].xmap[self.roach[i].ons]))
            data_list.append(self.roach[i].integrate_spectra(channel_list,n_channel_list,baseline_list,n_baseline_list,baseline_order,type=type))
        self.map_pixel_list = np.array(mp_list)
        self.map_t = np.array(t_list)
        self.map_x = np.array(x_list)
        self.map_y = np.array(y_list)
        self.map_p = np.array(p_list)
        self.map_n = np.array(n_list)
        self.map_bufpos = np.array(n_list)
        self.map_data = np.array(data_list)

    def create_map_grid_bufpos_from_grid(self,xgrid,ygrid,tole,channel_list,n_channel_list,baseline_list,n_baseline_list,baseline_order=0,pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],type=0):
        ngridx = len(xgrid)
        ngridy = len(ygrid)

        bufpos = 100
        for ipix in pixel_list:
            i = self.find_pixel_index(ipix)
            for iypt in range(ngridy):
                ypt = ygrid[iypt]
                for ixpt in range(ngridx):
                    xpt = xgrid[ixpt]
                    #grid_list = np.where( np.sqrt( (self.roach[i].xmap[self.roach[i].ons]-xpt)**2+(self.roach[i].ymap[self.roach[i].ons]-ypt)**2) < tole)
                    grid_list = np.where( np.sqrt( (self.roach[i].xmap-xpt)**2+(self.roach[i].ymap-ypt)**2) < tole)
                    self.roach[i].bufpos[grid_list] = bufpos
                    #print(xpt,ypt,bufpos)
                    bufpos = bufpos + 1
        #np.set_printoptions(threshold=np.nan)
        #print(self.roach[i].bufpos)
                    
    def create_map_grid_data_from_grid(self,xgrid,ygrid,tole,channel_list,n_channel_list,baseline_list,n_baseline_list,baseline_order=0,pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],type=0):
        ngridx = len(xgrid)
        ngridy = len(ygrid)

        t_list = []
        data_list = []
        test_x_list = []
        test_y_list = []
        bufpos_list = []
        x_list = []
        y_list = []
        n_list = []
        mp_list = []
        spectrum_list = []

        for ipix in pixel_list:
            i = self.find_pixel_index(ipix)
            mp_list.append(ipix)
            pt_list = []
            pn_list = []
            pbufpos_list = []
            px_list = []
            py_list = []
            ptest_x_list = []
            ptest_y_list = []
            pdata_list = []
            pspectrum_list = []
            for iypt in range(ngridy):
                ypt = ygrid[iypt]
                for ixpt in range(ngridx):
                    xpt = xgrid[ixpt]
                    #print(xpt,ypt)
                    grid_list = np.where( np.sqrt( (self.roach[i].xmap[self.roach[i].ons]-xpt)**2+(self.roach[i].ymap[self.roach[i].ons]-ypt)**2) < tole)
                    pt_list.append(np.mean(self.roach[i].spec_time[self.roach[i].ons][grid_list]))
                    pn_list.append(len(self.roach[i].xmap[self.roach[i].ons][grid_list]))
                    pbufpos_list.append(np.mean(self.roach[i].bufpos[self.roach[i].ons][grid_list]))
                    px_list.append(np.mean(self.roach[i].xmap[self.roach[i].ons][grid_list]))
                    py_list.append(np.mean(self.roach[i].ymap[self.roach[i].ons][grid_list]))

                    ptest_x_list.append(self.roach[i].xmap[self.roach[i].ons][grid_list])
                    ptest_y_list.append(self.roach[i].ymap[self.roach[i].ons][grid_list])


                    # z is the integrated intensity value for the averaged spectrum
                    zspec,z = self.roach[i].integrate_spectrum(grid_list,channel_list,n_channel_list,baseline_list,n_baseline_list,baseline_order,type)
                    pdata_list.append(z)
                    pspectrum_list.append(zspec)

            t_list.append(pt_list)
            n_list.append(len(xgrid)*len(ygrid))
            bufpos_list.append(pbufpos_list)
            x_list.append(px_list)
            y_list.append(py_list)
            test_x_list.append(ptest_x_list)
            test_y_list.append(ptest_y_list)

            
            data_list.append(pdata_list)
            spectrum_list.append(pspectrum_list)

        self.map_pixel_list = np.array(mp_list)
        self.map_t = np.array(t_list)
        self.map_bufpos = np.array(bufpos_list)
        self.map_x = np.array(x_list)
        self.map_y = np.array(y_list)
        self.map_test_x = np.array(test_x_list)
        self.map_test_y = np.array(test_y_list)
        self.map_n = np.array(n_list)
        self.map_data = np.array(data_list)
        self.map_spectra = np.array(spectrum_list)

    def create_map_grid_data(self,channel_list,n_channel_list,baseline_list,n_baseline_list,baseline_order,pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],type=0):

        t_list = []
        data_list = []
        bufpos_list = []
        x_list = []
        y_list = []
        n_list = []
        mp_list = []
        spectrum_list = []

        for ipix in pixel_list:
            i = self.find_pixel_index(ipix)
            mp_list.append(ipix)
            pt_list = []
            pn_list = []
            pbufpos_list = []
            px_list = []
            py_list = []
            pdata_list = []
            pspectrum_list = []
            # find the grid point data
            b = np.bincount(self.roach[i].bufpos)
            e = np.where(b>0)[0]
            h = np.where(e>=100)[0]
            hlen = len(h)
            for ih in range(hlen):
                grid_list = np.where(self.roach[i].bufpos[self.roach[i].ons] == e[h[ih]])
                pt_list.append(np.mean(self.roach[i].spec_time[self.roach[i].ons][grid_list]))
                pn_list.append(len(self.roach[i].xmap[self.roach[i].ons][grid_list]))
                pbufpos_list.append(np.mean(self.roach[i].bufpos[self.roach[i].ons][grid_list]))
                px_list.append(np.mean(self.roach[i].xmap[self.roach[i].ons][grid_list]))
                py_list.append(np.mean(self.roach[i].ymap[self.roach[i].ons][grid_list]))
                # z is the integrated intensity value for the averaged spectrum
                zspec,z = self.roach[i].integrate_spectrum(grid_list,channel_list,n_channel_list,baseline_list,n_baseline_list,baseline_order,type)
                pdata_list.append(z)
                pspectrum_list.append(zspec)

            t_list.append(pt_list)
            n_list.append(hlen)
            bufpos_list.append(pbufpos_list)
            x_list.append(px_list)
            y_list.append(py_list)
            data_list.append(pdata_list)
            spectrum_list.append(pspectrum_list)

        self.map_pixel_list = np.array(mp_list)
        self.map_t = np.array(t_list)
        self.map_bufpos = np.array(bufpos_list)
        self.map_x = np.array(x_list)
        self.map_y = np.array(y_list)
        self.map_n = np.array(n_list)
        self.map_data = np.array(data_list)
        self.map_spectra = np.array(spectrum_list)
    

    def create_grid(self):
        """ creates the map grid from map parameters
        """
        xepsilon = self.ifproc.xlength*0.001 # to allow arange to find the last point
        xstep = self.xstep*self.hpbw
        self.xgrid = np.arange(-self.ifproc.xlength/2.,self.ifproc.xlength/2.+xepsilon,xstep)
        self.nx = len(xgrid)
        yepsilon = self.ifproc.ylength*0.001 # to allow arange to find the last point
        ystep = self.ystep*self.hpbw
        self.ygrid = np.arange(-self.ifproc.ylength/2.,self.ifproc.ylength/2.+yepsilon,ystep)
        self.ny = len(ygrid)
        self.ngrid = self.nx*self.ny

    def find_grid_location(self,ix,iy):
        return(self.xgrid[ix],self.ygrid[iy])

    def find_grid_point_id(self,ix,iy):
        """ returns bufpos number of grid point ix,iy
        """
        index = iy*self.nx+ix
        return(index+100)



class SpecBankCal(SpecBank):
    """ this BASE class deals with a complete "bank" of spectra """
    def __init__(self,roach_files,ifproc_data,pixel_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],bdrop=50,edrop=100):
        """ inputs: 
        roach files  is list of files for the roach boards (nominally 4 of them)
        ifproc_data  is data from the corresponding ifproc file
        pixel_list   is a list of pixels to process
        """
        SpecBank.__init__(self,roach_files,ifproc_data,pixel_list=pixel_list)
    
        # go ahead and process the calibration files
        for ipix in range(self.npix):
            self.roach[ipix].compute_tsys_spectrum(bdrop,edrop)

    def test_cal(self,DATA_FILE):
        result = 0
        if DATA_FILE.line_rest_frequency != self.line_rest_frequency:
            print('Warning: Line Rest Frequency of Cal observation does not match map')
            result = result + 1
        if DATA_FILE.nchan != self.nchan:
            print('Warning: number of channels for Cal observation does not match map')
            result = result + 1
        if DATA_FILE.bandwidth != self.bandwidth:
            print('Warning: spectrometer bandwidth of Cal observation does not match map')
            result = result + 1
        return(result)
