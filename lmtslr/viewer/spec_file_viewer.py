import numpy as np
import matplotlib.pyplot as pl
import netCDF4

class SpecFileViewer():
    def __init__(self, netcdf_filename):
        nc = netCDF4.Dataset(netcdf_filename, 'r', format='NETCDF4')
        self.obsnum = nc.variables['Header.Obs.ObsNum'][0]
        self.source_name = netCDF4.chartostring(nc.variables['Header.Obs.SourceName'][:])

        self.x_position =  nc.variables['Header.Obs.XPosition'][0]
        self.y_position =  nc.variables['Header.Obs.YPosition'][0]

        self.nchan = nc.variables['Header.Line.NChannels'][0]
        self.chan = nc.variables['Header.Line.ChannelNumber'][:]
        self.cdelt = nc.variables['Header.SpectrumAxis.CDELT'][0]
        self.crpix = nc.variables['Header.SpectrumAxis.CRPIX'][0]
        self.crval = nc.variables['Header.SpectrumAxis.CRVAL'][0]
        self.ctype = netCDF4.chartostring(nc.variables['Header.SpectrumAxis.CTYPE'][:])
        self.caxis = nc.variables['Header.SpectrumAxis.CAXIS'][:]

        self.pixel = nc.variables['Data.Pixel'][:]
        self.xpos = nc.variables['Data.XPos'][:]
        self.ypos = nc.variables['Data.YPos'][:]
        self.rms = nc.variables['Data.RMS'][:]
        self.data = nc.variables['Data.Spectra'][:]
        
        nc.close()

    def sequoia_waterfall_plot(self,pixel_list,rms_cut,plot_range=[-1,1],figsize=8):
        fig1, ax1 = pl.subplots(4, 4, sharex='col', sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(figsize,figsize))
        fig1.text(0.02, 0.5, self.ctype, va='center', rotation='vertical')
        fig1.text(0.5, 0.1, 'Sample', ha='center')
        for the_pixel in pixel_list:
            pindex = np.where(self.pixel==the_pixel)[0]
            rindex = np.where(self.rms[pindex]<rms_cut)[0]
            ax1[np.mod(the_pixel,4), the_pixel//4].imshow(self.data[pindex[rindex]].transpose(),origin='lower',extent=[0.,float(len(rindex)),self.caxis[0],self.caxis[-1]],clim=plot_range,aspect='auto')
            ax1[np.mod(the_pixel,4), the_pixel//4].text(0.05*len(rindex),self.caxis[0] + 0.85*(self.caxis[-1]-self.caxis[0]), '%d'%(the_pixel))

    def pixel_waterfall_plot(self,the_pixel,rms_cut,plot_range=[-1,1]):
        pl.figure()
        pindex = np.where(self.pixel==the_pixel)[0]
        rindex = np.where(self.rms[pindex]<rms_cut)[0]
        pl.imshow(self.data[pindex[rindex]].transpose(),origin='lower',extent=[0.,float(len(rindex)),self.caxis[0],self.caxis[-1]],clim=plot_range,aspect='auto')
        pl.title('PIXEL: %d'%(the_pixel))
        pl.ylabel(self.ctype)
        pl.xlabel('Sample')
        pl.colorbar()

    def sequoia_rms_plot(self,pixel_list,rms_cut,plot_range=[0.,10.],figsize=8):
        fig2, ax2 = pl.subplots(4, 4, sharex='col', sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(figsize,figsize))
        fig2.text(0.02, 0.5, 'RMS', va='center', rotation='vertical')
        fig2.text(0.5, -0.1, 'Sample', ha='center')

        for the_pixel in pixel_list:
            pindex = np.where(self.pixel==the_pixel)[0]
            rindex = np.where(self.rms[pindex]<rms_cut)[0]
            ax2[np.mod(the_pixel,4), the_pixel//4].plot(self.rms[pindex[rindex]],'k.')
            #ax2[np.mod(the_pixel,4), the_pixel//4].text(0.05*len(rindex),plot_range[0] + 0.9*(plot_range[-1]-plot_range[0]), '%d'%(the_pixel))
            #ax2[np.mod(the_pixel,4), the_pixel//4].ylim(plot_range)

    def pixel_rms_plot(self,the_pixel,rms_cut,plot_range=[0.,10.]):
        pl.figure()
        pindex = np.where(self.pixel==the_pixel)[0]
        rindex = np.where(self.rms[pindex]<rms_cut)[0]
        pl.plot(self.rms[pindex[rindex]],'k.')
        pl.ylim(plot_range)
        pl.ylabel('RMS')
        pl.xlabel('Sample')
        pl.title('PIXEL: %d'%(the_pixel))
        
    def xy_position_plot(self):
        pl.figure()
        pl.plot(self.xpos,self.ypos,'k.')
        pl.xlabel('X')
        pl.ylabel('Y')

    def sequoia_rms_histogram(self,pixel_list,rms_cut,figsize=8):
        fig3, ax3 = pl.subplots(4, 4, sharex='col', sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(figsize,figsize))
        fig3.text(0.5, -0.1, 'RMS', ha='center')
        for the_pixel in pixel_list:
            pindex = np.where(self.pixel==the_pixel)[0]
            rindex = np.where(self.rms[pindex]<rms_cut)[0]
            ax3[np.mod(the_pixel,4),the_pixel//4].hist(self.rms[pindex[rindex]],bins=np.arange(0.,3.02,.02))

    def pixel_rms_histogram(self,the_pixel,rms_cut):
        pl.figure()
        pindex = np.where(self.pixel==the_pixel)[0]
        rindex = np.where(self.rms[pindex]<rms_cut)[0]
        pl.hist(self.rms[pindex[rindex]],bins=np.arange(0.,3.02,.02))
        pl.xlabel('RMS')
        pl.ylabel('N')
        pl.title('PIXEL: %d'%(the_pixel))
        
    def sequoia_mean_spectra_plot(self, pixel_list, rms_cut, figsize=8):
        fig4, ax4 = pl.subplots(4, 4, sharex='col', sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(figsize,figsize))
        fig4.text(0.5, -0.1, self.ctype, ha='center')
        for the_pixel in pixel_list:
            pindex = np.where(self.pixel==the_pixel)[0]
            rindex = np.where(self.rms[pindex]<rms_cut)[0]
            ax4[np.mod(the_pixel,4),the_pixel//4].plot(self.caxis,np.mean(self.data[pindex[rindex]],axis=0))
 

    def pixel_mean_spectrum_plot(self, the_pixel, rms_cut):
        pl.figure()
        pindex = np.where(self.pixel==the_pixel)[0]
        rindex = np.where(self.rms[pindex]<rms_cut)[0]
        pl.plot(self.caxis, np.mean(self.data[pindex[rindex]], axis=0))
        pl.xlabel(self.ctype)
        pl.ylabel('TA*')
        pl.title('PIXEL: %d'%(the_pixel))
                
