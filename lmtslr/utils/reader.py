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
        # reduce all spectra - uncalibrated - using IFProc Tsys
        calobsnum=I.calobsnum
        ifproc_cal_file = lookup_ifproc_file(calobsnum,path=path+'ifproc/')
        ICal = IFProcCal(ifproc_cal_file)
        ICal.compute_tsys()
        for ipix in range(S.npix):
            S.roach[ipix].reduce_ps_spectrum(type=2,normal_ps=True,calibrate=False,tsys_no_cal=ICal.tsys[ipix])

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
        calobsnum = I.calobsnum
        ifproc_cal_file = lookup_ifproc_file(calobsnum,path=path+'ifproc/')
        ICal = IFProcCal(ifproc_cal_file)
        ICal.compute_tsys()
        S.roach[0].reduce_ps_spectrum(type=2,normal_ps=False,calibrate=False, tsys_no_cal=ICal.tsys[0])
        S.roach[1].reduce_ps_spectrum(type=2,normal_ps=True,calibrate=False,tsys_no_cal=ICal.tsys[1])

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
        calobsnum = I.calobsnum
        ifproc_cal_file = lookup_ifproc_file(calobsnum,path=path+'ifproc/')
        ICal = IFProcCal(ifproc_cal_file)
        ICal.compute_tsys()
        for ipix in range(S.npix):
            S.roach[ipix].reduce_spectra(type=1,calibrate=False,tsys_no_cal=ICal.tsys[ipix])

    return I,S

def read_obsnum_otf_multiprocess(I,Ical,obsnum,list_of_pixels,bank,use_calibration,tsys=150.,path='/data_lmt/'):
    """ Reads the spectral line data from WARES spectrometer for a particular obsnum.
    
    This is a useful utility function for use in line processing. It combines a number
    of tasks that are commonly done together to read a set of spectral line data.
    
    This version created for multiprocessing.  IFProc data is passed as an argument

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

    # create the spec_bank object.  This reads all the roaches in the list "files"
    S = SpecBankData(files,I,pixel_list=list_of_pixels,bank=bank)

    # check whether to use calibration and open necessary file
    if use_calibration == True:
        S.cal_flag = False
        calobsnum = S.calobsnum
        cal_files,ncalfiles = lookup_roach_files(calobsnum,roach_list,path=path+'spectrometer/')
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
            S.roach[ipix].reduce_spectra(type=1,calibrate=False,tsys_no_cal=Ical.tsys[ipix])

    return S

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

