import os
import fnmatch

roach_pixels_all = [[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]]

def lookup_roach_files(obsnum,
                       roach_list=['roach0/', 'roach1/', 'roach2/', 'roach3/'],
                       path='/data_lmt/spectrometer/'):
    """ stand alone function to look up roach files which match a particular obsnum
        obsnum     is number to be matched
        roach_list is a list of the names of the directories for the search
        path       is the path to the roach directories
        returns: list of file names and number of files found
    """
    nroach = len(roach_list)
    filenames = []
    result = 0
    for i in range(nroach):
        for file in os.listdir(path+roach_list[i]):
            prefix = roach_list[i][:-1]
            if fnmatch.fnmatch(file,'%s_%d_*.nc'%(prefix, obsnum)):
                print('found %s'%(file))
                if  not 'allantest' in file:
                    print('append %s'%(file))
                    filenames.append(path+roach_list[i]+file)
                    result = result + 1
    if filenames == []:
        print('lookup_roach_files: no files for obsnum', obsnum)
    return(filenames,result)


def find_roach_from_pixel(pixel_id):
    """ find roach from pixel 
    """
    for i,l in enumerate(roach_pixels_all):
        if pixel_id in l:
            return [i]

    return []

def create_roach_list(pixel_list):
    """ build a list of roach boards to be read given a pixel list
    """
    rid = [0,0,0,0]
    for pixel_id in pixel_list:
        r = find_roach_from_pixel(pixel_id)
        if r != []:
            rid[r[0]] = 1
    roach_list = []
    for i in range(4):
        if rid[i] == 1:
            roach_list.append('roach%d/'%(i))
    return roach_list
