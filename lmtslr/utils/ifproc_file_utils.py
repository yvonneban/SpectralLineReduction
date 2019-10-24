import os
import fnmatch

def lookup_ifproc_file(obsnum,path='/data_lmt/ifproc/'):
    filename = ''
    for file in os.listdir(path):
        if fnmatch.fnmatch(file,'*_%06d_*.nc'%(obsnum)):
            print('found %s'%(file))
            filename = path+file
    if filename == '':
        print('lookup_ifproc_file: no file for obsnum ', obsnum)
        if 'lmttpm' not in path:
            print('look in lmttpm')
            return lookup_ifproc_file(obsnum,path='/data_lmt/lmttpm/')
    return(filename)
