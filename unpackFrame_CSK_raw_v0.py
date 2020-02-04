#!/usr/bin/env python3

# Minyan 2018-04-20, to support list of products

import isce
from isceobj.Sensor import createSensor
import shelve
import pickle
import argparse
import glob
from isceobj.Util import Poly1D
from isceobj.Planet.AstronomicalHandbook import Const
import os

from CSK_Utils import CSK_Utils

def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description='Unpack CSK data and store metadata in pickle file.')
    
    parser.add_argument('-i','--input', dest='h5dir', type=str,
            required=True, help='Input CSK directory, where all CSK*.h5 files are located (if multiple, merge)')

    parser.add_argument('-o', '--output', dest='slcdir', type=str,
            required=True, help='Output SLC directory')

    parser.add_argument('-name', '--name', dest='name', type=str,
            required=True, help='Name of the merged raw file, (name.raw)')

    parser.add_argument('-iseg', '--iseg', dest='iseg', type=int,
            required=True, help='number of segment of this track')

    #parser.add_argument('-usePickle', '--usePickle', dest='usePickle', action='store_true', default=False, help='To use Pickle instead of Shelve')

    #parser.add_argument('--ad', '--ad', dest='name', type=int,
    #        required=True, help='Asecending or Descending track')

    return parser.parse_args()


def unpack(h5dir, slcdir, name, iseg):
    '''
    Unpack HDF5 to binary SLC file.
    '''

    csk = CSK_Utils()

    if not os.path.isdir(slcdir):
        os.mkdir(slcdir)

    # h5dir is regex expression
    # dirlist is expanded directory
    dirlist = glob.glob(h5dir)
    dirlist.sort()

    # Cut a segment of the list to unpack
    dirlist = dirlist[iseg*csk.maxframenum:(iseg+1)*csk.maxframenum]

    # Find the h5 product list
    h5list=[]
    for dirname in dirlist:
        h5name = glob.glob(os.path.join(dirname,'*h5'))
        h5list = h5list + h5name

    if len(h5list)==0:
        return

    productdir = os.path.join(slcdir,name)
    if not os.path.isdir(productdir):
        os.mkdir(productdir)

    # Unpack here 
    obj = createSensor('COSMO_SKYMED')

    obj.hdf5FileList = h5list
    obj.output = os.path.join(productdir, name+'.raw')

    obj.extractImage()

    obj.frame.getImage().renderHdr()

    obj.extractDoppler()

    # Save the shelve and pickle file    
    pickName = os.path.join(productdir, 'raw')
    
    with shelve.open(pickName) as db:
        db['frame'] = obj.frame
    
    with open(pickName+'.pkl','wb') as f:
        pickle.dump(obj.frame, f)

    
    #cmd = 'fixImageXml.py -i ' + name + '.raw -f'
    #print(cmd)
    #os.system(cmd)


if __name__ == '__main__':
    '''
    Main driver.
    '''

    inps = cmdLineParse()

    # Remove the last "/"
    if inps.slcdir.endswith('/'):
        inps.slcdir = inps.slcdir[:-1]

    unpack(inps.h5dir, inps.slcdir, inps.name, inps.iseg)
