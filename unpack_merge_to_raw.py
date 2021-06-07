#!/usr/bin/env python3
import datetime
from datetime import date

import os
import sys
import glob

from collections import defaultdict

# ISCE
import isce
from isceobj.Sensor import createSensor
import shelve
import pickle
from isceobj.Util import Poly1D
from isceobj.Planet.AstronomicalHandbook import Const

def unpack_merge_to_raw(h5file_list, raw_data_dir, name, sate, remove=False):
    for h5file in h5file_list:
        print(h5file)
    print("raw_data_dir: ", raw_data_dir)
    print("name: ", name)

    raw_data_name = os.path.join(raw_data_dir, name+'.raw')

    redo = False
    if (not os.path.exists(raw_data_name) or redo == True) and remove==False:

        if not os.path.exists(raw_data_dir):
            os.mkdir(raw_data_dir)

        obj = createSensor('COSMO_SKYMED')
        obj.hdf5FileList = h5file_list
        obj.output = raw_data_name
        obj.extractImage()
        obj.frame.getImage().renderHdr()
        obj.extractDoppler()

        pickName = os.path.join(raw_data_dir, 'raw')

        with shelve.open(pickName) as db:
            db['frame'] = obj.frame
    
        with open(pickName+'.pkl','wb') as f:
            pickle.dump(obj.frame, f)
        print("raw object: ", dir(obj))

        # Save the key parameters such as azimuthPixelSize and rangePixelSize
        f = open(raw_data_dir + "/params.txt","w")
        f.write("azimuth pixel size: " + str(obj.azimuthPixelSize) + "\n")
        f.write("range pixel size: " + str(obj.rangePixelSize) + "\n")
        f.write("satellite name: " + sate + "\n")
        f.write("sampling rate: " + str(obj.samplingRate) + "\n")
        f.write("PRF: " + str(obj.PRF) + "\n")
        f.write("ground velocity: " + str(obj.VG) + "\n")
        f.write("satellite speed: " + str(obj.speed) + "\n")
        f.write("satellite height: " + str(obj.height) + "\n")
        f.close()
    else:
        print(raw_data_name, 'exists')

    if remove:
        if os.path.exists(raw_data_name):
            cmd = "rm " + raw_data_name
            print(cmd)

    return 0
