#!/usr/bin/env python3

import os
import sys

folder=sys.argv[1]

def main():

    for iw in ['IW1','IW2','IW3']:
        swath_folder = '/'.join([folder, iw]) 
        for filename in os.listdir(swath_folder):
            if filename.endswith('.vrt'):
                burst='.'.join(filename.split('.')[:-1])
                print(burst)
                inputfile='/'.join([swath_folder, burst+'.vrt'])
                outputfile='/'.join([swath_folder, burst])
                if not os.path.exists(outputfile):
                    #print(inputfile, outputfile)
                    cmd='gdal_translate -of ENVI -co INTERLEAVE=BIL ' + inputfile + ' ' + outputfile
                    print(cmd)
                    os.system(cmd)
                else:
                    print('skip: ', outputfile)



if __name__=='__main__':
    main()
