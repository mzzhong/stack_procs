#!/usr/bin/env python3

#Cunren Liang, 10-MAY-2018

import os
import sys
import datetime
import argparse


def runCmd(cmd):
    
    print("{}".format(cmd))
    status = os.system(cmd)
    if status != 0:
        raise Exception('error when running:\n{}\n'.format(cmd))


def runRecordCmd(cmd):

    with open('cmd.sh', 'a') as f:
        header  = '###################################################################\n'
        header += '#  commands of processing started at: {}\n'.format(datetime.datetime.now())
        header += '###################################################################\n'
        f.write(header)
        for x in cmd:
            f.write(x+'\n')
            runCmd(x)
        f.write('\n\n\n')


def cmdLineParse():
    '''
    Command line parser.
    '''

    parser = argparse.ArgumentParser(description='convert unwrapped interferogram in ISCE to GMT format')
    parser.add_argument('-unw', dest='unw', type=str, required=True,
            help = 'unwrapped and geocoded interferogram ISCE.')
    parser.add_argument('-wvl', dest='wvl', type=float, required=True,
            help = 'radar wavelength in meters.')
    parser.add_argument('-gv', dest='gv', type=int, default=0,
            help = 'gmt version. 0: < 5 (default). 1: >= 5.')  

    if len(sys.argv) <= 1:
        #print('')
        parser.print_help()
        sys.exit(1)
    else:
        return parser.parse_args()


if __name__ == '__main__':

    inps = cmdLineParse()
 
    unw = inps.unw
    wvl  = inps.wvl
    gv  = inps.gv

    if gv == 0:
        gmt_prefix = ''
    else:
        gmt_prefix = 'gmt '

    cmd = []
    #convert to deformation, sign not changed
    cmd.append('''imageMath.py --eval="a_0;a_1/2.0/PI*{wvl}/2.0" --a={unw} -t float -s BIL -o unw.geo'''.format(
            wvl = wvl,
            unw = unw
            ))
    #convert to gmt format
    cmd.append("isce2gis.py envi -i unw.geo")
    cmd.append("gdal_translate -of GMT -b 2 unw.geo unw0.grd")
    cmd.append("grdmath unw0.grd 0 NAN = unw.grd".format(gmt_prefix))

    #tidy up
    cmd.append('mkdir unw_isce unw_gmt')
    cmd.append("mv unw.geo unw.geo.xml unw.geo.vrt unw.geo.hdr unw_isce")
    cmd.append("mv unw.grd unw_gmt")
    cmd.append("rm unw0.grd* isce.log")

    #run commands
    runRecordCmd(cmd)


