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

    parser = argparse.ArgumentParser(description='convert LOS in ISCE to ENU unit vector')
    parser.add_argument('-los', dest='los', type=str, required=True,
            help = 'geocoded LOS file in ISCE.')
    parser.add_argument('-la', dest='la', type=int, default=0,
            help = 'convert to LOS or radar flying direction ENU unit vector. 0: LOS (default). 1: radar flying direction') 
    parser.add_argument('-rl', dest='rl', type=int, default=0,
            help = 'right or left looking. 0: right (default). 1: left. no need to specify if convert to LOS ENU unit vector')  
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
 
    los = inps.los
    rl  = inps.rl
    la  = inps.la
    gv  = inps.gv

    if gv == 0:
        gmt_prefix = ''
    else:
        gmt_prefix = 'gmt '

    cmd = []
    #convert LOS angles to LOS ENU unit vector
    if la == 0:
        #after converstion, invalid value in first and second bands: 0
        #                   invalid value in third band:             1
        cmd.append("imageMath.py --eval='sin(rad(a_0))*cos(rad(a_1+90));sin(rad(a_0))*sin(rad(a_1+90));cos(rad(a_0))' --a={los} -t float -s BIL -o enu.geo".format(
            los = los
            ))
    #convert azimuth angle to radar flying direction ENU unit vector
    else:
        if rl == 0:
            add90 = -90.0
        else:
            add90 = 90.0
        #after converstion, invalid value in first and second bands: 0
        #                   invalid value in third band:             1
        cmd.append("imageMath.py --eval='(a_0!=0)*cos(rad(a_1+90+({add90})));(a_0!=0)*sin(rad(a_1+90+({add90})));0+(a_0==0)' --a={los} -t float -s BIL -o enu.geo".format(
            add90 = add90,
            los = los
            ))

    #convert ENU vectors to GMT format
    cmd.append("isce2gis.py envi -i enu.geo")

    cmd.append("gdal_translate -of GMT -b 1 enu.geo east0.grd")
    cmd.append("{}grdmath east0.grd 0 NAN = east.grd".format(gmt_prefix))

    cmd.append("gdal_translate -of GMT -b 2 enu.geo north0.grd")
    cmd.append("{}grdmath north0.grd 0 NAN = north.grd".format(gmt_prefix))

    #note that above conversion sets "no data" to 1.0 in up component
    cmd.append("gdal_translate -of GMT -b 3 enu.geo up0.grd")
    cmd.append("{}grdmath up0.grd 1.0 NAN = up.grd".format(gmt_prefix))

    #tidy up
    cmd.append('mkdir enu_isce enu_gmt')
    cmd.append("mv enu.geo enu.geo.xml enu.geo.vrt enu.geo.hdr enu_isce")
    cmd.append("mv east.grd north.grd up.grd enu_gmt")
    cmd.append("rm east0.grd* north0.grd* up0.grd* isce.log")

    #run commands
    runRecordCmd(cmd)


