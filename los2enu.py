#!/usr/bin/env python3

# Minyan Zhong , June 2018

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
    parser.add_argument('-rl', dest='rl', type=int, default=0,
            help = 'right or left looking. 0: right (default). 1: left. no need to specify if convert to LOS ENU unit vector')  

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

    filedir = os.path.dirname(los)
    losfile = os.path.basename(los)

    enu_los_file = os.path.join(filedir,'enu_' + losfile)
    enu_azi_file = os.path.join(filedir,'enu_' + losfile.replace('los','azi'))

    print(enu_los_file)
    print(enu_azi_file)

    cmd = []
    #convert LOS angles to LOS ENU unit vector
    
    #after converstion, invalid value in first and second bands: 0
    #                   invalid value in third band:             1
    
    cmd.append("imageMath.py --eval='sin(rad(a_0))*cos(rad(a_1+90));sin(rad(a_0))*sin(rad(a_1+90));cos(rad(a_0))' --a={los} -t float -s BIL -o {output}".format(
            los = los, output = enu_los_file
            ))
    
    #convert azimuth angle to radar flying direction ENU unit vector
    if rl == 0:
        add90 = -90.0
    else:
        add90 = 90.0
    
    #after converstion, invalid value in first and second bands: 0
    #                   invalid value in third band:             1
    cmd.append("imageMath.py --eval='(a_0!=0)*cos(rad(a_1+90+({add90})));(a_0!=0)*sin(rad(a_1+90+({add90})));0+(a_0==0)' --a={los} -t float -s BIL -o {output}".format(
            add90 = add90,
            los = los,
            output = enu_azi_file
            ))
    
    #run commands
    runRecordCmd(cmd)


