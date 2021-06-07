#!/usr/bin/env python3

import os
import sys
import re

import numpy as np
import time
from CSK_Utils import CSK_Utils
import datetime
from datetime import date

import glob

from multiprocessing import Process
import argparse

from dense_offset import dense_offset

# Time and dates.
fmt = '%Y%m%d'
today  = datetime.datetime.now().strftime(fmt)

server = 'https://qc.sentinel1.eo.esa.int/'
queryfmt = '%Y-%m-%d'
datefmt = '%Y%m%dT%H%M%S'

S1Astart = '20140901'
S1Astart_dt = datetime.datetime.strptime(S1Astart, '%Y%m%d')


S1Bstart = '20160501'
S1Bstart_dt = datetime.datetime.strptime(S1Bstart, '%Y%m%d')

## Set up the project and runid
proj_name = "S1-Evans"
#proj_name = 'LA_Basin'

if proj_name == "S1-Evans":
    stack = 'tops'
    old_workdir = '/net/kraken/nobak/mzzhong/S1-Evans'
    workdir = '/net/kraken/nobak/mzzhong/S1-Evans-v2'

    # runid = 20180703
    # params: ww=256 wh=128 sw=10 sh=10 kw=128 kh=64
    #runid = 20180703

    # runid = 20200101
    # params: ww=256 wh=128 sw=10 sh=10 kw=128 kh=64
    # run with new GPU ampcor
    #runid = 20200101

    # runid = 20200102
    # params: ww=480 wh=128 sw=10 sh=10 kw=256 kh=64
    # run with new GPU ampcor
    #runid = 20200102

    # runid = 20200103
    # params: ww=480 wh=120 sw=10 sh=10 kw=240 kh=60
    # run with new gpu ampcor that fixes uncertainty estimation
    #runid = 20200103

    # runid = 20200104
    # params: ww=480 wh=120 sw=10 sh=10 kw=240 kh=60
    # run with new gpu ampcor thats allows large chip size 
    # and some small fixes in uncertainty estimation
    # e.g. change covariance invalid value from 99 to 0 
    # (2021.05.11)
    #runid = 20200104

    # runid = 20200105
    # params: ww=960 wh=240 sw=10 sh=10 kw=240 kh=60
    # run with new gpu ampcor (exp_6)
    runid = 20200105

elif proj_name == "RC":
    stack = 'tops_RC'
    workdir = '/net/kraken/nobak/mzzhong/RC_denseOffsets'
    runid = 20190725

elif proj_name == "LA_Basin":
    stack = 'tops_LA_Basin'
    workdir = '/net/kraken/nobak/mzzhong/LA_Basin'
    runid = 20201001
    # ww = 128 wh = 32 kw = 64 kh = 32 sw=10 sh = 10, os = 64
else:
    raise Exception("Unknown project name")

# Version
version='v12'

# Steps and processors.
steplist = ['init','create','unpack_slc_topo_master', 'average_baseline', 'geo2rdr_resample', 'extract_stack_valid_region', 'merge_master_slave_slc', 'dense_offsets', 'postprocess','geocode','plot_offsetfield','create_stack','show_doc']
nprocess = {}

nprocess['init'] = 1
nprocess['create'] = 1
nprocess['unpack_slc_topo_master'] = 8
nprocess['average_baseline'] = 12
nprocess['geo2rdr_resample'] = 3
nprocess['extract_stack_valid_region'] = 1
nprocess['merge_master_slave_slc'] = 1
nprocess['dense_offsets'] = 6
nprocess['postprocess'] = {'geometry':1, 'maskandfilter': 12}
nprocess['geocode'] = 8
nprocess['plot_offsetfield'] = 12
nprocess['create_stack'] = 1
nprocess['show_doc'] = 1

def createParser():

    parser = argparse.ArgumentParser( description='control the running of the stacks step list: [init(0), create(1), unpack_slc_topo_master(2), average_baseline(3), geo2rdr_resample(4), extract_stack_valid_region(5), merge_master_slave_slc(6), dense_offsets(7), postprocess(8), geocode(9),  plot_offsetfield(10), create_stack(11), show_doc(12)]')
    
    parser.add_argument('-t', '--tracks', dest='tracks',type=str,help='tracks to process(comma separated)',default=None)
    parser.add_argument('-p', '--pairs', dest='pairs',type=str,help='pairs to process(comma separated)',default=None)

    parser.add_argument('-first', dest='first',type=str,help='the first step',default=None)
    parser.add_argument('-last', dest='last',type=str,help='the last step',default=None)
    parser.add_argument('-do', dest='do',type=str,help='do this step',default=None)

    parser.add_argument('--s_date', dest='s_date', type=str, default=S1Astart, help='Start date')
    parser.add_argument('--e_date', dest='e_date', type=str, default=today, help='Stop date')

    parser.add_argument('-e', dest='exe',type=bool, help='execute the command or not', default=False)

    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)


def init(name):
    
    os.chdir(name)

    if name == 'track_37':
        
        stack_cmd = 'stackSentinel.py -s ../data_37/ -d ../DEM/Evans.dem -a ../AuxDir/ -o ../Orbits/ -b "-78.3 -73.5 -85 -65" -W offset -C geometry -p hh -m 20170904 --useGPU --start_date 2017-09-04'
        readme = open('readme','w')
        readme.write(stack_cmd)
        readme.close()

    else:
        print('unknown track')
        print(stop)

    os.chdir('..')

def create(name):

    os.chdir(name)
    cmd = 'cat '+ 'readme | bash'
    os.system(cmd)
    os.chdir('..')

def run_cmd(cmd,start=None,end=None, exe=False):
    
    #print(cmd)
    #time.sleep(np.random.randint(low=2, high=6))

    cmd = cmd[:-1]

    if start and end:
        cmd = cmd + ' -s Function-' + str(start) + ' -e Function-' + str(end)

    if exe:
        os.system(cmd)
    else:
        print(cmd)

    return 0

def check_exist(step,line):

    exist = False
    start = None
    end = None

    if step == 5:
        return (exist, start, end)

    config_file = line.split(' ')[-1]
    config_file = config_file[0:-1]

    f = open(config_file)
    if steplist[step] == 'unpack_slc_topo_master':
        print("check the unpacked slc")

        params = f.readlines()

        master = False
        for ip in range(len(params)):
            if params[ip][0:6] == 'outdir':
                outdir = params[ip].split()[-1]
                xml1 = outdir + '/IW1.xml'
                xml2 = outdir + '/IW2.xml'
                xml3 = outdir + '/IW3.xml'
            if params[ip][0:11] == 'geom_master':
                geom_master = params[ip].split()[-1]
                IW1 = geom_master + '/IW1'
                IW2 = geom_master + '/IW2'
                IW3 = geom_master + '/IW3'

                master = True

        # Need to count the number of files.
        if (os.path.exists(xml1) and os.path.exists(xml2) and os.path.exists(xml3)):
            if master == True and os.path.exists(IW1) and os.path.exists(IW2) and os.path.exists(IW3):
                exist = True
            if master == True and not (os.path.exists(IW1) and os.path.exists(IW2) and os.path.exists(IW3)):
                exist = False
                start = 2
                end = 2
            if master == False:
                exist = True

        else:
            exist = False
            

    if steplist[step] == 'average_baseline':

        #print('check the baseline')

        params = f.readlines()
        for ip in range(len(params)):
            if params[ip][0:4] == 'base':
                baseline_file = params[ip].split()[-1]
                
                exist = False
                if os.path.exists(baseline_file): 
                    # need to go into the txt file to check
                    num_lines = sum(1 for line in open(baseline_file))
                    if num_lines == 9:
                        exist = True


    if steplist[step] == 'geo2rdr_resample':

        exist1 = False
        exist2 = False

        # first count the total number of burst
        params = f.readlines()
        tot = {}
        # geordr step
        tot_azi = {}
        tot_rng = {}
        tot_slc = {}

        for ip in range(len(params)):
            if params[ip][0:6] == 'master':
                master = params[ip].split()[-1]

            if params[ip][0:5] == 'coreg':
                coregSLC = params[ip].split()[-1]
        
        for i in range(1,4):
            tot[i] = len(glob.glob(master+'/IW'+str(i)+'/*xml'))
            tot_azi[i] = len(glob.glob(coregSLC+'/IW'+str(i)+'/azimuth*off.xml'))
            tot_rng[i] = len(glob.glob(coregSLC+'/IW'+str(i)+'/range*off.xml'))
            tot_slc[i] = len(glob.glob(coregSLC+'/IW'+str(i)+'/*slc.xml'))

        xml1 = coregSLC + '/IW1.xml'
        xml2 = coregSLC + '/IW2.xml'
        xml3 = coregSLC + '/IW3.xml'

        print(tot,tot_azi,tot_rng,tot_slc)

        geo2rdr = False
        resample = False

        if (tot_azi == tot) and (tot_rng == tot):
            geo2rdr = True

        if (tot_azi == tot) and (tot_rng == tot) and (tot_slc == tot) and os.path.exists(xml1) and os.path.exists(xml2) and os.path.exists(xml3):
            resample = True

        if geo2rdr and resample:
            exist = True
        elif geo2rdr:
            exist = False
            start = 2
            end = 2
        else:
            exist = False

        print(geo2rdr, resample, exist)


    if steplist[step] == 'merge_master_slave_slc':
        
        print("check the merged slc")

        params = f.readlines()
        for ip in range(len(params)):
            if params[ip][0:7] == 'outfile':
                outdir = params[ip].split()[-1]
                xml1 = outdir + '.full.aux.xml'

                if os.path.exists(xml1):
                    exist = True
                # Note that topo product should be checked here ideally
        
    f.close()

    return (exist,start,end)

def main(iargs=None):

    ## Prepare the parameters.

    inps = cmdLineParse(iargs)

    # the tracks to process
    if inps.tracks is not None:
        tracks = inps.tracks.split(',')

    if inps.pairs is not None:
        pairs = inps.pairs.split(',')
    else:
        pairs = None

    # the steps
    if ( inps.first == None and inps.last == None ) and inps.do == None:
        print('please provide the steps')
        return 0

    elif inps.do:
        first_step = inps.do
        last_step = inps.do
    else:
        first_step = inps.first
        last_step = inps.last

    # use string step and number step
    if not (first_step.isdigit() and last_step.isdigit()):
        first = steplist.index(first_step)
        last = steplist.index(last_step)
    else:
        first = int(first_step)
        last = int(last_step)
    
    if first > last:
        print('the first step should not be later than the last step')
        return 0

    # Convert to datetime object.
    s_date = datetime.datetime.strptime(inps.s_date, fmt)
    e_date = datetime.datetime.strptime(inps.e_date, fmt)


    # Used for simple multithread-processing.
    count = 0
    jobs = []

    for step in range(first,last+1):
        stepname = steplist[step]

        print("step: ",stepname)

        if stepname == 'dense_offsets':

            offset = dense_offset(stack=stack, workdir=workdir, nproc = nprocess['dense_offsets'], runid = runid, version=version, exe = inps.exe)

        if stepname == 'postprocess':

            offset = dense_offset(stack=stack, workdir=workdir, nproc = nprocess['postprocess'], runid = runid, version=version, exe = inps.exe)

        if stepname == 'geocode':

            offset = dense_offset(stack=stack, workdir=workdir, nproc = nprocess['geocode'], runid=runid, version=version, exe = inps.exe)

        if stepname == 'plot_offsetfield':
            
            offset = dense_offset(stack=stack, workdir=workdir, nproc = nprocess['plot_offsetfield'],runid = runid, version=version, exe = inps.exe)

        if stepname == 'create_stack':

            offset = dense_offset(stack=stack, workdir=workdir, nproc = nprocess['create_stack'], runid=runid, version=version, exe = inps.exe)

        if stepname == 'show_doc':
            
            offset = dense_offset(stack=stack, workdir=workdir, nproc = nprocess['show_doc'],runid = runid, version=version, exe = False)


        # loop through the targes
        for i in tracks:
            name='track_' + str(i)

            if pairs:
                name = name+'_pairs' +'/' + pairs[0]
                print("track name: ", name)
 
            else:
                print("track name: ", name)

                # in case the folder doesn't exist
                if os.path.exists(name) == False:
                    continue

                print("current step: ", steplist[step])
                if stepname == 'init' or stepname == 'create':
                    print('step init and step create are temporarily forbidden')
                    print(stop)

                    func = globals()[steplist[step]]
                    p = Process(target=func, args=(name,))
                    jobs.append(p)
                    p.start()
                    count = count + 1
                    
                    # nprocess controller
                    if count == nprocess[steplist[step]]:
                        for ip in range(nprocess[steplist[step]]):
                            jobs[ip].join() #wait here
                        count = 0
                        jobs = []

                elif steplist[step] == 'dense_offsets':
                   
                    offset.initiate(trackname = name)
                    offset.run_offset_track()

                elif steplist[step] == 'postprocess':

                    offset.initiate(trackname = name)
                    offset.postprocess(s_date=s_date,e_date=e_date)


                elif steplist[step] == 'geocode':
                    
                    offset.initiate(trackname = name)
                    offset.geocode(s_date=s_date,e_date=e_date)

                elif steplist[step] == 'plot_offsetfield':
                    
                    offset.initiate(trackname = name)
                    offset.plot_offsetfield()

                elif steplist[step] == "create_stack":

                    offset.initiate(trackname = name)
                    offset.createOffsetFieldStack()

                elif steplist[step] == 'show_doc':
                    
                    offset.initiate(trackname = name)
                    offset.show_doc()
 
                else:
                    cmd_file = os.path.join(name,'run_files','run_'+str(step-1)+'_'+steplist[step])
                    f = open(cmd_file)
                    lines = f.readlines()

                    for ii in range(len(lines)):
                        print('===========================')

                        cmd = lines[ii]

                        # checking if results exist
                        (exist, start, end) = check_exist(step,cmd)

                        if exist:
                            print('exists, skip ',cmd)
                            continue

                        # If running with geo2rdr, check the data
                        if steplist[step] == 'geo2rdr_resample':
                            datestr = cmd.split('_')[-1][:-1]
                            print(datestr)
                            this_date = datetime.datetime.strptime(datestr, fmt)

                            if this_date < s_date or this_date>e_date:
                                print('out of the time window')
                                print('skip ',cmd)
                                continue
 
                        print('run ', cmd)
 
                        p = Process(target=run_cmd, args=(cmd,start,end,inps.exe))
                        jobs.append(p)
                        p.start()
                        count = count + 1

                        # nprocess controller
                        if count == nprocess[steplist[step]]:
                            while True:
                                for ip in range(nprocess[steplist[step]]):
                                    if not jobs[ip].is_alive():
                                        count = count - 1
                                        del jobs[ip]
                                        break
                                if count < nprocess[steplist[step]]:
                                    break

                    f.close()

if __name__=='__main__':
    main()
