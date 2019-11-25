#!/usr/bin/env python3

import os
import sys
import re

import numpy as np
import time
import datetime
from datetime import date

import glob

from multiprocessing import Process
import argparse

from dense_offset import dense_offset

proj = "CSK-Rutford"
#proj = "CSK-Evans"

# runid should be associated with a set up dense offset parameters
# version should be associated post-filtering
if proj == "CSK-Rutford":
    # Set up the project.
    stack = "stripmap"
    workdir = "/net/kraken/nobak/mzzhong/CSK-Rutford"
    
    # runid = 20190908
    # params: ww=128, wh=128, sw=20, sh=20, kw=64, kh=64
    #runid = 20190908

    # runid = 20190921
    # params: ww=128, wh=64, sw=20, sh=20, kw=64, kh=32
    runid = 20190921

    # runid = 20190925
    # params: ww=64, wh=64, sw=20, sh=20, kw=32, kh=32
    #runid = 20190925

elif proj=="CSK-Evans":
    # Set up the project.
    stack = "stripmap"
    workdir = "/net/kraken/nobak/mzzhong/CSK-Evans"
    runid = 20180712

steplist = ['init','create','crop','master','focus_split','geo2rdr_coarseResamp','dense_offset','postprocess','focus_all','geocode','plot_geocoded']
nprocess = {}
nprocess['init'] = 1
nprocess['create'] = 1
nprocess['crop'] = 8
nprocess['master'] = 1
nprocess['focus_split'] = 8
nprocess['geo2rdr_coarseResamp'] = 1
nprocess['dense_offset'] = 6
nprocess['postprocess'] = {'geometry':1, 'maskandfilter': 12}
nprocess['focus_all'] = 8
nprocess['geocode'] = 8
nprocess['plot_geocoded'] = 1

def createParser():
    parser = argparse.ArgumentParser( description='control the running of the stacks step list: [init(0), create(1), crop(2), master(3), focus_split(4), geo2rdr_coarseResamp(5), dense_offset(6), postprocess(7), focus_all(8), geocode(9), plot_geocoded(10)]')
    parser.add_argument('-fs', dest='fs',type=int,default=0,help='the first track')
    parser.add_argument('-ls', dest='ls',type=int,default=0,help='the last track')
    parser.add_argument('-is',dest="istack", type=int, default=0, help="the number of stack in this track")
    parser.add_argument('-t',dest="track",type=str, default="", help="the track to process")
    
    parser.add_argument('-first', dest='first',type=str,help='the first step',default=None)
    parser.add_argument('-last', dest='last',type=str,help='the last step',default=None)
    parser.add_argument('-do', dest='do',type=str,help='do this step',default=None)

    parser.add_argument('--onlyGeo2rdr', dest='onlyGeo2rdr',type=bool,default=False, help='only perform geo2rdr and skip resampling')
   
    parser.add_argument('-e', dest='exe',type=bool, help='execute the command or not', default=False)

    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)


def init(name, workdir):
    
    os.chdir(name)

    files = glob.glob(os.path.join('.','raw','2*'))
    master = files[0].split('/')[-1]
    print(master)

    stack_cmd = 'stackStripMap.py -s ./raw --bbox "-78.3 -73.5 -85 -65" -d ' + workdir + '/Evans_DEM_const/Evans_constant.dem -w . -m ' + master + ' -W slc --useGPU'
    readme = open('readme','w')
    readme.write(stack_cmd)
    readme.close()

    os.chdir('..')

def create(name):

    os.chdir(name)

    cmd = 'cat '+ 'readme | bash'
    os.system(cmd)

    os.chdir('..')


def run_cmd(cmd, start, end, exe=False):

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

    config_file = line.split(' ')[-1]
    config_file = config_file[0:-1]


    f = open(config_file)
    if steplist[step] == 'crop':
        "check the geom_master directory output: for example /net/kraken/nobak/mzzhong/CSK_ProcessDir2/track_200/raw_crop/20171213/20171213.raw"
        print('checking the cropped files')
        params = f.readlines()
        for ip in range(len(params)):
            if params[ip][0:6] == 'output':
                crop_dir = params[ip].split()[-1]

                date = crop_dir.split('/')[-1]

                rawxmlfile = crop_dir + '/' + date+'.raw.xml'
                rawfile = crop_dir + '/' + 'raw'
                
                if os.path.exists(rawxmlfile) and os.path.exists(rawfile):
                    exist = True
                    #print(stop)

    if steplist[step] == 'focus_split':
        
        print('check if focused')
        
        params = f.readlines()
        
        for ip in range(len(params)):
            if params[ip][0:5] == 'input':

                focus = params[ip].split()[-1]

                date = focus.split('/')[-1]

        slcfile = focus + '/' + date + '.raw.slc.xml'

        if os.path.exists(slcfile):
            exist = True

    if steplist[step] == 'master':
        
        print('checking the master focusing and geometry')
        
        params = f.readlines()
        
        for ip in range(len(params)):
            if params[ip][0:5] == 'input':

                focus = params[ip].split()[-1]

                date = focus.split('/')[-1]

            if params[ip][0:6] == 'output':
                geom_dir = params[ip].split()[-1]
                lat = geom_dir + '/lat.rdr.xml'
                lon = geom_dir + '/lon.rdr.xml'
                hgt = geom_dir + '/hgt.rdr.xml'
                los = geom_dir + '/los.rdr.xml'

        slcfile = focus + '/' + date + '.raw.slc.xml'

        if not os.path.exists(slcfile):
            exist = False

        elif os.path.exists(lat) and os.path.exists(lon) and os.path.exists(hgt) and os.path.exists(los):
            exist = True

        else:
            exist = False
            start = 2
            end = 2
    
    if steplist[step] == 'geo2rdr_coarseResamp':
        print('checking the coregistered slcs')

        params = f.readlines()
        for ip in range(len(params)):
            if params[ip][0:6] == 'outdir':

                offset = params[ip].split()[-1]
                azxml = offset + '/azimuth.off.xml'
                rngxml = offset + '/range.off.xml'

            if params[ip][0:5] == 'coreg':
                
                coreg_dir = params[ip].split()[-1]
                date = coreg_dir.split('/')[-1]
                slc = coreg_dir + '/' + date + '.slc.xml'

        if not (os.path.exists(azxml) and os.path.exists(rngxml)):
            exist = False

        elif os.path.exists(slc):
            exist = True
        else:
            exist = False
            start = 2
            end = 2

        # Skip the resampling process
        #if inps.onlyGeo2rdr == True:
        #    start =1
        #    end = 1
                
    f.close()
    
    return (exist, start, end)

def main(iargs=None):

    inps = cmdLineParse(iargs)

    # the tracks to process, first look at inps.track (for mode like CSK-Rutford)
    if inps.track!="":
        tracks = [int(track) for track in inps.track.split(',')]
    # then look at inps.fs and inps.ls (for mode like CSK-Evans)
    else:
        tracks = range(inps.fs,inps.ls+1)

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

    # Used for simple multihreading processing.
    jobs = []
    count = 0

    for step in range(first,last+1):

        stepname = steplist[step]
        
        # loop through the tracks
        if stepname == 'dense_offset':
            offset = dense_offset(stack=stack, workdir = workdir, nproc = nprocess['dense_offset'], runid=runid, exe=inps.exe)
        
        if stepname == 'postprocess':
            offset = dense_offset(stack=stack, workdir = workdir, nproc = nprocess['postprocess'], runid=runid, exe=inps.exe)

        if stepname == 'geocode':
            offset = dense_offset(stack=stack, workdir=workdir, nproc = nprocess['geocode'], runid=runid, exe = inps.exe)

        if stepname == 'plot_geocoded':
            offset = dense_offset(stack=stack, workdir=workdir, nproc = nprocess['plot_geocoded'], runid=runid, exe = inps.exe)

        for i in tracks:

            for j in range(1):

                if proj == "CSK-Rutford":
                    name = "track_" + str(i).zfill(3) + '_' + str(inps.istack)
                
                elif proj == "CSK-Evans":
                    name='track_' + str(i).zfill(2) + str(j)
                
                else:
                    raise Exception()
                
                print("track name: ", name)

                # in case the folder doesn't exist
                if os.path.exists(name) == False:
                    continue

                #print(step)
                #print(nprocess[steplist[step]])
                # controled by run_files

                print(stepname)
                if stepname == 'init':
                    print('forbidden')
                    print(stop)
                    init(name,workdir)

                elif stepname == 'create':

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

                elif stepname == 'dense_offset':
                    offset.initiate(trackname = name)
                    offset.run_offset_track()

                elif stepname == 'postprocess':
                    offset.initiate(trackname = name)
                    offset.postprocess()

                elif steplist[step] == 'geocode':
                    
                    offset.initiate(trackname = name)
                    offset.geocode()

                elif steplist[step] == 'plot_geocoded':
                    
                    offset.initiate(trackname = name)
                    last = (i == tracks[-1])
                    #offset.plot_geocoded(label='all', last = last)
                    offset.plot_geocoded(label='separate')

                elif stepname == "crop":
                    os.system("/net/kraken/nobak/mzzhong/CSK-Rutford/createFolder.py")

                elif stepname == 'master' or stepname == 'focus_split' or stepname=='geo2rdr_coarseResamp':

                    cmd_file = os.path.join(name,'run_files','run_'+str(step-1)+'_'+steplist[step])
                    f = open(cmd_file)
                    lines = f.readlines()
                    for ii in range(len(lines)):

                        exist, start, end = check_exist(step,lines[ii])
                        if exist:
                            print('skip: ',steplist[step])
                            continue

                        p = Process(target=run_cmd, args=(lines[ii],start,end,inps.exe))
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

                elif stepname == 'focus_all':
                    # focus all master and slave

                    print("focus all is closed")
                    print(stop)

                    # master step
                    step_master = steplist.index('master')
                    cmd_file = os.path.join(name,'run_files','run_'+str(step_master-1)+'_'+steplist[step_master])

                    f = open(cmd_file)
                    lines = f.readlines()
                    for ii in range(len(lines)):

                        exist, start, end = check_exist(step_master,lines[ii])
                        
                        # If only the second step, skip it, else 
                        if exist == False and start == 2 and end == 2:
                            exist = True
                        elif exist == False:
                            start = 1
                            end = 1

                        if exist:
                            print('skip: ',steplist[step_master])
                            continue

                        p = Process(target=run_cmd, args=(lines[ii],start,end,inps.exe))
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

                    #################################### 

                    # slave
                    step_slave = steplist.index('focus_split') 
                    cmd_file = os.path.join(name,'run_files','run_'+str(step_slave-1)+'_'+steplist[step_slave])
                    

                    f = open(cmd_file)
                    lines = f.readlines()
                    for ii in range(len(lines)):

                        exist, start, end = check_exist(step_slave,lines[ii])
                        if exist:
                            print('skip: ',steplist[step_slave])
                            continue

                        p = Process(target=run_cmd, args=(lines[ii],start,end,inps.exe))
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
