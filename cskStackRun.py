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

steplist = ['init','create','crop','master','focus_split','geo2rdr_coarseResamp','dense_offset','postprocess', 'check_dense_offset', 'focus_all','geocode','plot_offsetfield','create_stack']
nprocess = {}
nprocess['init'] = 1
nprocess['create'] = 1
nprocess['crop'] = 1
nprocess['master'] = 1
nprocess['focus_split'] = 8
nprocess['geo2rdr_coarseResamp'] = 1
nprocess['dense_offset'] = 6
nprocess['postprocess'] = {'geometry':1, 'maskandfilter': 12}
nprocess['check_dense_offset'] = 1
nprocess['focus_all'] = 8
nprocess['geocode'] = 8
nprocess['plot_offsetfield'] = 10
nprocess['create_stack'] = 1

def createParser():
    parser = argparse.ArgumentParser( description='control the running of the stacks step list: [init(0), create(1), crop(2), master(3), focus_split(4), geo2rdr_coarseResamp(5), dense_offset(6), postprocess(7), check_dense_offset(8), focus_all(9), geocode(10), plot_offsetfield(11), create_stack(12)]')
    
    parser.add_argument('-p','--proj', dest='project',type=str,help='project E(Evans) or R(Rutford)',required=True)
    
    parser.add_argument('-fs', dest='fs',type=int,default=0,help='the first track')
    parser.add_argument('-ls', dest='ls',type=int,default=0,help='the last track')
    parser.add_argument('-is',dest="istack", type=int, default=0, help="the number of stack in this track")
    parser.add_argument('-t',dest="track",type=str, default="", help="the track to process")
    
    parser.add_argument('-first', dest='first',type=str,help='the first step',default=None)
    parser.add_argument('-last', dest='last',type=str,help='the last step',default=None)
    parser.add_argument('-do', dest='do',type=str,help='do this step',default=None)

    parser.add_argument('--runid', dest='runid',type=int,help='runid for dense offset',default=None)
    parser.add_argument('--version', dest='version',type=str,help='version for dense offset processing',default=None)

    parser.add_argument('--onlyGeo2rdr', dest='onlyGeo2rdr',type=bool,default=False, help='only perform geo2rdr and skip resampling')
   
    parser.add_argument('-e', dest='exe',type=bool, help='execute the command or not', default=False)

    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)

def init(name, workdir):
    print("Prepare readme for: ", name)
    raw_folder = os.path.join(old_workdir, name, "raw")

    # Find master
    if proj == "CSK-Rutford":
        option = 1
    elif proj == "CSK-Evans":
        option = 2
    else:
        raise Exception()

    # choose the first one
    if option == 0:
        files = sorted(glob.glob(os.path.join(workdir, name, 'raw','2*')))
        master = files[0].split('/')[-1]

    # choose the same one as the old setup
    elif option == 1:
        f = open(os.path.join(old_workdir, name, "readme"))
        line = f.readlines()[0]
        f.close()
        elements = line.split()
        for i, element in enumerate(elements):
            if element == "-m" or element == "--master":
                master = elements[i+1]
                break
    # choose from the provided list
    elif option == 2:
        master = None
        f = open("master_list.txt")
        lines = f.readlines()
        f.close()
        for line in lines:
            try:
                line_name, line_master = line.split()
                if line_name == name:
                    master = line_master
            except:
                pass
    else:
        raise Exception()
    
    print("master: ", master)

    if proj == "CSK-Rutford":
        dem = "/net/kraken/nobak/mzzhong/CSK-Rutford-v2/DEM/Rutford_bedmachine_surface.dem"
    elif proj == "CSK-Evans":
        dem = "/marmot-nobak/mzzhong/CSK-Evans-v3/DEM/Evans_bedmachine_surface.dem"
    else:
        raise Exception()

    if proj == "CSK-Rutford":
        stack_cmd = 'stackStripMap.py -s {raw_folder} --bbox "-78.3 -73.5 -85 -65" -d {dem} -w . -m {master} -W slc --useGPU'.format(raw_folder=raw_folder, dem=dem, master=master)
    elif proj == "CSK-Evans":
        stack_cmd = 'stackStripMap.py -s ./raw --bbox "-78.3 -73.5 -85 -65" -d {dem} -w . -m {master} -W slc --useGPU'.format(dem=dem, master=master)
    
    print(stack_cmd)

    os.chdir(name)
    readme = open('readme','w')
    readme.write(stack_cmd + '\n')
    readme.close()

    os.chdir('..')

def create(name):
    os.chdir(name)
    cmd = 'cat '+ 'readme | bash'
    os.system(cmd)
    os.chdir('..')

def crop(name):
    raw_folder = os.path.join(old_workdir, name, "raw")
    cmd = "ln -s {} raw_crop".format(raw_folder)
    print(cmd)
    
    os.chdir(name)
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

    if inps.project.lower() in ["r" or "rutford"]:
        proj = "CSK-Rutford"
        full_track_list = [8,10,23,25,40,52,55,67,69,82,97,99,114,126,128,129,141,143,156,158,171,172,173,186,188,201,203,215,218,230,231,232]

    elif inps.project.lower() in ['e' or "evans"]:
        proj = "CSK-Evans"
        full_track_list = range(22)
    else:
        raise Exception()

    ################################################################ 
    # runid should be associated with a set up dense offset parameters
    # version should be associated post-filtering

    #### About runid ####
    # 960 x 480
    # runid = 20190900
    #params: ww=960, wh=480, sw=20, sh=20, kw=240, kh=120
    # runid = 20190900

    ############################
    # 480 x 240

    # runid = 20190901
    #params: ww=480, wh=240, sw=20, sh=20, kw=240, kh=120
    #runid = 20190901

    # runid = 201909010
    #params: ww=480, wh=240, sw=16, sh=16, kw=240, kh=120
    # for csk-r-v2, new dem, new ampcor
    #runid = 201909010

    # runid = 201909011
    #params: ww=480, wh=240, sw=16, sh=16, kw=120, kh=60
    #params: smaller skip size
    #runid = 201909011

    # runid = 201909012
    #params: ww=480, wh=240, sw=16, sh=16, kw=240, kh=120
    # Use the staging PyCuAmpcor (isce2_exp2) with removal of n2 and n4 in variance estimation
    #runid = 201909012

    # runid = 201909013
    # params: ww=480, wh=240, sw=16,sh =16, kw=240, kh=120
    # Use the staging PyCuAmpcor (isce3_exp3) with bug-fix in covariance calculation
    #runid = 201909013

    # runid = 201909014
    #params: ww=480, wh=240, sw=16, sh=16, kw=120, kh=60
    #params: smaller skip size
    # Run with the latest GPU
    #runid = 201909014

    #############################
    # runid = 20190904
    #params: ww=256, wh=128, sw=20, sh=20, kw=128, kh=64
    #runid = 20190904

    #############################
    # runid = 20190908
    # params: ww=128, wh=128, sw=20, sh=20, kw=64, kh=64
    #runid = 20190908

    #############################
    # runid = 20190921
    # params: ww=128, wh=64, sw=20, sh=20, kw=64, kh=32
    #runid = 20190921

    #############################
    # runid = 20190925
    # params: ww=64, wh=64, sw=20, sh=20, kw=32, kh=32
    #runid = 20190925


    ############# About version ################
    # v12 is standard version
    #version="v12"

    # v13 for 201909011 is to use doubled size median-filtering
    # from (5,5) to (9,9)
    #version="v13"

    # v14 for 201909011 and 201909014 is to increase the trim size from 22 to 25
    #version='v14'


    if proj == "CSK-Rutford":
        # Set up the project.
        stack = "stripmap"
    
        old_workdir = "/net/kraken/nobak/mzzhong/CSK-Rutford"
    
        #workdir = "/net/kraken/nobak/mzzhong/CSK-Rutford"
        workdir = "/net/kraken/nobak/mzzhong/CSK-Rutford-v2"

        runid = inps.runid

        version = inps.version
 
    elif proj=="CSK-Evans":

        # Set up the project.
        stack = "stripmap"
    
        # Old setup in 2018
        old_workdir = "/net/kraken/nobak/mzzhong/CSK-Evans"
        #runid = 20180712
    
        workdir = "/marmot-nobak/mzzhong/CSK-Evans-v3"

        runid = inps.runid

        version = inps.version
    else:
        raise Exception()
        
    ######################################################################

    # Find the tracks to process
    # full track
    if inps.track=="full":
        tracks = full_track_list
    # provided but not full
    elif inps.track!="":
        tracks = [int(track) for track in inps.track.split(',')]
    # then look at inps.fs and inps.ls (for mode like CSK-Evans)
    else:
        tracks = range(inps.fs,inps.ls+1)

    # Find the steps
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
            offset = dense_offset(stack=stack, workdir = workdir, nproc = nprocess['dense_offset'], runid=runid, version=version, exe=inps.exe)
        
        if stepname == 'postprocess':
            offset = dense_offset(stack=stack, workdir = workdir, nproc = nprocess['postprocess'], runid=runid, version=version, exe=inps.exe)

        if stepname == 'check_dense_offset':
            offset = dense_offset(stack=stack, workdir = workdir, nproc = nprocess['check_dense_offset'], runid=runid, version=version, exe=inps.exe)

        if stepname == 'geocode':
            offset = dense_offset(stack=stack, workdir=workdir, nproc = nprocess['geocode'], runid=runid, version=version, exe = inps.exe)

        if stepname == 'plot_offsetfield':
            offset = dense_offset(stack=stack, workdir=workdir, nproc = nprocess['plot_offsetfield'], runid=runid, version=version, exe = inps.exe)

        if stepname == 'create_stack':
            offset = dense_offset(stack=stack, workdir=workdir, nproc = nprocess['create_stack'], runid=runid, version=version, exe = inps.exe)

        for i in tracks:
            for j in range(1):
                name = "track_" + str(i).zfill(3) + '_' + str(inps.istack)
                
                print("track name: ", name)

                # in case the folder doesn't exist
                if os.path.exists(name) == False:
                    continue

                print(stepname)
                if stepname == 'init':
                    #print("not allowed to run")
                    #print(stop)
                    init(name, workdir)

                elif stepname == 'create':
                    # serial
                    create(name)

                    # multi-proc
                    #func = globals()[steplist[step]]
                    #p = Process(target=func, args=(name,))
                    #jobs.append(p)
                    #p.start()
                    #count = count + 1
                    #
                    ## nprocess controller
                    #if count == nprocess[steplist[step]]:
                    #    for ip in range(nprocess[steplist[step]]):
                    #        jobs[ip].join() #wait here
                    #    count = 0
                    #    jobs = []

                elif stepname == "crop":
                    crop(name)

                elif stepname == 'dense_offset':
                    offset.initiate(trackname = name)
                    offset.run_offset_track()

                elif stepname == 'postprocess':
                    offset.initiate(trackname = name)
                    offset.postprocess()

                elif stepname == 'check_dense_offset':
                    offset.initiate(trackname = name)
                    offset.run_check_offset_track()

                elif steplist[step] == 'geocode':
                    
                    offset.initiate(trackname = name)
                    offset.geocode()

                elif steplist[step] == 'plot_offsetfield':
                    
                    offset.initiate(trackname = name)
                    last = (i == tracks[-1])
                    offset.plot_offsetfield()

                elif steplist[step] == 'create_stack':
                    
                    offset.initiate(trackname = name)
                    offset.createOffsetFieldStack()

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
