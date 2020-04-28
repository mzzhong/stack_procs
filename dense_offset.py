#!/usr/bin/env python3
# Author: Minyan Zhong
import os
import sys
import numpy as np
import glob
import pathlib

import datetime

from datetime import date

import multiprocessing
from multiprocessing import Process

import subprocess

import pickle

from scipy import ndimage 
from scipy import signal

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm

import isce
import isceobj

import gdal

from medianfilter2d import medianfilter2d as mdf

import time

thrs_default=0.05
# v10: use _offset_filter_v2
# The mis-coregistration is assumed to be a constant, stationary is criterion is set to 0.05

#version = 'v11'
#thrs_default=0.1
# v11: use _offset_filter_v2
# The mis-coregistration is assumed to be a constant, stationary is criterion is set to 0.1

#version = 'v12'
#thrs_default=0.05
# v12: use _offset_filter_v2
# The mis-coregistration is assumed to be a plane (first order), stationary criterion is set to 0.1
# The mis-coregistration is derived from 20190901

nanvalue=-999

disp_temp_folder = "/net/kraken/nobak/mzzhong/CSK-Rutford/disp_temp_files"

#if len(os.listdir(disp_temp_folder))>0:
#    os.system("rm " + disp_temp_folder + "/*")

postproc_verbose = False
check_offset_verbose = True
 
def run_denseoffset_bash(cmd,exe=False,chmod=False):

    if exe:
        if type(cmd) == list:
            for this_cmd in cmd:
                subprocess.call(this_cmd, shell=True)

        else:
            #os.system("chmod u+x " + cmd)
            subprocess.call(cmd, shell=True)
    else:
        print(cmd)

class dense_offset_config():
    def __init__(self,ww=128, wh=128, sw=20, sh=20, kw=64, kh=64, runid=None, nwac=5, nwdc=5, oversample=32):

        self.outprefix = '_'.join(['cuampcor', 'ww'+str(ww), 'wh'+str(wh), 'os'+str(oversample)])
        
        print(self.outprefix)

        self.runid = runid
        self.ww = ww
        self.wh = wh
        self.sw = sw
        self.sh = sh
        self.kw = kw
        self.kh = kh

        self.oversample = oversample

        self.mm= max(max(self.sh,self.sw)*2, 50)
        self.gpuid = 0
        self.gross = 0
        self.deramp = 0

        self.nwac = nwac
        self.nwdc = nwdc

class dense_offset():

        def __init__(self, stack=None, workdir=None, nproc=None, runid=None, version=None, exe=False):

            self.stack = stack
            self.workdir = workdir
            self.nproc = nproc
            self.version = version
            self.exe = exe

            self.doc_count = 0
            self.jobs = []

            self.offsetFolder = 'cuDenseOffsets'

            self.filt_suffix = str(runid) + '_' + str(version)
            self.run_mode = "stack"
            self.rb_suffix = None 

            if stack == 'stripmap':
                self.slcs_folder = 'coregSLC/Coarse'
                self.master_folder = 'raw'
                self.slave_folder = self.slcs_folder

                self.master_suffix = '.raw.slc'
                self.slave_suffix = '.slc'

                self.maxday = 8

                # create the simplest doc object
                # 2019* CSK
                if runid==20190901:
                    self.doc = dense_offset_config(ww=480, wh=240, sw=20, sh=20, kw=240, kh=120, runid=runid)
                elif runid==20190904:
                    self.doc = dense_offset_config(ww=256, wh=128, sw=20, sh=20, kw=128, kh=64, runid=runid)
                elif runid==20190908:
                    self.doc = dense_offset_config(ww=128, wh=128, sw=20, sh=20, kw=64, kh=64, runid=runid)
                elif runid == 20190921:
                    self.doc = dense_offset_config(ww=128, wh=64, sw=20, sh=20, kw=64, kh=32, runid=runid)
                elif runid == 20190925:
                    self.doc = dense_offset_config(ww=64, wh=64, sw=20, sh=20, kw=32, kh=32, runid=runid)

                # Old setup for Evans
                elif runid == 20180712:
                    self.doc = dense_offset_config(ww=256, wh=256, sw=10, sh=10, kw=128, kh=128, runid=runid)

                else:

                    raise Exception("Undefined runid: ", runid)

                self.geometry = 'merged/geom_master'
                self.latname = 'lat.rdr'
                self.lonname = 'lon.rdr'
                self.losname = 'los.rdr'
                self.hgtname = 'hgt.rdr'

                self.geometry_suffix = '.rdr' 
 

            elif stack == 'tops':
                self.slcs_folder = 'merged/SLC'
                self.master_folder = self.slcs_folder
                self.slave_folder = self.slcs_folder

                self.master_suffix = '.slc.full'
                self.slave_suffix = '.slc.full'

                self.maxday = 12

                # create the simplest doc object
                # 2020* S1
                if runid == 20180703:

                    self.doc = dense_offset_config(ww=256, wh=128, sw=10, sh=10, kw=128, kh=64, runid=runid)

                elif runid == 20200101:
                    self.doc = dense_offset_config(ww=256, wh=128, sw=10, sh=10, kw=128, kh=64, runid=runid)

                elif runid == 20200102:
                    self.doc = dense_offset_config(ww=480, wh=128, sw=10, sh=10, kw=240, kh=64, runid=runid)

                else:
                    raise Exception("Undefined runid: ", runid)

                self.geometry = 'merged/geom_master'
                self.latname = 'lat.rdr.full'
                self.lonname = 'lon.rdr.full'
                self.losname = 'los.rdr.full'
                self.hgtname = 'hgt.rdr.full'

                self.geometry_suffix = '.rdr.full'

                # used to generate pixel size
                self.burst_xml = 'master/IW1.xml'

            elif stack == 'tops_RC':
                self.slcs_folder = 'merged/SLC'
                self.master_folder = self.slcs_folder
                self.slave_folder = self.slcs_folder

                self.master_suffix = '.slc.full'
                self.slave_suffix = '.slc.full'

                self.maxday = 16

                # create the simplest doc object
                ww = 64
                wh = 16
                oversample = 128
                self.doc = dense_offset_config(ww=ww, wh=wh, sw=20, sh=20, kw=ww//2, kh=wh//2, oversample=oversample)

                self.geometry = 'merged/geom_master'
                self.latname = 'lat.rdr.full'
                self.lonname = 'lon.rdr.full'
                self.losname = 'los.rdr.full'
                self.hgtname = 'hgt.rdr.full'

                self.geometry_suffix = '.rdr.full'

                # used to generate pixel size
                self.burst_xml = 'master/IW1.xml'

        
        def initiate_rb_pair(self, iter_num, rb_suffix):

            self.run_mode = "pair_rb"
            self.rb_iter_num = iter_num
            self.rb_suffix = rb_suffix

        def get_offset_geometry(self):

            doc = self.doc

            # Point to the geometry files.
            latfile = os.path.join(self.runfolder, self.geometry, 'lat' + self.geometry_suffix)
            lonfile = os.path.join(self.runfolder, self.geometry, 'lon' + self.geometry_suffix)
            losfile = os.path.join(self.runfolder, self.geometry, 'los' + self.geometry_suffix)
            hgtfile = os.path.join(self.runfolder, self.geometry, 'hgt' + self.geometry_suffix)

            doc.latfile = latfile 
            doc.lonfile = lonfile  
            doc.losfile = losfile  
            doc.hgtfile = hgtfile  

            # Offsetfield size.
            numWinDown = (doc.azsize - doc.mm*2 - doc.sh*2 - doc.wh) // doc.kh
            numWinAcross = (doc.rngsize - doc.mm*2 - doc.sw*2 - doc.ww) // doc.kw

            #numWinDown = 704

            # Update doc.
            doc.numWinDown = numWinDown
            doc.numWinAcross = numWinAcross

            print("Window width, Window height, Skip width, Skip height")
            print(doc.ww, doc.wh, doc.kw, doc.kh)
            print("Number of winows in down direction, Number of window in across direction")
            print(doc.numWinDown,doc.numWinAcross)

            # Check if the files exist.
            filelist = ['lat','lon','hgt','los']
            exist = True
            for geofile in filelist:
                offset_geofile = os.path.join(self.runfolder, self.geometry, geofile + '_offset_' + str(doc.runid) + '.rdr.xml')
                print(offset_geofile)

                if not os.path.exists(offset_geofile):
                    exist = False

            # Allocate new geometry files            
            lon = np.zeros(shape=(numWinDown,numWinAcross),dtype=np.float64)
            lat = np.zeros(shape=(numWinDown,numWinAcross),dtype=np.float64) 
            hgt = np.zeros(shape=(numWinDown,numWinAcross),dtype=np.float64)

            inc = np.zeros(shape=(numWinDown,numWinAcross),dtype=np.float32)
            azi = np.zeros(shape=(numWinDown,numWinAcross),dtype=np.float32)

            # Changed to track down the miscoregistration issue 
            # The issuse should be resolved after setting the starting pixel at (doc.mm, doc.mm)

            # option 1 oldest
            #centerOffsetHgt = doc.sh + doc.kh//2-1
            #centerOffsetWidth = doc.sw + doc.kw//2-1

            # option 2 used for RC earthquake
            #centerOffsetHgt = doc.sh + doc.wh//2-1
            #centerOffsetWidth = doc.sw + doc.ww//2-1

            # option 3
            centerOffsetHgt = doc.wh//2-1
            centerOffsetWidth = doc.ww//2-1

            print("numWinDown: ", numWinDown)
            print("nunWinAcross: ", numWinAcross)
            print("skip in width: ", doc.kw)
            print("skip in height: ", doc.kh)

            # Subset the original geometry file.
            downs = []
            for iwin in range(numWinDown):
                print(iwin)

                # The starting pixel is at doc.mm
                down = doc.mm + doc.kh * iwin  + centerOffsetHgt

                downs.append(down)
                
                off = down * doc.rngsize
                off2 = down * doc.rngsize * 2 # BIL
    
                # The starting pixel is at doc.mm
                range_indices = doc.mm + np.arange(numWinAcross) * doc.kw + centerOffsetWidth

                #print('margin size: ', doc.mm)
                #print('centerOffsetWidth: ', centerOffsetWidth)
                #print('full range size: ', doc.rngsize)
                #print(range_indices)
                #print(stop)
    
                latline = np.memmap(filename=latfile,dtype='float64',offset=8*off,shape=(doc.rngsize))
                lonline = np.memmap(filename=lonfile,dtype='float64',offset=8*off,shape=(doc.rngsize))
                hgtline = np.memmap(filename=hgtfile,dtype='float64',offset=8*off,shape=(doc.rngsize))
                 
                losline = np.memmap(filename=losfile,dtype='float32',offset=4*off2,shape=(doc.rngsize*2))

                #lat[iwin,:] = latline[start:end:doc.kw]
                #lon[iwin,:] = lonline[start:end:doc.kw]
                #hgt[iwin,:] = hgtline[start:end:doc.kw]

                lat[iwin,:] = latline[range_indices]
                lon[iwin,:] = lonline[range_indices]
                hgt[iwin,:] = hgtline[range_indices]

                incline = losline[0:doc.rngsize]
                aziline = losline[doc.rngsize:doc.rngsize*2]

                inc[iwin,:] = incline[range_indices]
                azi[iwin,:] = aziline[range_indices]


            # Give the arrays to doc project.
            doc.lat = lat
            doc.lon = lon
            doc.inc = inc
            doc.azi = azi

            # Concantenante los and azi (bsq).
            print(inc)
            print(azi)

            # Write to disk in isce format.
            for geofile in filelist:

                offset_geofile = os.path.join(self.runfolder, self.geometry, geofile + '_offset_' + str(doc.runid) + '.rdr')
                print(offset_geofile)

                if geofile == 'lat':
                    doc.offsetLatFile = offset_geofile

                    driver = gdal.GetDriverByName( 'ENVI' )
                    dst_ds = driver.Create(offset_geofile, lat.shape[1], lat.shape[0], 1, 6 )
                    dst_ds.GetRasterBand(1).WriteArray(lat, 0 ,0 )
                    dst_ds = None

                    # write to disk 
                    #write_cmd = geofile + '.tofile("' + offset_geofile + '")'
                    #exec(write_cmd)


                elif geofile == 'lon':
                    doc.offsetLonFile = offset_geofile
                    driver = gdal.GetDriverByName( 'ENVI' )
                    dst_ds = driver.Create(offset_geofile, lon.shape[1], lon.shape[0], 1, 6 )
                    dst_ds.GetRasterBand(1).WriteArray(lon, 0 ,0 )
                    dst_ds = None

                elif geofile == 'hgt':
                    doc.offsetHgtFile = offset_geofile
                    driver = gdal.GetDriverByName( 'ENVI' )
                    dst_ds = driver.Create(offset_geofile, hgt.shape[1], hgt.shape[0], 1, 6 )
                    dst_ds.GetRasterBand(1).WriteArray(hgt, 0 ,0 )
                    dst_ds = None


                elif geofile == 'los':
                    doc.offsetLosFile = offset_geofile
                    driver = gdal.GetDriverByName( 'ENVI' )
                    dst_ds = driver.Create(offset_geofile, inc.shape[1], inc.shape[0], 2, 6 )
                    dst_ds.GetRasterBand(1).WriteArray(inc, 0 ,0 )
                    dst_ds.GetRasterBand(2).WriteArray(azi, 0 ,0 )
                    dst_ds = None


                # Create xml.
                if geofile == 'los':
                    # create xml
                    outImg = isceobj.createImage()
                    outImg.setDataType('FLOAT')
                    outImg.setFilename(offset_geofile)
                    outImg.setBands(2)
                    outImg.scheme = 'BSQ'
                    outImg.setWidth(numWinAcross)
                    outImg.setLength(numWinDown)
                    outImg.setAccessMode('read')
                    outImg.renderHdr()
                    #outImg.renderVRT()
                    pass


                else:
                    # create xml 
                    outImg = isceobj.createImage()
                    outImg.setDataType('FLOAT')
                    outImg.setFilename(offset_geofile)
                    outImg.setBands(1)
                    outImg.scheme = 'BSQ'
                    outImg.setWidth(numWinAcross)
                    outImg.setLength(numWinDown)
                    outImg.setAccessMode('read')
                    outImg.renderHdr()
                    #outImg.renderVRT()
                    pass

            return 0


        def get_slc_size(self,xmlfile):
            
            import xml.etree.ElementTree as ET

            doc = self.doc

            tree = ET.parse(xmlfile)

            root = tree.getroot()

            for child in root:

                #print(child.attrib)
                if 'name' in child.attrib.keys() and child.attrib['name'] == 'coordinate1':
                    for grandchild in child:
                        if 'name' in grandchild.attrib.keys() and grandchild.attrib['name'] == 'size':
                            for element in grandchild:
                                if element.tag == 'value':
                                    doc.rngsize = int(element.text)

                if 'name' in child.attrib.keys() and child.attrib['name'] == 'coordinate2':
                    for grandchild in child:
                        if 'name' in grandchild.attrib.keys() and grandchild.attrib['name'] == 'size':
                            for element in grandchild:
                                if element.tag == 'value':
                                    doc.azsize = int(element.text)

            return 0

        def get_pixel_size(self):

            from iscesys.Component.ProductManager import ProductManager as PM
           
            doc = self.doc
 
            if self.stack.startswith('tops'):
                pm = PM()
                pm.configure()
                obj = pm.loadProduct(os.path.join(self.runfolder,self.burst_xml))
                burst = obj.bursts[0]
                Vs = np.linalg.norm(burst.orbit.interpolateOrbit(burst.sensingMid, method='hermite').getVelocity())
    
                print ('Velocity: ', Vs)
                print ('az Time interval: ', burst.azimuthTimeInterval)
                
                # Need to do correction for footprint!!!
                altitude = 693
                radius = 6371 

                doc.azPixelSize=Vs*burst.azimuthTimeInterval * radius / (radius + altitude)
                print ('azPixelSize: ', doc.azPixelSize)
                
                doc.rngPixelSize=burst.rangePixelSize
                print ('rngPixelSize: ', doc.rngPixelSize)

            elif self.stack == 'stripmap':

                param_file = os.path.join(self.trackfolder, 'params.txt')

                f = open(param_file,'r')
                params = f.readlines()
                doc.azPixelSize = round(float(params[0].split()[-1]),4)
                doc.rngPixelSize = round(float(params[1].split()[-1]),4)
                f.close()
                print ('azPixelSize: ', doc.azPixelSize)
                print ('rngPixelSize: ', doc.rngPixelSize)

            return 0

        def reference(self):
            
            from grossOffsets import grossOffsets

            doc = self.doc
            
            objGrossOff = grossOffsets(figpath = os.path.join(self.runfolder, self.offsetFolder), runid=doc.runid)

            mode = 'exterior'
            
            if mode == 'interior':
                objGrossOff.setMode(mode=mode)
                objGrossOff.setXSize(doc.rngsize)
                objGrossOff.setYize(doc.azsize)
                objGrossOff.setMargin(doc.mm)
                objGrossOff.setWinSizeHgt(doc.wh)
                objGrossOff.setWinSizeWidth(doc.ww)
                objGrossOff.setSearchSizeHgt(doc.sh)
                objGrossOff.setSearchSizeWidth(doc.sw)
                objGrossOff.setSkipSizeHgt(doc.kh)
                objGrossOff.setSkipSizeWidth(doc.kw)
                
                objGrossOff.setLatFile(doc.latfile)
                objGrossOff.setLonFile(doc.lonfile)
                objGrossOff.setLosFile(doc.losfile)

            elif mode == 'exterior':
                objGrossOff.setMode(mode=mode)
                objGrossOff.setOffsetLat(doc.lat)
                objGrossOff.setOffsetLon(doc.lon)
                objGrossOff.setOffsetInc(doc.inc)
                objGrossOff.setOffsetAzi(doc.azi)

                objGrossOff.setNumWinDown(doc.numWinDown)
                objGrossOff.setNumWinAcross(doc.numWinAcross)
            else:
                raise Exception('Wrong gross offset mode')

            # Setting the pixel size
            objGrossOff.setPixelSize(azPixelSize = doc.azPixelSize, rngPixelSize = doc.rngPixelSize)
            
            # Set the days
            objGrossOff.setbTemp(1)

            # Run to get reference offsetfield
            doc.grossDown, doc.grossAcross = objGrossOff.runGrossOffsets()

            return 0

        def show_doc(self):

            doc = self.doc
            print(doc.azsize)
            print(doc.rngsize)

        def prepare_master(self, run_folder):

            f = open(os.path.join(self.runfolder,'readme'))
            params = f.readlines()
            f.close()
            elements = params[0].split(' ')

            masterdate = None
            for i, element in enumerate(elements):
                if element == "-m":
                    masterdate = elements[i+1]
                    break

            if masterdate:

                # For stripmap, create the symbolic link to master
                if self.stack == 'stripmap':
    
                    masterfolder = os.path.join(self.runfolder, self.slcs_folder, masterdate)
    
                    # Create symbolic link to the master image data
                    if not os.path.exists(masterfolder):
                        os.mkdir(masterfolder)
                        linked_prefix = os.path.join(self.runfolder, self.master_folder, masterdate, masterdate)
                        linking_prefix = os.path.join(self.runfolder, self.slcs_folder, masterdate, masterdate)
     
                        cmd1 = 'ln -s ' + linked_prefix + '.raw.slc ' + linking_prefix + '.slc'
                        cmd2 = 'ln -s ' + linked_prefix + '.raw.slc.xml ' + linking_prefix + '.slc.xml'
                        cmd3 = 'ln -s ' + linked_prefix + '.raw.slc.vrt ' + linking_prefix + '.slc.vrt'

                        os.system(cmd1)
                        os.system(cmd2)
                        os.system(cmd3)

            else:
                raise Exception("Undefined master")

            return masterdate

        def initiate(self, trackname):

            doc = self.doc
            self.trackname = trackname
            self.tracknumber = int(trackname.split('_')[1])
            self.trackfolder = self.workdir + '/' + self.trackname

            # runfolder
            self.runfolder = self.trackfolder

            # create the offset folder 
            if not os.path.exists(self.trackfolder + '/' + self.offsetFolder):
                os.mkdir(self.trackfolder + '/' + self.offsetFolder)

            # Find the master date
            self.masterdate = self.prepare_master(run_folder = self.trackfolder)

            # Collect dates
            slcs = glob.glob(self.trackfolder + '/' + self.slcs_folder + '/2*')
            #for slc in slcs:
            #    if slc.split('/')[-1]=='20171119_backup':
            #        print(slc)
            #        raise Exception('Wrong Name')
            intdates = [int(slc.split('/')[-1]) for slc in slcs if len(slc.split('/')[-1]) == 8 ]
            
            intdates.sort() 
            self.obdates =[ date(int(str(intdate)[0:4]), int(str(intdate)[4:6]), int(str(intdate)[6:8])) for intdate in intdates]
            #print(self.obdates)

            # Create the all possible offsetfields.
            self.offsetfields = []
            for date1 in self.obdates:
                for date2 in self.obdates:
                    if (date2-date1).days >=1 and (date2-date1).days<=self.maxday:
    
                        date1str = date1.strftime('%Y%m%d')
                        date2str = date2.strftime('%Y%m%d')

                        self.offsetfields.append([date1str,date2str,date1,date2])

            # Preparation for offset field (pixelsize, geometry, etc)
            self.preparation_offsetfield(mode="stack")

        def initiate_track_pair(self, trackname, pairname, mode):

            doc = self.doc

            ## Initiate the basics
            self.trackname = trackname
            self.trackfolder = os.path.join(self.workdir, self.trackname) 
            self.pairname = pairname
            self.pairfolder = os.path.join(self.workdir, self.trackname, 'pairs', self.pairname)

            # run folder
            if mode=="pair":
                self.runfolder = self.pairfolder
            elif mode=="track":
                self.runfolder = self.trackfolder
            else:
                raise Exception("Undefined mode")

            # create the offset folder 
            pathlib.Path(self.runfolder + '/' + self.offsetFolder).mkdir(exist_ok=True)

            # Find the master date
            self.masterdate = self.prepare_master(run_folder = self.runfolder)

            date1str, date2str = pairname.split('_')
            date1 = datetime.datetime.strptime(date1str,"%Y%m%d").date()
            date2 = datetime.datetime.strptime(date2str,"%Y%m%d").date()

            self.offsetfields=[]
            self.offsetfields.append([date1str,date2str,date1,date2])

            self.preparation_offsetfield(mode=mode)

        def preparation_offsetfield(self, mode=None):

            ## offsetfield object initiation.
            # Load the existed result, if it exist. 
            # Point to the offset pickle file.

            # Control Redoing the preparation or not. 
            redo_preparation_offsetfield = 0
            
            doc = self.doc

            self.doc_pkl = os.path.join(self.runfolder, self.offsetFolder, str(doc.runid) + '.pkl')

            if os.path.exists(self.doc_pkl) and redo_preparation_offsetfield==0:
                print("Load doc object:", self.doc_pkl)
                with open(self.doc_pkl,'rb') as f:
                    self.doc = pickle.load(f)
                    doc = self.doc

                # Modification when necessary
                #doc.outprefix="cuampcor"
                # Save it
                #with open(self.doc_pkl ,'wb') as f:
                #    pickle.dump(doc,f)

            else:
                # Choose the xml file of the master and find its size
                xmls = glob.glob(self.runfolder + '/' + self.slcs_folder + '/' + self.masterdate + '/' +'*.xml')
 
                # Get the SLC image size
                self.get_slc_size(xmls[0])
            
                # Get the offset field pixel size
                self.get_pixel_size()
                print("Azimuth pixel size, Range pixel size")
                print(doc.azPixelSize, doc.rngPixelSize)

                # Get offset field geometry
                self.get_offset_geometry()
             
                # Obtain reference velocity
                if self.stack in ["stripmap","tops"]:
                    self.reference()
                else:
                    raise Exception("undefined stack")

                ##### The last chance to change doc object ######

                # Set runtime gpu params
                doc.nwac = 170
                doc.nwdc = 1

                # Save it
                with open(self.doc_pkl ,'wb') as f:
                    pickle.dump(doc,f)

            #print(stop)

            return 0

        def run_offset_track(self):

            doc = self.doc
            
            # Esitmate dense offsets for one track.
            for i_off, offsetfield in enumerate(self.offsetfields):
    
                date1str = offsetfield[0]
                date2str = offsetfield[1]

                offset_folder = os.path.join(self.runfolder,self.offsetFolder,date1str+'_'+date2str)
                
                if not os.path.exists(offset_folder):
                    os.mkdir(offset_folder)
                    print(offset_folder)

                # Generate runfile for dense_offset estimation.
                master = os.path.join(self.runfolder,self.slcs_folder,date1str, date1str + self.slave_suffix)
                slave = os.path.join(self.runfolder,self.slcs_folder,date2str, date2str + self.slave_suffix )

                offset_outprefix = os.path.join(offset_folder, doc.outprefix)
                offset_outsuffix = '_run_' + str(doc.runid)

                # Create runfile for this offset field
                run_file, offsetfield_filename  = self._cuampcor_runfile(offset_folder, master, slave, offset_outprefix, offset_outsuffix)

                redo=0
                if os.path.exists(offsetfield_filename) and redo==0:
                    print("existing: ", offsetfield_filename, "skip")
                else:
                    print("work on: ", offsetfield_filename)
                    #continue
                    func = globals()['run_denseoffset_bash']
                    p = Process(target=func, args=('bash '+ run_file, self.exe)) #default exe is False
                    self.jobs.append(p)
                    p.start()
                    self.doc_count = self.doc_count + 1
                    #time.sleep(1)
    
                    # nprocess controller
    
                    if self.doc_count == self.nproc or i_off == len(self.offsetfields)-1:
                        #for ip in range(self.nproc):
                        for ip in range(len(self.jobs)):
                            self.jobs[ip].join() #wait here
    
                        self.doc_count = 0
                        self.jobs = []

        def run_offset_pair_iter(self, iter_num, exe=False):

            doc = self.doc
            date1str = self.offsetfields[0][0]
            date2str = self.offsetfields[0][1]
            offset_folder = os.path.join(self.runfolder,self.offsetFolder,date1str+'_'+date2str)
            if not os.path.exists(offset_folder):
                os.mkdir(offset_folder)
                print(offset_folder)

            # Find the two slcs
            master = os.path.join(self.runfolder,self.slcs_folder,date1str, date1str + self.slave_suffix)
            if iter_num==0:
                slave = os.path.join(self.runfolder,self.slcs_folder,date2str, date2str + self.slave_suffix)
            else:
                # Rubber Sheeting obtained new slcs
                slave = os.path.join(self.runfolder, self.slcs_folder, date2str, date2str + '_' + "rb_" + str(iter_num-1) + self.slave_suffix)

            # Prefix & suffix
            rb_suffix = self.rb_suffix
            offset_outprefix = os.path.join(offset_folder,doc.outprefix)
            offset_outsuffix = '_run_' + str(doc.runid)+ '_' + rb_suffix

            # Generate runfile for dense_offset estimation.
            # Create runfile for this offset field
            run_file, offsetfield_filename  = self._cuampcor_runfile(offset_folder, master, slave, offset_outprefix, offset_outsuffix)

            redo_offset_pair_iter = 0
            if os.path.exists(offsetfield_filename) and redo_offset_pair_iter==0:
                print("skip: ", offsetfield_filename)
            else:
                print("run for: ", offsetfield_filename)
                run_denseoffset_bash(cmd="bash "+run_file, exe=exe)

            return 0

        def _cuampcor_runfile(self, offset_folder, master, slave, prefix, suffix): 

            doc = self.doc

            # set parameters
            # Get GPU id
            gpuid = self.doc_count % self.nproc + (8-self.nproc)
            print(gpuid)
            
            gpuAmpcor_driver = "/net/kamb/ssd-tmp1/mzzhong/insarRoutines/stack_procs/cuDenseOffsets_dev.py" 

            run_file = offset_folder + '/gpurun' + '_' + str(doc.runid) + '.sh'
            f = open(run_file,'w')
            f.write('#!/bin/bash\n')
            f.write('master=' + master+ '\n')
            f.write('slave=' + slave+ '\n')

            #f.write('fixImageXml.py -f -i $master'+ '\n')
            #f.write('fixImageXml.py -f -i $slave'+ '\n')

            f.write('ww=' + str(doc.ww)+ '\n')
            f.write('wh=' + str(doc.wh)+ '\n')
            f.write('sw=' + str(doc.sw)+ '\n')
            f.write('sh=' + str(doc.sh)+ '\n')
            f.write('kw=' + str(doc.kw)+ '\n')
            f.write('kh=' + str(doc.kh)+ '\n')
            f.write('mm=' + str(doc.mm)+ '\n')
            f.write('gross=' + str(doc.gross)+ '\n')
            f.write('gpuid=' + str(gpuid)+ '\n')
            f.write('deramp=' + str(doc.deramp)+ '\n')
            f.write('oo=' + str(doc.oversample)+ '\n')
            #f.write('nwac=' + str(doc.nwac)+ '\n')
            #f.write('nwdc=' + str(doc.nwdc)+ '\n')
            f.write('nwac=' + str(170)+ '\n')
            f.write('nwdc=' + str(1)+ '\n')
            f.write('outprefix='+prefix+ '\n')
            f.write('outsuffix='+suffix+ '\n')

            f.write('rm ' + '$outprefix$outsuffix' + '*\n')
            f.write(gpuAmpcor_driver + ' --master $master --slave $slave --ww $ww --wh $wh --sw $sw --sh $sh --mm $mm --kw $kw --kh $kh --gross $gross --outprefix $outprefix --outsuffix $outsuffix --deramp $deramp --gpuid $gpuid --nwac $nwac --nwdc $nwdc --oo $oo \n')

            f.close()

            # Check if the offset field has been calculated
            offsetfield_filename = prefix + suffix + '.bip'

            return (run_file, offsetfield_filename)


        def run_MaskAndFilter_track(self):

            doc = self.doc

            snr_threshold = 6
            filter_winsize = 8
            
            for offsetfield in self.offsetfields:
    
                date1str = offsetfield[0]
                date2str = offsetfield[1]

                offset_folder = os.path.join(self.trackfolder,self.offsetFolder,date1str+'_'+date2str)

                offset_outprefix = os.path.join(offset_folder,doc.outprefix)

                offsetfile = offset_outprefix + '_run_' + str(doc.runid) + '.bip'
                snrfile = offset_outprefix + '_run_' + str(doc.runid) + '_snr.bip'

                #print(offsetfile)
                #print(snrfile)

                run_file = offset_folder + '/mask_filter' + '_' + str(doc.runid) + '.sh'
                f = open(run_file,'w')
                f.write('#!/bin/bash\n')
                f.write('offsetfile=' + offsetfile + '\n')
                f.write('snrfile=' + snrfile+ '\n')
                f.write('snr_threshold=' + str(snr_threshold) + '\n')
                f.write('filter_winsize=' + str(filter_winsize)+ '\n')
                f.write('output_dir=' + offset_folder + '\n')

                f.write('MaskAndFilter.py -d $offsetfile -s $snrfile -n $filter_winsize -t $snr_threshold -o $output_dir' + '\n')
                f.close()

                func = globals()['run_denseoffset_bash']
                p = Process(target=func, args=('bash '+run_file,))
                self.jobs.append(p)
                p.start()
                self.doc_count = self.doc_count + 1

                # nprocess controller

                if self.doc_count == self.nproc['maskandfilter']:
                    for ip in range(self.nproc['maskandfilter']):
                        print(ip)
                        self.jobs[ip].join() #wait here

                    self.doc_count = 0
                    self.jobs = []

            return 0


        def fill_in_holes(self,data, hole_size):

            doc = self.doc

            # Find all the nan values.
            ind = np.where(np.isnan(data))
            num_nan = len(ind[0])

            allvalues = np.zeros(shape=(num_nan, hole_size*hole_size))
            holevalue = np.zeros(shape=(num_nan,))

            # Find the neighboring values of each nan.
            count = 0
            for i in range(hole_size):
                for j in range(hole_size):
                    #print(i,j)

                    # Row
                    ind_i = ind[0] + i - hole_size//2
                    ind_i[ind_i<0] = 0
                    ind_i[ind_i >= doc.numWinDown] = doc.numWinDown-1

                    # Column
                    ind_j = ind[1] + j - hole_size//2
                    ind_j[ind_j<0]=0
                    ind_j[ind_j >= doc.numWinAcross] = doc.numWinAcross-1

                    allvalues[:,count] = data[(ind_i,ind_j)] 
                    
                    count = count + 1

            # Count the number of nans in the neighborhood.
            count_nans = np.isnan(allvalues).sum(axis=1)

            # Enough supporting values in the neighborhood or not.
            ratio = 1/2
            invalid_points = count_nans > np.rint(hole_size * hole_size * ratio)

            # Find the median of neighboring values to fill in.
            holevalues = np.nanmedian(allvalues, axis=1)
    
            # Set the not well supported values back to nan.
            holevalues[invalid_points] = np.nan

            # Fill in the nans.
            data = np.copy(data)
            data[ind] = holevalues

            # Apply median filter again.
            kernel_size = (5,5)
            #data = ndimage.median_filter(input=data,size=kernel_size,mode='nearest')
            #data = signal.medfilt(data,kernel_size=kernel_size)
            data = mdf(data,kernel_size=kernel_size)
            #data = cv2.medianBlur(data,kernel_size[0])

            return data

        def _get_refer_mask(self, az, refer_az, rng, refer_rng, choice=None):

            doc = self.doc

            # Important:
            # threshold:
            # Choice 1:  hard-coded default: 1 m/d
            # Choice 2: difference and ratio to the orignal value

            if choice is None:
                choice = 2
            
            plot = False

            # Choice 1: Apply the same threshold for both az and rng
            if choice == 1:

                thrs = 1 # m/d

                # Difference m/d
                dif_az_day = az - refer_az
                dif_rng_day = rng - refer_rng

                # Select any points with at least one measurement greater than the threshold
                mask = np.logical_or(np.abs(dif_az_day)>thrs, np.abs(dif_rng_day)>thrs)

            # Choice 2: Apply difference threshold for azimuth and range. The default sets both of them to 1 m/d
            elif choice == 2:

                ########## Used for Ridgecrest Earthquake. Enforce no filtering based on reference, because reference is simply zero.
                #thrs_known_az = 1000
                #thrs_known_rng = 1000

                # Used for Glacier applications, set the threshold as 1 m/d
                thrs_known_az = 1
                thrs_known_rng = 1

                dif_az_day = az - refer_az
                dif_rng_day = rng - refer_rng

                mask = np.logical_or(np.abs(dif_az_day)>thrs_known_az,  np.abs(dif_rng_day)>thrs_known_rng)

                if plot:
                    fig = plt.figure(1,figsize=(15,10))
                    #ax = fig.add_subplot(131)
                    #ax.imshow(unknown_refer_inds,cmap='coolwarm')
                    #ax.set_title('unknown')

                    #ax = fig.add_subplot(132)
                    #ax.imshow(known_refer_inds,cmap='coolwarm')
                    #ax.set_title('known')

                    ax = fig.add_subplot(133)
                    ax.imshow(mask,cmap='coolwarm')
                    ax.set_title('final')
                    fig.savefig('mask.png')
                    print(stop)

            # Choice 3
            # For stationary points: apply for strict threshold
            # For other points, same as choice 2. The default is 1
            elif choice == 3:

                dif_az_day = az - refer_az
                dif_rng_day = rng - refer_rng

                thrs_sta_az = 0.2
                thrs_sta_rng = 0.2

                thrs_mov_az = 1
                thrs_mov_rng = 1

                ## Find stationary points ##
                # thrs_sta = 0.05
                #ind_az = np.abs(refer_az)<thrs_sta
                #ind_rng = np.abs(refer_rng)<thrs_sta
                #ind_stationary = np.logical_and(ind_az,ind_rng)

                # Updated (2020.03.13)
                ind_stationary = self._get_stationary_points(refer_az, refer_rng)

                # Pick out the stationary points with misfit larger than the threshold
                stationary_mask = np.logical_and(np.logical_or(np.abs(dif_az_day)>thrs_sta_az,  np.abs(dif_rng_day)>thrs_sta_rng), ind_stationary)

                ## Find moving points ##
                ind_moving = np.invert(ind_stationary)

                # Pick out the moving points with misfit larger than the threshold
                moving_mask = np.logical_and(np.logical_or(np.abs(dif_az_day)>thrs_mov_az,  np.abs(dif_rng_day)>thrs_mov_rng), ind_moving)

                ## Joint the bad stationary points and moving points together ###
                mask = np.logical_or(stationary_mask, moving_mask)

            # Choice 4 use ratio instead of absolute value
            elif choice == 4:

                # Difference m/d
                dif_az_day = az - refer_az
                dif_rng_day = rng - refer_rng

                # Ratio of variation.
                var_az_day = np.abs(az - refer_az) / np.max(np.abs(refer_az),0.001)
                var_rng_day = np.abs(rng - refer_rng) / np.max(np.abs(refer_rng),0.001)
    
                # Difference.
                dif_rng_day = rng - refer_rng
                dif_rng_day = rng - refer_rng

                # Control and moving points.
                thrs = 0.05
                ind_az = np.abs(refer_az)<thrs
                ind_rng = np.abs(refer_rng)<thrs

                control_points = np.logical_and(ind_az, ind_rng)

                moving_points = np.invert(control_points)

                # Mask control points by difference.
                control_points = np.logical_and(np.abs(dif_az_day)>thrs, np.abs(dif_rng_day)>thrs)

            return mask

        def _get_stationary_points(self, refer_az, refer_rng, thrs=thrs_default):

            # Find "Controlling Points". Stationary places. Smaller than thrs m/d
            # Note that np.nan values correspond to False, so they are excluded in the controlling points
            ind_az = np.abs(refer_az)<thrs
            ind_rng = np.abs(refer_rng)<thrs

            # Consider the points smaller than the threshold and validly obtained
            # ind_stationary is 2D boolean matrix
            ind_stationary = np.logical_and(ind_az,ind_rng)
           
            return ind_stationary

        def _get_misCoreg(self, az, rng, refer_az, refer_rng, ind_con_pts, order=1):

            doc = self.doc

            # Set the area that I want to exclude manually (ad hoc to Ridgecreast Earthquake)
            if self.stack == "tops_RC":
                if self.trackname == "track_64":
                    ind_con_pts[750:1700,644:] = False
                    # Only keep SW2
                    # Remove SW1
                    ind_con_pts[:, :644] = False
                    # Remove SW3
                    ind_con_pts[:, 1369:] = False

                elif self.trackname == "track_7101":
                    ind_con_pts[1600:2600, 700:1700] = False

            if postproc_verbose: 
                # Show the mask on stationary points, invalid values (NaN) are excluded.
                plt.figure(102,figsize=(10,10))
                plt.imshow(ind_con_pts)
                plt.title("The mask on stationary points, invalid values (NaN) are excluded")
                plt.savefig("102.png")

            # Creating 2D meshgrids
            ind_d = np.arange(0, doc.numWinDown)
            ind_a = np.arange(0, doc.numWinAcross)
            aa,dd = np.meshgrid(ind_a,ind_d)

            # Available points for calculating miscoregistration (1D flattened array)
            ind_d_mis = dd[ind_con_pts]
            ind_a_mis = aa[ind_con_pts]

            print('ind_d_mis: ', ind_d_mis)
            print('ind_a_mis: ', ind_a_mis)

            # Corresponding miscoregistration (1D flattened array)
            az_value_mis = az[ind_con_pts] - refer_az[ind_con_pts]
            rng_value_mis = rng[ind_con_pts] - refer_rng[ind_con_pts]

            # Set up the 2D miscoregistration matrix
            az_mat_mis = np.full(shape=(doc.numWinDown, doc.numWinAcross), fill_value=np.nan)
            az_mat_mis[ind_d_mis, ind_a_mis] = az_value_mis

            rng_mat_mis = np.full(shape=(doc.numWinDown, doc.numWinAcross), fill_value=np.nan)
            rng_mat_mis[ind_d_mis, ind_a_mis] =  rng_value_mis

            # Miscoregistration
            # Some notes on choosing the order
            # For Ridgecrest earthquake, order = 1 for A64 and D71 1km x 1 km median filter case, D71 500m x 500m median filter case
            # For Ridgecrest earthquake, order = 0 for A64, 500m x 500m median filter case

            # For glacier applications, This is not intensively testes yet.
            # Start with order=0, 2019.09.10 
            if order==0:

                az_mat_mis_pred = np.nanmedian(az[ind_con_pts] - refer_az[ind_con_pts])
                rng_mat_mis_pred = np.nanmedian(rng[ind_con_pts] - refer_rng[ind_con_pts])

                print('Constant miscoregistration correction:')
                print(az_mat_mis_pred, rng_mat_mis_pred)

            elif order ==1:
    
                # 2D polyfit
                # Azimuth
                num = len(az_value_mis)
    
                # Represent the dimension for bias G.shape = (num, 3) coordinates: [1, x, y]
                col1 = np.full(shape=(num,), fill_value=1)
                G = np.stack([col1, ind_a_mis, ind_d_mis], axis=1)
    
                (c,a,b) = np.linalg.lstsq(G, az_value_mis, rcond=None)[0]
                az_mat_mis_pred = c + a * aa + b * dd
    
                # Range
                num = len(rng_value_mis)
                col1 = np.full(shape=(num,), fill_value=1)
                G = np.stack([col1,ind_a_mis,ind_d_mis],axis=1)
    
                (c,a,b) = np.linalg.lstsq(G, rng_value_mis, rcond=None)[0]
                rng_mat_mis_pred = c + a * aa + b * dd
    
                # Plot
                plot_the_plane = False
                if plot_the_plane:
                    fig = plt.figure(200, figsize=(10,10))
                    
                    vmin = np.min(az_mat_mis_pred)
                    vmax = np.max(az_mat_mis_pred)
                    
                    ax = fig.add_subplot(121)
                    f1 = ax.imshow(az_mat_mis,cmap=plt.cm.jet,vmin=vmin,vmax=vmax)
        
                    ax = fig.add_subplot(122)
                    f1 = ax.imshow(az_mat_mis_pred,cmap=plt.cm.jet,vmin=vmin,vmax=vmax)
                    fig.colorbar(f1)
        
                    fig.savefig('200.png', format='png')

            # Save
            save_the_plane=False
            
            if doc.runid == 20190901:
                save_the_plane = True

            if save_the_plane:
                pkl_file = os.path.join(doc.offset_folder, "misCoreg_" + str(doc.runid) + '.pkl')
                misCoreg_pred = {}
                misCoreg_pred["az_mat_mis"] = az_mat_mis_pred
                misCoreg_pred["rng_mat_mis"] = rng_mat_mis_pred
                with open(pkl_file,"wb") as f:
                    pickle.dump(misCoreg_pred, f)

            return (az_mat_mis_pred,rng_mat_mis_pred)

        def _load_misCoreg(self):

            doc = self.doc

            misCoreg_pkl = os.path.join(doc.offset_folder, "misCoreg_20190901.pkl")
            with open(misCoreg_pkl,"rb") as f:
                misCoreg_pred = pickle.load(f)
                az_mat_mis_pred_low = misCoreg_pred["az_mat_mis"]
                rng_mat_mis_pred_low = misCoreg_pred["rng_mat_mis"]

                if type(az_mat_mis_pred_low)==np.float64:
                    az_mat_mis_pred, rng_mat_mis_pred = az_mat_mis_pred_low, rng_mat_mis_pred_low
                else:
                    row, col = az_mat_mis_pred_low.shape
                    zoom_ratio = max(doc.numWinDown/row, doc.numWinAcross/col)
                    az_mat_mis_pred = ndimage.zoom(az_mat_mis_pred_low, zoom=zoom_ratio, order=1, mode="nearest")
                    rng_mat_mis_pred = ndimage.zoom(rng_mat_mis_pred_low, zoom=zoom_ratio, order=1, mode="nearest")

                    new_row, new_col = az_mat_mis_pred.shape

                    if new_row > doc.numWinDown:
                        start_ind = (new_row - doc.numWinDown)//2
                        az_mat_mis_pred = az_mat_mis_pred[start_ind:start_ind+doc.numWinDown, :]
                        rng_mat_mis_pred = rng_mat_mis_pred[start_ind:start_ind+doc.numWinDown, :]

                    if new_col > doc.numWinAcross:
                        start_ind = (new_col - doc.numWinAcross)//2
                        az_mat_mis_pred = az_mat_mis_pred[:, start_ind:start_ind + doc.numWinAcross]
                        rng_mat_mis_pred = rng_mat_mis_pred[:, start_ind:start_ind + doc.numWinAcross]

                    #print(az_mat_mis_pred.shape, rng_mat_mis_pred.shape)
                    #print(doc.numWinDown, doc.numWinAcross)

            return az_mat_mis_pred, rng_mat_mis_pred


        def _get_median_filter_window(self):

            doc = self.doc

            # Used for Ridgecreast earthquake 1kmx1km
            #med2_kernel_size = (9,15)

            # (5,7)
            # Used for Ridgecreast earthquake 500mx500m
            # Used for CSK-Rutford 128 x 128 window
            # Size: range x azimuth
            if doc.runid==20190901:
                return (5,5)
            # winsize = 240m x 240m, step 120m x 120m
            elif doc.runid ==20190904:
                return (5,5)

            # winsize 120m x 240 m, step 60 m x 120 m
            elif doc.runid==20190908:
                return (5,3)

            # winsize 120m x 120m, step size 60 m x 60 m
            elif doc.runid == 20190921:
                return (5,5)

            # winsize 60m x 120m, step size 30 m x 60 m
            elif doc.runid == 20190925:
                return (9,5)

            # winsize 500m x 1800m step size 250 m x 900 m 
            elif doc.runid == 20200101:
                return (9,5)

            # winsize 1000m x 1800m step size 500 m x 900 m
            elif doc.runid == 20200102:
                return (9,5)

            else:
                raise Exception("Need size of median filter")

        def _get_nan_expand_window_size(self):

            doc = self.doc

            if doc.runid==20190901:
                return 7
            elif doc.runid==20190904:
                return 7
            elif doc.runid==20190908:
                return 7

            # winsize 120m x 120m, step size 60 m x 60 m
            elif doc.runid == 20190921:
                return 7

            # winsize 60m x 120m, step size 30 m x 60 m
            elif doc.runid == 20190925:
                return 7
            elif doc.runid == 20200101:
                return 7
            elif doc.runid == 20200102:
                return 7
            else:
                raise Exception("Need size of Nan expand window")

        def _generic_median_filter(self, data, data_snr, mask=None, label=None, refer=None):

            doc = self.doc

            ## Get the window sizes
            median_filter_window = self._get_median_filter_window()
            nan_expand_window_size = self._get_nan_expand_window_size()
           
            ## Start ###
            med1_data = np.copy(data)

            # Step 0: mask out using snr
            do_step_0 = False
            
            if do_step_0:
                med1_data[data_snr<4] = np.nan

            # Step 1:
            do_step_1 = True
            if do_step_1:
                # Mask out where the reference is nan.
                med1_data[np.isnan(refer)] = np.nan
                
                # Using the mask, set masked value to nan
                if mask is not None: 
                    med1_data[mask] = np.nan

            # Step 2: normal median filter
            do_step_2 = True
            if do_step_2:
                med2_data = np.copy(med1_data)
                
                med2_kernel_size = median_filter_window
                #print(med2_kernel_size)

                # Median filters
                # option 1:
                #med2_data = ndimage.median_filter(input=med2_data,size=med2_kernel_size,mode='nearest')
                
                # option 2:
                #med2_data = signal.medfilt(med1_data,kernel_size=med2_kernel_size)
                
                # option 3:
                med2_data = mdf(med1_data, kernel_size=med2_kernel_size)

                # option 4:
                #med2_data = cv2.medianBlur(med1_data,med2_kernel_size[0])

            else:
                med2_data = np.copy(med1_data)

            # Step 3: iteratively fill in small holes, fill in + median filter
            do_step_3 = False

            if do_step_3:
                iteration = 3
                hole_size = 7
                med3_data = np.copy(med2_data)
                for i in range(iteration):
                    med3_data = self.fill_in_holes(data=med3_data, hole_size = hole_size)
                    #med3_data = self.fill(med3_data)

            else:
                med3_data = np.copy(med2_data)


            ## Step 4: set neighboring values of nans as nan.
            do_step_4 = True
            if do_step_4:
    
                med4_kernel_size = nan_expand_window_size

                med4_data = np.copy(med3_data)
                ind = np.where(np.isnan(med4_data))
    
                for i in range(med4_kernel_size):
                    for j in range(med4_kernel_size):
                        ind_i = ind[0] + i - med4_kernel_size//2
                        ind_i[ind_i < 0] = 0
                        ind_i[ind_i >= doc.numWinDown]=doc.numWinDown-1
    
                        ind_j = ind[1] + j - med4_kernel_size//2
                        ind_j[ind_j < 0] = 0
                        ind_j[ind_j >= doc.numWinAcross]=doc.numWinAcross-1
    
                        med4_data[(ind_i,ind_j)] = np.nan
            else:
                med4_data = np.copy(med3_data)

            # Step 5: remove the values on the margin, if interpolation is performed.
            if do_step_3:
                bb = (iteration + 1) * (hole_size // 2)
                med4_data[:iteration,:] = np.nan
                med4_data[-iteration:,:] = np.nan
                med4_data[:,:iteration] = np.nan
                med4_data[:,-iteration:] = np.nan

            return med4_data

        def _display_az_rng(self, azOff, rngOff, azOff_re, rngOff_re, azOff_cmp, rngOff_cmp, title= None):

            doc = self.doc
            pad = 0.2

            # Ranges of values.
            # Plot azimuth offset.
            vmin = np.nanmin(azOff_re)
            vmax = np.nanmax(azOff_re)

            vmin = - max(abs(vmin), abs(vmax))
            vmax =   max(abs(vmin), abs(vmax))

            vmin10 = np.floor(vmin*10)/10-pad
            vmax10 = np.ceil(vmax*10)/10+pad

            fig = plt.figure(figsize=(18,6))

            frac = 0.07
            padbar = 0.1
            tickstep=0.4

            ax = fig.add_subplot(161)
            ax.set_title('Predicted azimuth offset') 
            im = ax.imshow(azOff_re, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            fig.colorbar(im,fraction=frac, pad=padbar, orientation='horizontal',ticks=np.arange(np.round(vmin10),np.round(vmax10)+1,tickstep),label='m')
            #fig.colorbar(im,fraction=0.07, orientation='horizontal',label='meter')

            ax = fig.add_subplot(162)
            ax.set_title('Raw azimuth offset') 
            im = ax.imshow(azOff_cmp, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            fig.colorbar(im,fraction=frac, pad=padbar, orientation='horizontal',ticks=np.arange(np.round(vmin10),np.round(vmax10)+1,tickstep),label='m')

            ax = fig.add_subplot(163)
            ax.set_title('Filtered azimuth offset')
            im = ax.imshow(azOff, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            fig.colorbar(im,fraction=frac, pad=padbar, orientation='horizontal',ticks=np.arange(np.round(vmin10),np.round(vmax10)+1,tickstep),label='m')

            # Plot range offset.
            vmin = np.nanmin(rngOff_re)
            vmax = np.nanmax(rngOff_re)

            vmin = - max(abs(vmin), abs(vmax))
            vmax =   max(abs(vmin), abs(vmax))

            vmin10 = np.floor(vmin*10)/10-pad
            vmax10 = np.ceil(vmax*10)/10+pad

            ax = fig.add_subplot(164)
            ax.set_title('Predicted range offset')
            im = ax.imshow(rngOff_re, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            fig.colorbar(im,fraction=frac, pad=padbar, orientation='horizontal',ticks=np.arange(np.round(vmin10),np.round(vmax10)+1,tickstep),label='m')

            ax = fig.add_subplot(165)
            ax.set_title('Raw range offset') 
            im = ax.imshow(rngOff_cmp, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            fig.colorbar(im,fraction=frac, pad=padbar, orientation='horizontal',ticks=np.arange(np.round(vmin10),np.round(vmax10)+1,tickstep),label='m')

            ax = fig.add_subplot(166)
            ax.set_title('Filtered range offset')          
            im = ax.imshow(rngOff, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            fig.colorbar(im,fraction=frac, pad=padbar, orientation='horizontal',ticks=np.arange(np.round(vmin10),np.round(vmax10)+1,tickstep),label='m')

            # Save the figure
            figdir = os.path.join(os.path.join(self.workdir,'figs', self.trackname,'radar'))
            try:
                if not os.path.exists(figdir):
                    os.makedirs(figdir)
            except:
                pass

            # PDF
            fig.savefig(os.path.join(figdir,'offset_' + title + ".pdf"), format='pdf',bbox_inches='tight')
            # PNG
            fig.savefig(os.path.join(figdir,'offset_' + title + ".png"), format='png',bbox_inches='tight')

            # Cut a line through to find the boundary of Swath, for Ridgecrest S1 only.
            #fig  = plt.figure(2, figsize=(10,10))
            #ax = fig.add_subplot(111)
            #line = azOff[2000,:]
            #ax.plot(line)
            #fig.savefig(figdir + '/'+'line.png')
            #print(np.arange(len(line))[line<-8])

            plt.close(fig)

            return 0

        def _save_az_rng(self, azOffset_filtered, azOffsetName, rngOffset_filtered, rngOffsetName, final_mask=None):

            doc = self.doc

            # Convert back to pixels in real interim days.
            azOffset_filtered_ori = azOffset_filtered / doc.azPixelSize * doc.interim
            rngOffset_filtered_ori = rngOffset_filtered / doc.rngPixelSize * doc.interim

            redo = 1

            # Remove the existed ones.
            os.system('rm -r ' + doc.offset_folder + '/filt*' + self.filt_suffix + '.off')

            ######################
            # Save azimuth offset.

            nbands = 1
            bandType = 6

            if not os.path.exists(azOffsetName) or redo == 1:
               
                driver = gdal.GetDriverByName( 'ENVI' )
                dst_ds = driver.Create(azOffsetName, azOffset_filtered_ori.shape[1], azOffset_filtered_ori.shape[0], nbands, bandType )
                dst_ds.GetRasterBand(1).WriteArray( azOffset_filtered_ori, 0 ,0 )
                dst_ds = None

                outImg = isceobj.createImage()
                outImg.setDataType('Float')
                outImg.setFilename(os.path.join(doc.offset_folder, azOffsetName))
                outImg.setBands(1)
                outImg.scheme = 'BIP'
                outImg.setWidth(doc.numWinAcross)
                outImg.setLength(doc.numWinDown)
                outImg.setAccessMode('read')
                outImg.renderHdr()
                outImg.renderVRT()
           
            ####################
            # Save range offset.

            nbands = 1
            bandType = 6

            if not os.path.exists(rngOffsetName) or redo == 1:
            
                driver = gdal.GetDriverByName( 'ENVI' )
                dst_ds = driver.Create(rngOffsetName, rngOffset_filtered_ori.shape[1], rngOffset_filtered_ori.shape[0], nbands, bandType )
                dst_ds.GetRasterBand(1).WriteArray( rngOffset_filtered_ori, 0 ,0 )
                dst_ds = None


                outImg = isceobj.createImage()
                outImg.setDataType('Float')
                outImg.setFilename(os.path.join(doc.offset_folder, rngOffsetName))
                outImg.setBands(1)
                outImg.scheme = 'BIP'
                outImg.setWidth(doc.numWinAcross)
                outImg.setLength(doc.numWinDown)
                outImg.setAccessMode('read')
                outImg.renderHdr()
                outImg.renderVRT()

            return 0

        #################### 2019.09.28 Developing for better filtering ####################
        def _offset_filter_v2(self, offsetfield, offsetfile, snrfile, offsetLosFile, title, redo=True):

            # Filter the offset fields.
            # The intermediate unit is meter/day.
            # The final unit is still in pixel.
            doc = self.doc

            if self.rb_suffix:
                self.filt_suffix = self.filt_suffix + '_' + self.rb_suffix
            else:
                self.filt_suffix = self.filt_suffix

            azOffsetName = os.path.join(doc.offset_folder, 'filtAzimuth_' + self.filt_suffix  + '.off')

            rngOffsetName = os.path.join(doc.offset_folder, 'filtRange_' + self.filt_suffix + '.off')

            redo_this = redo
            if os.path.exists(azOffsetName) and os.path.exists(rngOffsetName) and redo_this==False:
                print("filtered offsetfield files exist")
                return

            # Read in raw offset fields
            ds = gdal.Open(offsetfile)
            azOffset = ds.GetRasterBand(1).ReadAsArray()
            rngOffset = ds.GetRasterBand(2).ReadAsArray()

            # Read in SNR
            ds = gdal.Open(snrfile)
            data_snr = ds.GetRasterBand(1).ReadAsArray()

            # Read in incidence angle
            ds = gdal.Open(offsetLosFile)
            inc = ds.GetRasterBand(1).ReadAsArray()

            # Generatre a mask for invalid values at margins (mostly useful for S1, GPU ampcor gives (-4, -4) to invalid cross-correlation))
            # Alternative way is to use incidence angle in offsetLosFile, but could be wrong if offsetLosFile doesn't match offset field
            # Mask: (True: invalid, False: valid)
            mask_of_invalid = np.logical_or(np.logical_and(azOffset==-4, rngOffset==-4), inc==0)

            # Convert the unit from pixel to meter per day
            azOffset = azOffset * doc.azPixelSize / doc.interim
            rngOffset = rngOffset * doc.rngPixelSize / doc.interim

            # Get reference offset field
            refer_azOffset, refer_rngOffset, vmin_az, vmax_az, vmin_rng, vmax_rng = self._obtain_reference()
 
            # Make a copy of the offset fields for later use
            azOffset_copy = azOffset.copy()
            rngOffset_copy = rngOffset.copy()

            ###  Filtering ###
            # Generate mask by comparison to reference. Values deviate too much from reference set as True. Use for deriving mis-coreg.
            mask_by_reference = self._get_refer_mask(azOffset, refer_azOffset, rngOffset, refer_rngOffset, choice=2)
            # Combine the mask of both invalid values and erroreous estimations
            mask = np.logical_or(mask_of_invalid, mask_by_reference)

            if postproc_verbose:
                # Show the mask based on deviation from reference
                plt.figure(2, figsize=(10,10))
                plt.imshow(mask_by_reference)
                plt.title("The first mask: mask based on deviation from reference offset fields")
                plt.savefig('2.png')
     
                plt.figure(3, figsize=(10,10))
                plt.imshow(mask)
                plt.title("Combined mask of invalid values (-4, -4) and values of large deviations")
                plt.savefig('3.png')
    
                # Show the original unfiltered and uncoregistered azimuth offset
                fig = plt.figure(4, figsize=(10,10))
                ax1 = fig.add_subplot(121)
                im1 = ax1.imshow(azOffset,vmin=vmin_az, vmax=vmax_az, cmap=cm.jet)
                fig.colorbar(im1)
    
                ax2 = fig.add_subplot(122)
                im2 = ax2.imshow(rngOffset,vmin=vmin_rng, vmax=vmax_rng, cmap=cm.jet)
                fig.colorbar(im2)
    
                plt.title("Original Ampcor azimuth offset and range offset")
                plt.savefig('4.png')

            ############ Run filtering ###############################
            # 1) set the masked value to np.nan and 2) run 2-D median filter across the fields
            filtered_azOffset = self._generic_median_filter(azOffset, data_snr, mask=mask, label='az',refer=refer_azOffset)
            filtered_rngOffset = self._generic_median_filter(rngOffset, data_snr, mask=mask, label='rng',refer=refer_rngOffset)

            # Manually remove and first and second swath (ad hoc to Ridgecreast Earthquake)
            if self.stack == "tops_RC":
                if self.trackname == "track_64":
                    #filtered_azOffset[:,:644] = np.nan
                    #filtered_azOffset[:,1369:] = np.nan
                    #filtered_rngOffset[:,:644] = np.nan
                    #filtered_rngOffset[:,1369:] = np.nan
                    pass

                elif self.trackname == "track_7101":
                    pass

            if postproc_verbose:
                # Show the filtered range and azimuth offset
                fig = plt.figure(5, figsize=(10,10))
                ax1 = fig.add_subplot(121)
                im1 = ax1.imshow(filtered_azOffset,vmin=vmin_az, vmax=vmax_az, cmap=cm.jet)
                fig.colorbar(im1)
    
                ax2 = fig.add_subplot(122)
                im2 = ax2.imshow(filtered_rngOffset,vmin=vmin_rng, vmax=vmax_rng, cmap=cm.jet)
                fig.colorbar(im2)
                
                plt.title('Filtered but uncoregistered')
                plt.savefig('5.png')


            ## Miscoregistration estimation
            # Step 1: Update mask for invalid values
            # So far, mask_by_reference only contains bad value from get_refer_mask using choice = 2
            # Mask has mask_by_reference plus the invalid (-4, -4) points
            # The median filter may introduce more nan values into the offset fields.
            # The NaN in azOffset filtered is a superset of Mask and Mask_by_reference
            # Include them into mask_by_reference

            # Old code:
            mask_by_reference[np.isnan(filtered_azOffset)] = True
            mask = mask_by_reference

            # New code
            #mask = np.isnan(filtered_azOffset)

            if postproc_verbose:
                # Show the mask corresponding the Nan value before miscoregistration correction
                plt.figure(figsize=(10,10))
                im = plt.imshow(mask)
                plt.title("NaN value before miscoregistration correction")
                plt.savefig('9.png')

#            # Step 2: Estimate how complex the miscoregistration can be based on the valid control points
#            good_ratio_seg1, good_ratio_seg2 = self._check_filtered_offsetfield(title, filtered_azOffset, filtered_rngOffset, refer_azOffset, refer_rngOffset)
#
#            if good_ratio_seg1>=1/3 and good_ratio_seg2>=1/3:
#                misCoreg_order = 1
#            else:
#                misCoreg_order = 0
#
#            # Adjustment
#            if self.version=="v10" or self.version=="v101" or self.version=="v11":
#                misCoreg_order = 0
#
#            if doc.runid!=20190901 and self.version=="v12":
#                external_misCoreg = True
#
#            # Only version 12 uses miscoregistration derived from lower resolution results
#            else:
#                external_misCoreg = False

            external_misCoreg = False
            misCoreg_order = 0
            # Step 3: Get miscoregistration
            if external_misCoreg == False:
                # Only non-masked values (non nan value, valid estimate) should be used for calculating miscoregistration
                # Find "Controlling Points". Stationary places. Smaller than thrs m/d (Update 2020.03.14)
                ind_stationary = self._get_stationary_points(refer_azOffset, refer_rngOffset)
                # Remve invalid values from controlling points
                ind_control_points = np.logical_and(ind_stationary, np.invert(mask))
                # Calucate the micoregistration
                az_mis, rng_mis = self._get_misCoreg(filtered_azOffset, filtered_rngOffset, refer_azOffset, refer_rngOffset, ind_control_points, order = misCoreg_order)
            else:
                print(stop)
                az_mis, rng_mis = self._load_misCoreg()

            print('Estimated miscoregistration (azimuth & range): ', az_mis, rng_mis)

            ### End of micoregistration estimation   

            # For rubbersheeting, there is no miscoregistration
            if self.run_mode == "pair_rb":
                print("pair rb mode, no miscoregistration correction")
                miscoreg_correction = False

                # record the miscoregistration
                self.az_mis = az_mis
                self.rng_mis = rng_mis

            else:
                ################## Start a new round of filtering ############################
                # Mis-coregistration correction.

                miscoreg_correction = True
                
                if miscoreg_correction:
                    if self.stack == "stripmap":
                        print("Perform miscoregistration correction")
                        azOffset = azOffset_copy - az_mis
                        rngOffset = rngOffset_copy - rng_mis

                    elif self.stack == "tops":
                        print("Skip miscoregistration")
                        azOffset = azOffset_copy
                        rngOffset = rngOffset_copy

                # Get new mask by reference (choice = 3!)
                mask_by_reference = self._get_refer_mask(azOffset, refer_azOffset, rngOffset, refer_rngOffset, choice=3)
                # Combine the mask of both invalid values and erroreous estimations
                mask = np.logical_or(mask_of_invalid, mask_by_reference)
    
                # Show the new mask
                if postproc_verbose:
                    # Show the mask based on deviation from reference
                    plt.figure(12, figsize=(10,10))
                    plt.imshow(mask_by_reference)
                    plt.title("The first mask after miscoreg correction: mask based on deviation from reference offset fields (strict or stationary points)")
                    plt.savefig('12.png')
         
                    plt.figure(13, figsize=(10,10))
                    plt.imshow(mask)
                    plt.title("Combined mask of invalid values (-4, -4) and values of large deviations after micoregistration correction")
                    plt.savefig('13.png')
     
                # Filtering
                filtered_azOffset = self._generic_median_filter(azOffset, data_snr, mask=mask, label='az',refer=refer_azOffset)
                filtered_rngOffset = self._generic_median_filter(rngOffset, data_snr, mask=mask, label='rng',refer=refer_rngOffset)
    
                if postproc_verbose:
                    # Show the filtered and coregistered offset fields
                    fig = plt.figure(10, figsize=(10,10))
                    ax1 = fig.add_subplot(121)
                    im1 = ax1.imshow(filtered_azOffset,vmin=vmin_az, vmax=vmax_az, cmap=cm.jet)
                    fig.colorbar(im1)
        
                    ax2 = fig.add_subplot(122)
                    im2 = ax2.imshow(filtered_rngOffset,vmin=vmin_rng, vmax=vmax_rng, cmap=cm.jet)
                    fig.colorbar(im2)
                    
                    plt.title('Filtered and coregistered')
                    plt.savefig('10.png')

            ################## End of the second round of filtering

            ############# SAVE to ENVI files ############################
            # Convert the unit from meter per day to pixel.
            # Save the offset fields.
            self._save_az_rng(filtered_azOffset, azOffsetName, filtered_rngOffset, rngOffsetName)

            ############ Plot the results ##############################
            #self._plot_offsetfield(offsetfield)
            
            return 0

        def offset_filter(self, s_date=None, e_date=None):

            doc = self.doc

            # Perform segmentation
            self._segmentation()

            for count, offsetfield in enumerate(self.offsetfields):

                date1str, date2str, date1, date2 = offsetfield

                # Skip the dates outside of the range.
                if s_date is not None and date1 < s_date.date():
                    continue

                if e_date is not None and date2 > e_date.date():
                    continue

                doc.interim = (date2-date1).days

                doc.offset_folder = os.path.join(self.trackfolder,self.offsetFolder, date1str+'_'+date2str)

                offset_outprefix = os.path.join(doc.offset_folder,doc.outprefix)

                # Prepare the corresponding files.
                offsetfile = offset_outprefix + '_run_' + str(doc.runid) + '.bip'
                snrfile = offset_outprefix + '_run_' + str(doc.runid) + '_snr.bip'

                offsetLosFile = doc.offsetLosFile

                # Check if the files exist
                if not os.path.exists(offsetfile) or not os.path.exists(snrfile):
                    print("The offset field file doesn't exist")
                    print("skip ", date1str+'_'+date2str)
                    continue
                else:
                    print("Offset field file exists")

                title = date1str + '_' + date2str

                # check if we should do it
                exist = False

                if exist == False:
                    print('Should work on', title)
                    
                    if self.exe == True:

                        print('work on it!')
                        
                        func = self._offset_filter_v2

                        p = Process(target=func, args=(offsetfield, offsetfile, snrfile, offsetLosFile, title,))
                        self.jobs.append(p)
                        p.start()
                        self.doc_count = self.doc_count + 1

                        # nprocess controller
                        if self.doc_count == self.nproc['maskandfilter'] or count==len(self.offsetfields)-1:
                            #for ip in range(self.nproc['maskandfilter']):
                            for ip in range(len(self.jobs)):
                                print(ip)
                                self.jobs[ip].join() #wait here

                            self.doc_count = 0
                            self.jobs = []

            return 0

        def offset_filter_pair_iter(self, iter_num, exe=False, redo=True):

            doc = self.doc

            date1str, date2str, date1, date2 = self.offsetfields[0]

            doc.interim = (date2-date1).days

            doc.offset_folder = os.path.join(self.runfolder, self.offsetFolder, date1str+'_'+date2str)

            # Prefix & suffix
            rb_suffix = self.rb_suffix 
            offset_outprefix = os.path.join(doc.offset_folder,doc.outprefix)
            offset_outsuffix = '_run_' + str(doc.runid)+ '_' + rb_suffix

            # Prepare the corresponding files.
            offsetfile = offset_outprefix + offset_outsuffix + '.bip'
            snrfile = offset_outprefix + offset_outsuffix + '_snr.bip'
            offsetLosFile = doc.offsetLosFile

            # Check if the files exist
            if not os.path.exists(offsetfile) or not os.path.exists(snrfile):
                print("The offset field file doesn't exist")
                print("skip ", date1str+'_'+date2str)
                return

            title = date1str + '_' + date2str

            # check if we should do it
            exist = False

            if exist == False:
                print('Should work on', title)
                if exe == True:
                    print('work on it!')
                    self._offset_filter_v2(offsetfile, snrfile, offsetLosFile, title, redo=redo)

            return 0


        def postprocess(self, s_date=None, e_date=None):

            self.offset_filter(s_date,e_date)
            
            return 0

        ############### End of offset filter #################################

        ## For checking the filtered offset field ##### 2020.03.10 ############
        def _check_filtered_offsetfield(self, title, filtered_azOffset=None, filtered_rngOffset=None, refer_azOffset=None, refer_rngOffset=None):

            # Check the filtered offset field
            # The unit is meter/day.
            doc = self.doc
            stationary_segments = doc.stationary_segments

            if filtered_azOffset is None or filtered_rngOffset is None:

                # Check filtered offset field files
                if self.rb_suffix:
                    self.filt_suffix = self.filt_suffix + '_' + self.rb_suffix
                else:
                    self.filt_suffix = self.filt_suffix
    
                azOffsetName = os.path.join(doc.offset_folder, 'filtAzimuth_' + self.filt_suffix  + '.off')
                print('azOffsetName: ', azOffsetName)
    
                rngOffsetName = os.path.join(doc.offset_folder, 'filtRange_' + self.filt_suffix + '.off')
                print('rngOffsetName: ',rngOffsetName)
    
                if not os.path.exists(azOffsetName) or not os.path.exists(rngOffsetName):
                    print("filtered offsetfield doesn't exist")
                    return
                else:
                    print("filtered offsetfield exists")
    
                # Read in filtered offset fields
                # Note that the offset field contains invalid values
                ds = gdal.Open(azOffsetName)
                filtered_azOffset = ds.GetRasterBand(1).ReadAsArray()
                ds = gdal.Open(rngOffsetName)
                filtered_rngOffset = ds.GetRasterBand(1).ReadAsArray()
    
                # Convert observation from pixel to meter per day.
                filtered_azOffset = filtered_azOffset * doc.azPixelSize / doc.interim
                filtered_rngOffset = filtered_rngOffset * doc.rngPixelSize / doc.interim
    
            # Get valid estimate from offset field
            ind_invalid_estimate = np.isnan(filtered_azOffset)
            ind_valid_estimate = np.invert(ind_invalid_estimate)

            # Reference
            if refer_azOffset is None or refer_rngOffset is None:
                # Find reference offset field
                refer_azOffset, refer_rngOffset, _ , _ , _ , _ = self._obtain_reference()

            # Get stationary points from reference
            ind_stationary = self._get_stationary_points(refer_azOffset, refer_rngOffset)

            # Find the maximum and minimum of the measurement on the stationary points
            #max_az_offset = np.nanmax(filtered_azOffset[ind_stationary])
            #min_az_offset = np.nanmin(filtered_azOffset[ind_stationary])
            #max_rng_offset = np.nanmax(filtered_rngOffset[ind_stationary])
            #min_rng_offset = np.nanmin(filtered_rngOffset[ind_stationary])
            #print(max_az_offset, min_az_offset, max_rng_offset, min_rng_offset)

            # Find and show control points: stationary and valid
            ind_control_points = np.logical_and(ind_stationary, ind_valid_estimate)
            schematic = np.zeros(shape=filtered_azOffset.shape)
            schematic[ind_valid_estimate] = 2
            schematic[ind_control_points] = 1
             
            # Calculate statistics
            # Segment 1
            total_num_seg1 = np.nansum(stationary_segments==1)
            total_num_contp_seg1 = np.nansum(np.logical_and(stationary_segments==1, schematic==1))
            good_ratio_seg1  = round(total_num_contp_seg1 / total_num_seg1,2)
            
            # Segment 2
            total_num_seg2 = np.nansum(stationary_segments==2)
            total_num_contp_seg2 = np.nansum(np.logical_and(stationary_segments==2, schematic==1))
            good_ratio_seg2  = round(total_num_contp_seg2 / total_num_seg2,2)

            print("M2")
            print(total_num_seg1, total_num_seg2)
            print(total_num_contp_seg1, total_num_contp_seg2)

            if check_offset_verbose:
                fig = plt.figure(1, figsize=(10,10))
                ax = fig.add_subplot(111)
                ax.imshow(schematic)
                ax.set_title(" ".join(["seg1:",str(good_ratio_seg1),"  seg2: ",str(good_ratio_seg2)]))
                fig_name = os.path.join(disp_temp_folder, title + "_" + str(doc.runid) + ".png")
                fig.savefig(fig_name)
 
            return good_ratio_seg1, good_ratio_seg2

        def _obtain_reference(self):

            doc = self.doc

            # Reference
            # For glacier applications, load the predicted offset fields.
            if self.stack in ["tops","stripmap"] and (self.run_mode=="stack" or self.run_mode=="pair_rb" and self.rb_iter_num==0):
                refer_azOffset = doc.grossDown
                refer_rngOffset = doc.grossAcross
            # For tectonic geodesy (e.g. Ridgecrest Earthquake), the reference is simply zero.
            else:
                refer_azOffset = np.zeros(shape=azOffset.shape)
                refer_rngOffset = np.zeros(shape=rngOffset.shape)

            # Convert reference from "pixel per day" to "meter per day".
            refer_azOffset = refer_azOffset * doc.azPixelSize
            refer_rngOffset = refer_rngOffset * doc.rngPixelSize

            # Get vmin and vmax
            vmin_az = np.nanmin(refer_azOffset)
            vmax_az = np.nanmax(refer_azOffset)

            vmin_az = - max(abs(vmin_az), abs(vmax_az))
            vmax_az =   max(abs(vmin_az), abs(vmax_az))

            vmin_rng = np.nanmin(refer_rngOffset)
            vmax_rng = np.nanmax(refer_rngOffset)

            vmin_rng = - max(abs(vmin_rng), abs(vmax_rng))
            vmax_rng =   max(abs(vmin_rng), abs(vmax_rng))

            return refer_azOffset, refer_rngOffset, vmin_az, vmax_az, vmin_rng, vmax_rng

        def _segmentation(self):

            # use KMeans for clustering
            from sklearn.cluster import KMeans

            doc = self.doc

            # Get reference offset field
            refer_azOffset, refer_rngOffset, _, _ ,_ ,_ = self._obtain_reference()

            # Get stationary points from reference
            ind_stationary = self._get_stationary_points(refer_azOffset, refer_rngOffset)

            # Get the shape of image
            nRow, nCol = ind_stationary.shape

            # Segmentation of stationary points
            X = [[i,j] for i in range(nRow) for j in range(nCol) if ind_stationary[i,j]]

            # KMeans clustering
            random_state = 1
            #kmeans = KMeans(n_clusters = track_to_numSeg[self.trackname], random_state=random_state, max_iter=600).fit(X)
            kmeans = KMeans(n_clusters = 2, random_state=random_state, max_iter=300).fit(X)

            # Save the segmentation
            stationary_segments = np.zeros(shape=ind_stationary.shape)
            for ix in range(len(X)):
                stationary_segments[X[ix][0], X[ix][1]] = kmeans.labels_[ix] + 1

            # Show the segmentation
            fig = plt.figure(100, figsize=(10,10))
            ax = fig.add_subplot(111)
            ax.imshow(stationary_segments)
            figName = disp_temp_folder + '/' + "track_" + str(self.tracknumber).zfill(3) + '_' + str(random_state) + '_segments.png'
            fig.savefig(figName)

            # Save it to the offset field object
            doc.stationary_segments = stationary_segments

            return 0

        def run_check_offset_track(self, s_date=None, e_date=None):

            doc = self.doc

            # Perform segmentation
            self._segmentation()

            count = 0
            for offsetfield in self.offsetfields:

                count = count + 1

                date1str = offsetfield[0]
                date2str = offsetfield[1]

                date1 = offsetfield[2]
                date2 = offsetfield[3]

                # Skip the dates outside of the range.
                if s_date is not None and date1 < s_date.date():
                    continue

                if e_date is not None and date2 > e_date.date():
                    continue

                doc.interim = (date2-date1).days

                doc.offset_folder = os.path.join(self.trackfolder, self.offsetFolder, date1str+'_'+date2str)

                title = date1str + '_' + date2str

                # check if we should do it
                exist = False

                if exist == False:

                    print('Should work on', title)
                    
                    # ad hoc control
                    #if self.exe == True or title == '20170613_20170625':
                    #if self.exe == True or title == '20171030_20171111':
                    #if self.exe == True or count == 1:

                    if self.exe == True:

                        print('work on it!')
                        
                        func = self._check_filtered_offsetfield

                        p = Process(target=func, args=(title,))
                        self.jobs.append(p)
                        p.start()
                        self.doc_count = self.doc_count + 1

                        # nprocess controller
                        if self.doc_count == self.nproc or count==len(self.offsetfields):
                            for ip in range(len(self.jobs)):
                                print(ip)
                                self.jobs[ip].join() #wait here

                            self.doc_count = 0
                            self.jobs = []

            return 0

        def _get_bbox(self):

            doc = self.doc

            latfile = doc.offsetLatFile
            lonfile = doc.offsetLonFile

            ds=gdal.Open(latfile)
            b=ds.GetRasterBand(1)
            data=b.ReadAsArray()
            
            minlat = np.min(data[np.nonzero(data)])
            maxlat = np.max(data[np.nonzero(data)])
            
            ds=gdal.Open(lonfile)
            b=ds.GetRasterBand(1)
            data=b.ReadAsArray()
         
            minlon = np.min(data[np.nonzero(data)])
            maxlon = np.max(data[np.nonzero(data)])

            # Lat is normal, Lon is scaled to be smaller
            # Lat step is normal, Lon step should be larger
            # Rounding of Lat should be finer than Lon
            # coe_lat > coe_lon
            # The larger coe is, the finer rounding is

            if self.stack in ["tops","stripmap"]:
                # Round lon to 0.1 degree
                #coe_lon = 2
                coe_lon = 10

                # Round lat to 0.02 degree
                #coe_lat = 10
                coe_lat = 50
            else:
                
                coe_lon = 10
                coe_lat = 10
                #raise Exception("Need lat lon resolution")
        
            minlat = np.ceil(minlat * coe_lat)/coe_lat
            minlon = np.ceil(minlon * coe_lon)/coe_lon
        
            maxlat = np.floor(maxlat * coe_lat)/coe_lat
            maxlon = np.floor(maxlon * coe_lon)/coe_lon
        
            SNWE = str(minlat) + " " + str(maxlat) + " " + str(minlon) + " " +str(maxlon)            
            
            return SNWE

        def _get_trim_size(self):

            doc = self.doc
            if doc.runid==20190901:
                return 45
            elif doc.runid == 20190904:
                return 25
            elif doc.runid == 20190921:
                return 15
            else:
                raise Exception("Trim size undefined")

        def _trim_borders(self, filename):

            trim_size = self._get_trim_size()

            ds = gdal.Open(filename)
            print("filename: ", filename)
            data_b1 = ds.GetRasterBand(1).ReadAsArray()
            data_b2 = ds.GetRasterBand(2).ReadAsArray()

            # Trim the value based on nodata, which is 0 by default

            # Find the borders
            nRow, nCol = data_b1.shape
            mask = np.zeros(shape=data_b1.shape)
            for j in range(nCol):
                data_idx = np.where(data_b1[:,j]!=0)
                first_idx = data_idx[0][0]
                last_idx = data_idx[0][-1]
                print(first_idx, last_idx)
                mask[first_idx : first_idx + trim_size, j] = 1
                mask[last_idx - trim_size + 1 : last_idx+1, j] = 1

            # Remove the borders
            data_b1[mask==1] = 0
            data_b2[mask==1] = 0

            # Write data
            dirname = os.path.dirname(filename)
            basefilename, ending = os.path.basename(filename).split(".")
            new_filename = os.path.join(dirname, basefilename+'_new'+"."+ending)
            print("new_filename: ", new_filename)

            nbands = 2
            bandType = 6

            driver = gdal.GetDriverByName( 'ENVI' )
            dst_ds = driver.Create(filename, data_b1.shape[1], data_b2.shape[0], nbands, bandType )
            dst_ds.GetRasterBand(1).WriteArray( data_b1, 0 ,0 )
            dst_ds.GetRasterBand(2).WriteArray( data_b2, 0 ,0 )
            dst_ds = None

            #print(stop)
            return 0

        def _geocode(self,s_date=None,e_date=None):
            
            doc = self.doc

            latfile = doc.offsetLatFile
            lonfile = doc.offsetLonFile
            losfile = doc.offsetLosFile

            # Figure out the bbox
            bbox_SNWE = self._get_bbox()

            print("bbox inferred from latlon files: ")
            bbox_SNWE = '"' + bbox_SNWE + '"'
            print(bbox_SNWE)

            if self.stack in ["tops","stripmap"]:
                
                # Define grid, Rutford Lat ~= 78.5, r = 1270 km
                # R/r = 6371 / 1270 = 5
                if doc.runid in [20190921, 20190904, 20190901]:
                    # Define grid to be 110m x 110m
                    # Lat uses R, Lon uses r
                    # 120/((2*math.pi*6371000)/360)
                    # Need five significant number to determine a grid point for slc pixel resolution
                    lat_step = 0.001
                    # Lon uses r
                    lon_step = 0.005

                # S1
                elif doc.runid == 20180703:
                    lat_step = 0.005
                    lon_step = 0.02

                # new S1 256 x 128
                elif doc.runid == 20200101:
                    lat_step = 0.005
                    lon_step = 0.02

                # new S1 480 x 128
                elif doc.runid == 20200102:
                    lat_step = 0.005
                    lon_step = 0.02

                else:
                    # Define grid, approximate 500 x 500 m (important) for Evans
                    #lon_step = 0.02
                    #lat_step = 0.005
                    raise Exception("Need lat lon grid size")  

            # Ridegecrest 200m x 200m
            else:
                lon_step = 0.002
                lat_step = 0.002

            # Geocode LOS file.
            print('Geocoding LOS file ...')
            cmd1 = 'geocodeGdal.py -l ' + latfile + ' -L ' + lonfile + ' -x ' + str(lon_step) + ' -y ' + str(lat_step) + ' -f ' + losfile + ' -b ' + bbox_SNWE + ' --suffix ' + self.version
            print(cmd1)

            filedir = os.path.dirname(losfile)
            losfilename, ending = os.path.basename(losfile).rsplit(".",1)
            gc_losfile = os.path.join(filedir,'gc_' + losfilename + '_' + self.version)+ "." + ending

            # Convert to observational vectors
            print('Convert LOS to observational vectors ...')
            cmd2 = 'los2enu.py -los {losfile}'.format(losfile = gc_losfile)
            print(cmd2)

            forceDoLos = False
            if self.exe or forceDoLos:
                # Delete the old files.
                os.system('rm ' + self.trackfolder + '/' + self.geometry + '/gc*' + str(doc.runid) + '*' + self.version + '*' )
                os.system('rm ' + self.trackfolder + '/' + self.geometry + '/temp*')

                # Geocoding
                os.system(cmd1)

                # Trim the los file to remove the extrapolation on the borders
                self._trim_borders(gc_losfile)

                # Conversion
                os.system(cmd2)
            else:
                print("Not executed")

            # Wait here until finish
            time.sleep(1)

            # generate the vectors in enu coordinates
            for count, offsetfield in enumerate(self.offsetfields):

                #if count>0:
                #    break

                date1str, date2str, date1, date2 = offsetfield

                # Skip the dates outside of the range.
                if s_date is not None and date1 < s_date.date():
                    continue

                if e_date is not None and date2 > e_date.date():
                    continue

                doc.interim = (date2-date1).days

                doc.offset_folder = os.path.join(self.trackfolder,self.offsetFolder,date1str+'_'+date2str)

                # The filtered offsetfields to be geocoded
                azOffset = os.path.join(doc.offset_folder, 'filtAzimuth_' + self.filt_suffix + '.off')
                rngOffset = os.path.join(doc.offset_folder, 'filtRange_' + self.filt_suffix + '.off')

                offset_outprefix = os.path.join(doc.offset_folder,doc.outprefix)

                # The covariance file to be geocoded
                covfile = offset_outprefix + '_run_' + str(doc.runid) + '_cov.bip'

                # The snr file to be geocoded
                snrfile = offset_outprefix + '_run_' + str(doc.runid) + '_snr.bip'

                if not os.path.exists(azOffset) or not os.path.exists(rngOffset):
                    print("The filtered offset field file doesn't exist")
                    print("skip ", date1str+'_'+date2str)
                    continue

                # The geocoded filenames
                gc_azOffset = os.path.join(doc.offset_folder, 'gc_filtAzimuth_' + self.filt_suffix + '.off')
                gc_rngOffset = os.path.join(doc.offset_folder, 'gc_filtRange_' + self.filt_suffix + '.off')


                #resamplingMethod = "near"
                resamplingMethod = "bilinear"

                cmd1 = 'geocodeGdal.py -l ' + latfile + ' -L ' + lonfile + ' -x ' + str(lon_step) + ' -y ' + str(lat_step) + ' -f ' + azOffset + ' -b ' + bbox_SNWE + ' -r ' + resamplingMethod

                cmd2 = 'geocodeGdal.py -l ' + latfile + ' -L ' + lonfile + ' -x ' + str(lon_step) + ' -y ' + str(lat_step) + ' -f ' + rngOffset + ' -b ' + bbox_SNWE + ' -r ' + resamplingMethod

                # geocode covariance file
                cmd3 = 'geocodeGdal.py -l ' + latfile + ' -L ' + lonfile + ' -x ' + str(lon_step) + ' -y ' + str(lat_step) + ' -f ' + covfile + ' -b ' + bbox_SNWE + ' -r ' + resamplingMethod + ' --suffix ' + self.version

                # geocode SNR file
                cmd4 = 'geocodeGdal.py -l ' + latfile + ' -L ' + lonfile + ' -x ' + str(lon_step) + ' -y ' + str(lat_step) + ' -f ' + snrfile + ' -b ' + bbox_SNWE + ' -r ' + resamplingMethod + ' --suffix ' + self.version


                # Parallel computing is problematic (seems to be fixed).

                forceDo = False
                #if not os.path.exists(gc_azOffset) or not os.path.exists(gc_rngOffset):
                #    print("The geocoded files do not exist!")
                #    print("work on :", offsetfield)
                #    forceDo = True

                # check if we should do it
                if self.exe or forceDo:
                    # Remove the old files
                    os.system('rm ' + doc.offset_folder + '/gc*' + str(doc.runid) + "*" + self.version + "*")
                    os.system('rm ' + doc.offset_folder + '/*temp*')

                    # Run the new commands
                    func =  globals()['run_denseoffset_bash']
                    p = Process(target=func, args=([cmd1,cmd2,cmd3,cmd4], (self.exe or forceDo)))
                    #p = Process(target=func, args=([cmd4], (self.exe or forceDo)))

                    self.jobs.append(p)
                    p.start()
 
                    self.doc_count = self.doc_count + 1

                    # nprocess controller
                    if self.doc_count >= self.nproc or count+1==len(self.offsetfields):
                        #for ip in range(self.nproc):
                        for ip in range(len(self.jobs)):
                            print(ip)
                            self.jobs[ip].join() #wait here

                        self.doc_count = 0
                        self.jobs = []

                else:
                    print(cmd1)
                    print(cmd2)
                    print(cmd3)
                    print(cmd4)
                    print("Not executed")

            return 0

        def geocode(self, s_date=None, e_date=None):
            self._geocode(s_date, e_date)

        ################## End of geocoding ##############################

        def _plot_offsetfield(self, offsetfield):

            doc = self.doc

            date1str, date2str, date1, date2 = offsetfield

            title = date1str + '_' + date2str

            doc.interim = (date2 - date1).days

            doc.offset_folder = os.path.join(self.trackfolder,self. offsetFolder, date1str + '_' + date2str)

            # Obtain filtered offset field
            azOffsetName = os.path.join(doc.offset_folder,'filtAzimuth_' + self.filt_suffix + '.off')
            rngOffsetName = os.path.join(doc.offset_folder, 'filtRange_' + self.filt_suffix + '.off')

            ds = gdal.Open(azOffsetName)
            filtered_azOffset = ds.GetRasterBand(1).ReadAsArray()

            ds = gdal.Open(rngOffsetName)
            filtered_rngOffset = ds.GetRasterBand(1).ReadAsArray()

            # Convert from the unit from pixel to meter per day
            filtered_azOffset = filtered_azOffset * doc.azPixelSize / doc.interim
            filtered_rngOffset = filtered_rngOffset *  doc.rngPixelSize / doc.interim

            # Get reference offset field
            refer_azOffset, refer_rngOffset, vmin_az, vmax_az, vmin_rng, vmax_rng = self._obtain_reference()

            # Get raw offset field
            offset_outprefix = os.path.join(doc.offset_folder, doc.outprefix)
            offsetfile = offset_outprefix + '_run_' + str(doc.runid) + '.bip'

            ds = gdal.Open(offsetfile)
            raw_azOffset = ds.GetRasterBand(1).ReadAsArray()
            raw_rngOffset = ds.GetRasterBand(2).ReadAsArray()

            # Convert from the unit from pixel to meter per day
            raw_azOffset = raw_azOffset * doc.azPixelSize / doc.interim
            raw_rngOffset = raw_rngOffset *  doc.rngPixelSize / doc.interim

            # Save the file as objects to disk and plot them later
            filtered_offset_folder = disp_temp_folder
            
            filtered_offset = {}
            filtered_offset['filtered_azOffset'] = filtered_azOffset
            filtered_offset['filtered_rngOffset'] = filtered_rngOffset
            filtered_offset['refer_azOffset'] = refer_azOffset
            filtered_offset['refer_rngOffset'] = refer_rngOffset
            filtered_offset['raw_azOffset'] = raw_azOffset
            filtered_offset['raw_rngOffset'] = raw_rngOffset

            figdir = os.path.join(self.workdir, 'figs', self.trackname,'radar')
            pathlib.Path(figdir).mkdir(exist_ok=True)

            filtered_offset['figdir'] = figdir
            
            # suffix
            filtered_offset['title'] = title + '_' + str(doc.runid) + '_' + self.version
            # Save to pickle file
            pkl_name = os.path.join(filtered_offset_folder, "track_" + self.trackname + "_" +title + ".pkl")
            with open(pkl_name, "wb") as f:
                pickle.dump(filtered_offset, f)

            # Plot it
            os.system("/net/kamb/ssd-tmp1/mzzhong/insarRoutines/stack_procs/display_az_rng.py " + pkl_name)

            # Remove the pickle file
            os.system("rm " + pkl_name)

            return 0 

        def plot_offsetfield(self, option="rdr"):

            if not option in ["rdr","geo"]:
                print("Wrong option provided: ", option)
                raise Exception()

            doc = self.doc

            for count, offsetfield in enumerate(self.offsetfields):

                print("Should work on: ", offsetfield)

                # check if we should do it
                if self.exe:

                    print("Plot it")
                    
                    func =  self._plot_offsetfield
                    p = Process(target=func, args=(offsetfield, ))
                    self.jobs.append(p)
                    p.start()
 
                    self.doc_count = self.doc_count + 1

                    # nprocess controller
                    if self.doc_count >= self.nproc or count+1==len(self.offsetfields):

                        for ip in range(len(self.jobs)):
                            print(ip)
                            self.jobs[ip].join() #wait here

                        self.doc_count = 0
                        self.jobs = []

            return 0

        ############### End of plotting offset field ##########################

        ############### Methods below are for data extration ###################

        # Create a stack for easy data access
        def createOffsetFieldStack(self):

            # Control redo the offset field stack
            self.offsetFieldStack_pkl = os.path.join(self.runfolder, self.offsetFolder, "offsetFieldStack_" + str(self.doc.runid) + "_" + self.version + ".pkl")

            redo_offsetFieldStack = 1
            
            if os.path.exists(self.offsetFieldStack_pkl) and redo_offsetFieldStack==0:
                print("offset field stack exists")
                with open(self.offsetFieldStack_pkl,"rb") as f:
                    self.offsetFieldStack = pickcle.load(f)
                    offsetFieldStack = self.offsetFieldStack

            else:

                print("Generating offset field stack")
                offsetFieldStack = {}

                # Load in the los file.
                geo_losFile = os.path.join(self.trackfolder,self.geometry, 'gc_los_offset_' + str(self.doc.runid) + '_' + self.version + '.rdr')
                print("geo_losFile: ",geo_losFile)
                ds = gdal.Open(geo_losFile)
                losField = ds.GetRasterBand(1).ReadAsArray()

                # Put los into the dictionary
                offsetFieldStack["los"] = losField

                # Put geometry into the dictionary
                geo_losVrtFile = geo_losFile + '.vrt'
                dataset = gdal.Open(geo_losVrtFile)
                geoTransform = dataset.GetGeoTransform()

                lon0 = geoTransform[0]
                lon_interval = geoTransform[1]

                lat0 = geoTransform[3]
                lat_interval = geoTransform[5]

                xsize = dataset.RasterXSize
                ysize = dataset.RasterYSize

                offsetFieldStack["coord"] = {}
                offsetFieldStack["coord"]["lon0"] = lon0
                offsetFieldStack["coord"]["lon_interval"] = lon_interval
                offsetFieldStack["coord"]["xsize"] = xsize
                offsetFieldStack["coord"]["lat0"] = lat0
                offsetFieldStack["coord"]["lat_interval"] = lat_interval
                offsetFieldStack["coord"]["ysize"] = ysize

                # Put offsetfield into the dictionary
                for offsetfield in self.offsetfields:
                    print(offsetfield)
    
                    date1str = offsetfield[0]
                    date2str = offsetfield[1]
    
                    date1 = offsetfield[2]
                    date2 = offsetfield[3]
    
                    title = date1str+'_'+date2str
    
                    offset_folder = os.path.join(self.trackfolder,self.offsetFolder,date1str+'_'+date2str)
        
                    geo_rngOffsetFile = os.path.join(offset_folder, 'gc_filtRange_' + self.filt_suffix  + '.off')
                    geo_azOffsetFile = os.path.join(offset_folder, 'gc_filtAzimuth_' + self.filt_suffix  + '.off')

                    if not os.path.exists(geo_rngOffsetFile) or not os.path.exists(geo_azOffsetFile):
                        print("File doesn't exist, skip: ", title)
                        print(geo_rngOffsetFile)
                        print(geo_azOffsetFile)
                        raise Exception()
                    
                    ds = gdal.Open(geo_rngOffsetFile)
                    rngOffsetField = ds.GetRasterBand(1).ReadAsArray()
                    
                    ds = gdal.Open(geo_azOffsetFile)
                    azOffsetField = ds.GetRasterBand(1).ReadAsArray()

                    offsetFieldStack[(title, "rng")] = rngOffsetField
                    offsetFieldStack[(title, "azi")] = azOffsetField

                print("Save offset field stack")
                with open(self.offsetFieldStack_pkl,"wb") as f:
                    pickle.dump(offsetFieldStack, f)

        def gc_lon_lat_axis(self):

            doc = self.doc

            gc_losfile = os.path.join(self.trackfolder,self.geometry, 'gc_los_offset_' + str(doc.runid) + '.rdr')

            gc_losvrtfile = gc_losfile + '.vrt'
            dataset = gdal.Open(gc_losvrtfile)
            geoTransform = dataset.GetGeoTransform()

            lon0 = geoTransform[0]
            lon_interval = geoTransform[1]

            lat0 = geoTransform[3]
            lat_interval = geoTransform[5]

            xsize = dataset.RasterXSize
            ysize = dataset.RasterYSize

            lon_list = np.linspace(lon0, lon0 + lon_interval*(xsize-1), xsize)
            lat_list = np.linspace(lat0, lat0 + lat_interval*(ysize-1), ysize)

            lon_list = np.round(lon_list * 1000)/1000
            lat_list = np.round(lat_list * 1000)/1000

            doc.lon_list = lon_list
            doc.lat_list = lat_list

        def int5d_to_float(self,x):
            return [num/(10**5) for num in x]

        # Find the index of the point in the geocoded offset field
        def point_index(self, point):

            doc = self.doc

            if not (hasattr(doc,'lon_list') and hasattr(doc,'lat_list')):
                self.gc_lon_lat_axis()

            lon_list = doc.lon_list
            lat_list = doc.lat_list

            lon, lat = point

            # Here, we need to convert lon, lat to float
            lon, lat = self.int5d_to_float([lon, lat])

            if len(np.where(lon_list == lon)[0])==1:
                ind_x = np.where(lon_list == lon)[0][0]
            else:
                ind_x = None

            if len(np.where(lat_list == lat)[0])==1:
                ind_y = np.where(lat_list == lat)[0][0]
            else:
                ind_y = None

            return (ind_x, ind_y)

        # 2019.12.05
        def point_set_index(self, point_set):

            # Obtain coordinates
            doc = self.doc

            if not (hasattr(doc,'lon_list') and hasattr(doc,'lat_list')):
                self.gc_lon_lat_axis()

            lon_list = doc.lon_list
            lat_list = doc.lat_list

            ind_set = {}
            # Lon and Lat are independent
            for point in point_set:

                lon, lat = point

                # Here, we need to convert lon, lat to float
                lon, lat = self.int5d_to_float([lon, lat])

                if len(np.where(lon_list == lon)[0])==1:
                    ind_x = np.where(lon_list == lon)[0][0]
                else:
                    ind_x = None
    
                if len(np.where(lat_list == lat)[0])==1:
                    ind_y = np.where(lat_list == lat)[0][0]
                else:
                    ind_y = None

                ind_set[point] = (ind_x, ind_y)

            return ind_set

        # Updated version of point_set_index (2019.12.05)
        def point_set_index_v2(self, point_set, coord=None):

            # Obtain coordinates
            doc = self.doc

            if coord is None:
                
                raise Exception("Coordinates of the stack is required")
                gc_losfile = os.path.join(self.trackfolder,self.geometry, 'gc_los_offset_' + str(doc.runid) + '.rdr')
   
                gc_losvrtfile = gc_losfile + '.vrt'

                dataset = gdal.Open(gc_losvrtfile)
                geoTransform = dataset.GetGeoTransform()
    
                lon0 = geoTransform[0]
                lon_interval = geoTransform[1]
    
                lat0 = geoTransform[3]
                lat_interval = geoTransform[5]
    
                xsize = dataset.RasterXSize
                ysize = dataset.RasterYSize

            else:
                lon0 = coord["lon0"]
                lon_interval = coord["lon_interval"]
    
                lat0 = coord["lat0"]
                lat_interval = coord["lat_interval"]
    
                xsize = coord["xsize"]
                ysize = coord["ysize"]

            print("Hey: ", lon0, lon_interval, lat0, lat_interval)

            # Lon and Lat are independent
            # Get array
            lon_arr = []
            lat_arr = []
            for point in point_set:
                lon, lat = point
                lon_arr.append(lon)
                lat_arr.append(lat)

            # Convert to float
            lon_arr = self.int5d_to_float(lon_arr)
            lat_arr = self.int5d_to_float(lat_arr)

            # Convert to numpy array
            lon_arr = np.asarray(lon_arr)
            lat_arr = np.asarray(lat_arr)

            # Calculate the index
            lon_arr_ind = np.round((lon_arr - lon0)/lon_interval)
            #lon_arr_ind = lon_arr_ind.astype("int")

            lat_arr_ind = np.round((lat_arr - lat0)/lat_interval)
            #lat_arr_ind = lat_arr_ind.astype("int")

            # Remove the invalid
            #print(lon_arr_ind)
            lon_arr_ind[lon_arr_ind<0] = np.nan
            lon_arr_ind[lon_arr_ind>=xsize] = np.nan
            #print(lon_arr_ind)

            #print(lat_arr_ind)
            lat_arr_ind[lat_arr_ind<0] = np.nan
            lat_arr_ind[lat_arr_ind>=ysize] = np.nan
            #print(lat_arr_ind)

            # Given to the array
            ind_set = {}
            for i, point in enumerate(point_set):
                if not np.isnan(lon_arr_ind[i]) and not np.isnan(lat_arr_ind[i]):
                    ind_set[point] = (int(lon_arr_ind[i]), int(lat_arr_ind[i]))
                else:
                    ind_set[point] = (None,None)

            return ind_set

        # Important: the beginning of data and offset field logistics
        def extract_offset_set_series(self, point_set=None, test_point=None, dates=None, offsetFieldStack=None, sate=None, track_num=None):

            # Strategy:
                # Open each offsetfield once, find offsets for all points.
                # Offsetfield outer loop.
                # Points inner loop.

            # Obtain the index of each point in this offsetFieldStack
            # Old        
            # ind_set = {}
            # for point in point_set:
            #     ind_set[point] = self.point_index(point)
            # New
            #ind_set = self.point_set_index(point_set)
            ind_set = self.point_set_index_v2(point_set,  offsetFieldStack["coord"])

            # Initialization.
            pairs_set = {}
            offsets_set = {}
            for point in point_set:
                pairs_set[point] = []
                offsets_set[point] = []

            # Obtain losField
            if offsetFieldStack is None:
                # Look at the offsetfields
                # Load in the los file.
                geo_losFile = os.path.join(self.trackfolder,self.geometry, 'gc_los_offset_' + str(self.doc.runid)+'.rdr')            
                #print(geo_losFile)
                ds = gdal.Open(geo_losFile)
                losField = ds.GetRasterBand(1).ReadAsArray()

            else:
                losField = offsetFieldStack['los']

            # Loop through the offset fields
            for offsetfield in self.offsetfields:

                date1str = offsetfield[0]
                date2str = offsetfield[1]

                date1 = offsetfield[2]
                date2 = offsetfield[3]

                # Only use dates in allowed range.
                if not ((date1 in dates) and (date2 in dates)):
                    continue

                title = date1str+'_'+date2str

                # Read in the original data
                if offsetFieldStack is None:

                    offset_folder = os.path.join(self.trackfolder,self.offsetFolder,date1str+'_'+date2str)
    
                    geo_rngOffsetFile = os.path.join(offset_folder,'gc_filtRange_' + self.filt_suffix  + '.off')
                    geo_azOffsetFile = os.path.join(offset_folder,'gc_filtAzimuth_' + self.filt_suffix  + '.off')
    
                    # Make sure the offsetfield exist.
                    if not os.path.exists(geo_rngOffsetFile) or not os.path.exists(geo_azOffsetFile):
                        continue
                    
                    ds = gdal.Open(geo_rngOffsetFile)
                    rngOffsetField = ds.GetRasterBand(1).ReadAsArray()
                    
                    ds = gdal.Open(geo_azOffsetFile)
                    azOffsetField = ds.GetRasterBand(1).ReadAsArray()

                # Read from the stack 2019.12.21
                else:
                    try:
                        rngOffsetField = offsetFieldStack[(title,"rng")]
                        azOffsetField = offsetFieldStack[(title,"azi")]
                    except:
                        print(offsetfield, track_num)
                        raise Exception("geocoded offset field missing!")

                ########################################################

                for point in point_set:

                    # Get the index
                    ind_x = ind_set[point][0]
                    ind_y = ind_set[point][1]

                    # Available in this track offsetfield, based on grid_point_set
                    
                    if ind_x is not None and ind_y is not None: 

                        # 2020.03.20
                        # However, it is possible and different resolution offset field dataset are not consistent
                        if ind_y < rngOffsetField.shape[0] and ind_x < rngOffsetField.shape[1]:
                            rngOffset = rngOffsetField[ind_y, ind_x]
                            azOffset = azOffsetField[ind_y, ind_x]
                            los = losField[ind_y, ind_x]

                            if point == test_point and track_num==10 and title=="20131115_20131123":
                                print("point: ", point)
                                print("rngOffset, azOffset: ", rngOffset, azOffset)
                                print("ind_y, ind_x: ", ind_y, ind_x)

                            # The value is valid.
                            # los > 5 is able to remove border/filled-in values
                            if (not np.isnan(rngOffset) and not np.isnan(azOffset) and los > 5):
    
                                # Convert offset from pixel to meter.
                                rngOffset = rngOffset * self.doc.rngPixelSize
                                azOffset = azOffset * self.doc.azPixelSize
    
                                # Comply to the los vector definition.
                                # Pointed from ground to satellite.
                                rngOffset = -rngOffset
    
                                # Save them
                                pairs_set[point].append([date1,date2])
                                offsets_set[point].append([rngOffset, azOffset])
                        else:
                            print("Hey, index is out of range: ", ind_y, ind_x, rngOffsetField.shape, point, sate, track_num) 
                            raise Exception("Index out of range")

                    #if point == test_point and track_num==10 and title=="20131115_20131123":
                    #    print("rngOffset, azOffset: ", rngOffset, azOffset)
                        #print(stop)

                # End of this offsetfield.

            return (pairs_set, offsets_set)

## End of dense_offset object
