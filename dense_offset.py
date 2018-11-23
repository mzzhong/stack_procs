#!/usr/bin/env python3
import os
import numpy as np
import glob

import datetime

from datetime import date

from multiprocessing import Process

import subprocess

import pickle

from scipy import ndimage 
from scipy import signal

import matplotlib.pyplot as plt
from matplotlib import cm

import isce
import isceobj

import gdal

from medianfilter2d import medianfilter2d as mdf

import cv2
 
def run_denseoffset_bash(cmd,exe=False):
    
    if exe:
        subprocess.call(cmd, shell=True)
    else:
        print(cmd)

class dense_offset_config():
    def __init__(self,ww=128, wh=128, sw=20, sh=20, kw=64, kh=64, nwac=5, nwdc=5):
        
        self.outprefix = 'cuampcor'
        self.runid = 1
        self.ww = ww
        self.wh = wh
        self.sw = sw
        self.sh = sh
        self.kw = kw
        self.kh = kh

        self.mm= max(max(self.sh,self.sw)*2, 50)
        self.gpuid = 0
        self.gross = 0
        self.deramp = 0

        self.nwac = nwac
        self.nwdc = nwdc

class dense_offset():

        def __init__(self, stack=None, workdir=None, nproc=None, exe=False):

            self.stack = stack
            self.workdir = workdir
            self.nproc = nproc
            self.exe = exe

            self.doc_count = 0
            self.jobs = []

            self.offsetFolder = 'cuDenseOffsets'

            if stack == 'stripmap':
                self.slcs_folder = 'coregSLC/Coarse'
                self.master_folder = 'raw'
                self.slave_folder = self.slcs_folder

                self.master_suffix = '.raw.slc'
                self.slave_suffix = '.slc'

                self.maxday = 8 

                # create the simplest doc object
                self.doc = dense_offset_config(ww=256, wh=256, sw=20, sh=20, kw=128, kh=128)

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
                self.doc = dense_offset_config(ww=256, wh=128, sw=20, sh=20, kw=128, kh=64)

                self.geometry = 'merged/geom_master'
                self.latname = 'lat.rdr.full'
                self.lonname = 'lon.rdr.full'
                self.losname = 'los.rdr.full'
                self.hgtname = 'hgt.rdr.full'

                self.geometry_suffix = '.rdr.full'

                # used to generate pixel size
                self.burst_xml = 'master/IW1.xml'


        def get_offset_geometry(self):

            doc = self.doc

            # Point to the geometry files.
            latfile = os.path.join(self.trackfolder, self.geometry, 'lat' + self.geometry_suffix)
            lonfile = os.path.join(self.trackfolder, self.geometry, 'lon' + self.geometry_suffix)
            losfile = os.path.join(self.trackfolder, self.geometry, 'los' + self.geometry_suffix)
            hgtfile = os.path.join(self.trackfolder, self.geometry, 'hgt' + self.geometry_suffix)

            doc.latfile = latfile 
            doc.lonfile = lonfile  
            doc.losfile = losfile  
            doc.hgtfile = hgtfile  

            # Offsetfield size.
            numWinDown = (doc.azsize - doc.mm*2 - doc.sh*2 - doc.wh) // doc.kh
            numWinAcross = (doc.rngsize - doc.mm*2 - doc.sw*2 - doc.ww) // doc.kw

            # Update doc.
            doc.numWinDown = numWinDown
            doc.numWinAcross = numWinAcross

            print(doc.ww, doc.wh, doc.kw, doc.kh)
            print(doc.numWinDown,doc.numWinAcross)

            # Check if the files exist.
            filelist = ['lat','lon','hgt','los']
            exist = True
            for geofile in filelist:
                offset_geofile = os.path.join(self.trackfolder, self.geometry, geofile + '_offset_' + str(doc.runid) + '.rdr.xml')
                print(offset_geofile)

                if not os.path.exists(offset_geofile):
                    exist = False

            # Allocate new geometry files            
            lon = np.zeros(shape=(numWinDown,numWinAcross),dtype=np.float64)
            lat = np.zeros(shape=(numWinDown,numWinAcross),dtype=np.float64) 
            hgt = np.zeros(shape=(numWinDown,numWinAcross),dtype=np.float64)

            inc = np.zeros(shape=(numWinDown,numWinAcross),dtype=np.float32)
            azi = np.zeros(shape=(numWinDown,numWinAcross),dtype=np.float32)

            centerOffsetHgt = doc.sh + doc.kh//2-1
            centerOffsetWidth = doc.sw + doc.kw//2-1

            print(numWinDown)
            print(numWinAcross)
            print(doc.kw)
            print(doc.kh)

            # Subset the original geometry file.
            for iwin in range(numWinDown):
                print(iwin)
                
                down = doc.mm + doc.kh * iwin  + centerOffsetHgt
                
                off = down * doc.rngsize
                off2 = down * doc.rngsize * 2 # BIL
    
                start = doc.mm + centerOffsetWidth
                end = doc.mm + doc.kw * numWinAcross
    
                latline = np.memmap(filename=latfile,dtype='float64',offset=8*off,shape=(doc.rngsize))
                lonline = np.memmap(filename=lonfile,dtype='float64',offset=8*off,shape=(doc.rngsize))
                hgtline = np.memmap(filename=hgtfile,dtype='float64',offset=8*off,shape=(doc.rngsize))
                 
                losline = np.memmap(filename=losfile,dtype='float32',offset=4*off2,shape=(doc.rngsize*2))

                lat[iwin,:] = latline[start:end:doc.kw]
                lon[iwin,:] = lonline[start:end:doc.kw]
                hgt[iwin,:] = hgtline[start:end:doc.kw]
               
                incline = losline[0:doc.rngsize]
                aziline = losline[doc.rngsize:doc.rngsize*2] 
                inc[iwin,:] = incline[start:end:doc.kw]
                azi[iwin,:] = aziline[start:end:doc.kw]

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

                offset_geofile = os.path.join(self.trackfolder, self.geometry, geofile + '_offset_' + str(doc.runid) + '.rdr')
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
 
            if self.stack == 'tops':
                pm = PM()
                pm.configure()
                obj = pm.loadProduct(os.path.join(self.trackfolder,self.burst_xml))
                burst = obj.bursts[0]
                Vs = np.linalg.norm(burst.orbit.interpolateOrbit(burst.sensingMid, method='hermite').getVelocity())
    
                print ('Velocity: ', Vs)
                print ('az Time interval: ', burst.azimuthTimeInterval)
                doc.azPixelSize=Vs*burst.azimuthTimeInterval
                print ('azPixelSize: ', doc.azPixelSize)
                doc.rngPixelSize=burst.rangePixelSize
                print ('rngPixelSize: ', doc.rngPixelSize)

            elif self.stack == 'stripmap':

                doc.azPixelSize = 2.286
                doc.rngPixelSize = 0.930
                print ('azPixelSize: ', doc.azPixelSize)
                print ('rngPixelSize: ', doc.rngPixelSize)

            return 0

        def reference(self):
            
            from grossOffsets import grossOffsets

            doc = self.doc
            
            objGrossOff = grossOffsets(figpath = os.path.join(self.trackfolder, self.offsetFolder))

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

            # setting the pixel size
            objGrossOff.setPixelSize(azPixelSize = doc.azPixelSize, rngPixelSize = doc.rngPixelSize)
            
            # set the days
            objGrossOff.setbTemp(1)

            # run to get reference offsetfield
            doc.grossDown, doc.grossAcross = objGrossOff.runGrossOffsets()

            #print(doc.grossDown)
            #print(np.nanmax(doc.grossDown), np.nanmin(doc.grossDown))
            #print(stop)

            return 0

        def show_doc(self):

            doc = self.doc
            
            print(doc.azsize)
            print(doc.rngsize)


        def initiate(self, trackname, runid):

            ## Initiate the basics
            
            self.trackname = trackname
            self.trackfolder = self.workdir + '/' + self.trackname
            
            # create the offset folder     
            if not os.path.exists(self.trackfolder + '/' + self.offsetFolder):
                os.mkdir(self.trackfolder + '/' + self.offsetFolder)

            # For stripmap, create the symbolic link to master
            if self.stack == 'stripmap':
                f = open(os.path.join(self.trackfolder,'run_files','run_2_master'))
                params = f.readlines()
                masterdate = params[0].split('_')[-1][:-1]
                print(masterdate)
            
                masterfolder = os.path.join(self.trackfolder, self.slcs_folder, masterdate)
                if not os.path.exists(masterfolder):
                    os.mkdir(masterfolder)
                    linked_prefix = os.path.join(self.trackfolder, self.master_folder, masterdate, masterdate)
                    linking_prefix = os.path.join(self.trackfolder, self.slcs_folder, masterdate, masterdate)
 
                    cmd1 = 'ln -s ' + linked_prefix + '.raw.slc ' + linking_prefix + '.slc'
                    cmd2 = 'ln -s ' + linked_prefix + '.raw.slc.xml ' + linking_prefix + '.slc.xml'
                    cmd3 = 'ln -s ' + linked_prefix + '.raw.slc.vrt ' + linking_prefix + '.slc.vrt'

                    print(cmd1)
                    print(cmd2)
                    print(cmd3)

                    os.system(cmd1)
                    os.system(cmd2)
                    os.system(cmd3)

            # Collect dates
            slcs = glob.glob(self.trackfolder + '/' + self.slcs_folder + '/2*')
            intdates = [int(slc.split('/')[-1]) for slc in slcs]
            
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

            #print(self.offsetfields)

            
            ## offsetfield object initiation.
            # Load the existed result, if it exist. 
            # Point to the offset pickle file.
            self.doc_pkl = os.path.join(self.trackfolder,self.offsetFolder,str(runid) + '.pkl')
            
            redo = 0
            if os.path.exists(self.doc_pkl) and redo==0:
                with open(self.doc_pkl,'rb') as f:
                    self.doc = pickle.load(f)
                    doc = self.doc
 
            else:

                ## Update offsetfield object.
                doc = self.doc
                doc.runid = runid

                # Arbitrarily choose a SLC Xml file
                xmls = glob.glob(self.trackfolder + '/' + self.slcs_folder + '/2*/*'+self.slave_suffix+'.xml')
                print(xmls)
 
                # Get the SLC image size
                self.get_slc_size(xmls[0])
            
                # Get the offsetfield pixel size
                
                self.get_pixel_size()

                # Obtain the offset geometry
                self.get_offset_geometry()
             
                # Job-level initiation 
                # doc.winsize = ...

                # obtain reference velocity
                self.reference()

                # runtime gpu params
                doc.nwac = 10
                doc.nwdc = 10

                # save it
                with open(self.doc_pkl ,'wb') as f:
                    pickle.dump(doc,f)

            return 0

        def run_offset_track(self):

            doc = self.doc
            
            # Esitmate dense offsets for one track.
            for offsetfield in self.offsetfields:
    
                date1str = offsetfield[0]
                date2str = offsetfield[1]

                offset_folder = os.path.join(self.trackfolder,self.offsetFolder,date1str+'_'+date2str)
                
                if not os.path.exists(offset_folder):
                    os.mkdir(offset_folder)
                    print(offset_folder)
                    #print(stop)

                # Generate runfile for dense_offset estimation.
                master = os.path.join(self.trackfolder,self.slcs_folder,date1str, date1str + self.slave_suffix)
                slave = os.path.join(self.trackfolder,self.slcs_folder,date2str, date2str + self.slave_suffix )

                # set parameters
                # the last gpus
                doc.gpuid = self.doc_count % self.nproc + (8-self.nproc)
                print(doc.gpuid)
                #print(offset_folder)
                #print(doc.outprefix)
                
                #doc.outprefix = offset_folder + doc.outprefix

                offset_outprefix = os.path.join(offset_folder,doc.outprefix)

                #print(offset_outprefix)

                run_file = offset_folder + '/gpurun' + '_' + str(doc.runid) + '.sh'
                f = open(run_file,'w')
                f.write('#!/bin/bash\n')
                f.write('date1=' + date1str + '\n')
                f.write('date2=' + date2str+ '\n')
                f.write('master=' + master+ '\n')
                f.write('slave=' + slave+ '\n')

                f.write('fixImageXml.py -f -i $master'+ '\n')
                f.write('fixImageXml.py -f -i $slave'+ '\n')

                f.write('ww=' + str(doc.ww)+ '\n')
                f.write('wh=' + str(doc.wh)+ '\n')
                f.write('sw=' + str(doc.sw)+ '\n')
                f.write('sh=' + str(doc.sh)+ '\n')
                f.write('kw=' + str(doc.kw)+ '\n')
                f.write('kh=' + str(doc.kh)+ '\n')
                f.write('mm=' + str(doc.mm)+ '\n')
                f.write('gross=' + str(doc.gross)+ '\n')
                f.write('gpuid=' + str(doc.gpuid)+ '\n')
                f.write('deramp=' + str(doc.deramp)+ '\n')

                f.write('nwac=' + str(doc.nwac)+ '\n')
                f.write('nwdc=' + str(doc.nwdc)+ '\n')

                f.write('outprefix='+str(offset_outprefix)+ '\n')

                f.write('runid='+str(doc.runid)+'\n')  
                f.write('outsuffix=' + '_run_$runid'+'\n')

                #f.write('rm ' + '$outprefix$outsuffix' + '*\n')
                f.write('cuDenseOffsets.py --master $master --slave $slave --ww $ww --wh $wh --sw $sw --sh $sh --mm $mm --kw $kw --kh $kh --gross $gross --outprefix $outprefix --outsuffix $outsuffix --deramp $deramp --gpuid $gpuid --nwac $nwac --nwdc $nwdc'+ '\n')

                f.close()

                func = globals()['run_denseoffset_bash']
                p = Process(target=func, args=('bash '+run_file,self.exe))
                self.jobs.append(p)
                p.start()
                self.doc_count = self.doc_count + 1

                # nprocess controller

                if self.doc_count == self.nproc:
                    for ip in range(self.nproc):
                        self.jobs[ip].join() #wait here

                    self.doc_count = 0
                    self.jobs = []

            return 0

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

            ind = np.where(np.isnan(data))
            num_nan = len(ind[0])

            allvalues = np.zeros(shape=(num_nan,hole_size*hole_size))
            holevalue = np.zeros(shape=(num_nan,))

            count = 0
            for i in range(hole_size):
                for j in range(hole_size):
                    #print(i,j)
                    ind_i = ind[0] + i - hole_size//2
                    ind_i[ind_i<0]=0
                    ind_i[ind_i>=doc.numWinDown]=doc.numWinDown-1

                    ind_j = ind[1] + j - hole_size//2
                    ind_j[ind_j<0]=0
                    ind_j[ind_j>=doc.numWinAcross]=doc.numWinAcross-1

                    allvalues[:,count] = data[(ind_i,ind_j)] 
                    count = count + 1

            count_nans = np.isnan(allvalues).sum(axis=1)
            ratio = 1/2
            invalid_points = count_nans > np.rint(hole_size * hole_size * ratio)

            # fill in
            holevalues = np.nanmedian(allvalues, axis=1)
    
            # invalid points
            holevalues[invalid_points] = np.nan

            # median filter
            data = np.copy(data)
            data[ind] = holevalues

            kernel_size = (5,5)
            #data = ndimage.median_filter(input=data,size=kernel_size,mode='nearest')
            #data = signal.medfilt(data,kernel_size=kernel_size)
            data = mdf(data,kernel_size=kernel_size)
            #data = cv2.medianBlur(data,kernel_size[0])

            return data

        def _get_refer_mask(self,az,refer_az, rng, refer_rng, thrs=None):

            doc = self.doc

            # difference meter per day
            dif_az_day = az - refer_az
            dif_rng_day = rng - refer_rng
 
            # important
            # threshold 1 meter per day

            if thrs is None:
                thrs = 1

            mask = np.logical_or(np.abs(dif_az_day)>thrs, np.abs(dif_rng_day)>thrs)

            return mask

        def _get_misCoreg(self,az,refer_az, rng, refer_rng, mask, order=1):

            doc = self.doc

            # Unmoved places, smaller than thr m/d.
            thr = 0.05
            
            ind_az = np.abs(refer_az)<thr
            ind_rng = np.abs(refer_rng)<thr

            # Consider the points smaller than the threshold and validly obtained.
            ind_stationary = np.logical_and(np.logical_and(ind_az,ind_rng),np.invert(mask))

            ind_d = np.arange(0,doc.numWinDown)
            ind_a = np.arange(0,doc.numWinAcross)

            aa,dd = np.meshgrid(ind_a,ind_d)

            # Available points for calculating micoregistration
            ind_d_mis = dd[ind_stationary]
            ind_a_mis = aa[ind_stationary]

            az_value_mis = az[ind_stationary] - refer_az[ind_stationary]
            rng_value_mis = rng[ind_stationary] - refer_rng[ind_stationary]

            az_mat_mis = np.full(shape=(doc.numWinDown, doc.numWinAcross), fill_value=np.nan)
            az_mat_mis[ind_d_mis, ind_a_mis] = az_value_mis

            rng_mat_mis = np.full(shape=(doc.numWinDown, doc. numWinAcross), fill_value=np.nan)
            rng_mat_mis[ind_d_mis, ind_a_mis] =  rng_value_mis

            # 2D polyfit
            # Azimuth
            num = len(az_value_mis)
            col1 = np.full(shape=(num,), fill_value=1)
            G = np.stack([col1,ind_a_mis,ind_d_mis],axis=1)

            (c,a,b) = np.linalg.lstsq(G, az_value_mis, rcond=None)[0]
            az_mat_mis_pred = a * aa + b*dd + c

            # Range
            num = len(rng_value_mis)
            col1 = np.full(shape=(num,), fill_value=1)
            G = np.stack([col1,ind_a_mis,ind_d_mis],axis=1)

            (c,a,b) = np.linalg.lstsq(G, rng_value_mis, rcond=None)[0]
            rng_mat_mis_pred = a * aa + b*dd + c

            #fig = plt.figure(figsize=(10,10))
            
            #vmin = np.min(rng_mat_mis_pred)
            #vmax = np.max(rng_mat_mis_pred)
            
            #ax = fig.add_subplot(121)
            #f1 = ax.imshow(rng_mat_mis,cmap=plt.cm.jet,vmin=vmin,vmax=vmax)

            #ax = fig.add_subplot(122)
            #f1 = ax.imshow(rng_mat_mis_pred,cmap=plt.cm.jet,vmin=vmin,vmax=vmax)
            #fig.colorbar(f1)

            #fig.savefig('test2.png', format='png')
 
            #if order==0:
            #    az_mis = np.nanmedian(az[ind_stationary] - refer_az[ind_stationary])
            #    rng_mis = np.nanmedian(rng[ind_stationary] - refer_rng[ind_stationary])

            #print(az_mis,rng_mis)

            return (az_mat_mis_pred,rng_mat_mis_pred)
 
        def _run_offset_filter(self,data,data_snr,mask=None, label=None, refer=None):

            doc = self.doc
           
            med1_data = np.copy(data)

            # Step 0: mask out using snr
            med1_data[data_snr<4] = np.nan

            # Step 1: mask out where the reference is nan
            med1_data[np.isnan(refer)] = np.nan

            if mask is not None: 
                med1_data[mask] = np.nan

            # Step 2: normal median filter
            med2_data = np.copy(data)
            med2_kernel_size = (7,7)

            # median filters
            #med2_data = ndimage.median_filter(input=med2_data,size=med2_kernel_size,mode='nearest') 
            #med2_data = signal.medfilt(med1_data,kernel_size=med2_kernel_size)
            med2_data = mdf(med1_data,kernel_size=med2_kernel_size)
            #med2_data = cv2.medianBlur(med1_data,med2_kernel_size[0])

            # step 3: iteratively fill in small holes, fill in + median filter
            iteration = 3
            hole_size = 7
            med3_data = np.copy(med2_data)
            for i in range(iteration):
                med3_data = self.fill_in_holes(data=med3_data, hole_size = hole_size)
                #med3_data = self.fill(med3_data)

            ## step 4: set neighboring values to nans to be nan
            med4_kernel_size = 7
            med4_data = np.copy(med3_data)
            ind = np.where(np.isnan(med4_data))

            for i in range(med4_kernel_size):
                for j in range(med4_kernel_size):
                    ind_i = ind[0] + i - med4_kernel_size//2
                    ind_i[ind_i<0]=0
                    ind_i[ind_i>=doc.numWinDown]=doc.numWinDown-1

                    ind_j = ind[1] + j - med4_kernel_size//2
                    ind_j[ind_j<0]=0
                    ind_j[ind_j>=doc.numWinAcross]=doc.numWinAcross-1

                    med4_data[(ind_i,ind_j)] = np.nan

            # finalize, remove the boundary
            #med4_data[np.isnan(refer)] = np.nan
            bb = (iteration + 1) * (hole_size // 2)
            med4_data[:iteration,:] = np.nan
            med4_data[-iteration:,:] = np.nan
            med4_data[:,:iteration] = np.nan
            med4_data[:,-iteration:] = np.nan

            return med4_data

        def _display_az_rng(self, azOff, rngOff, azOff_re, rngOff_re, azOff_cmp, rngOff_cmp, title= None):

            doc = self.doc

            pad = 1

            # Ranges of values.
            # Plot azimuth offset.
            vmin = np.nanmin(azOff_re)
            vmax = np.nanmax(azOff_re)

            vmin10 = np.floor(vmin*10)/10-pad
            vmax10 = np.ceil(vmax*10)/10+pad

            #vmin10 = -0.05
            #vmax10 = 0.05

            #print(vmin10, vmax10)
            #print(stop)

            self.fig = plt.figure(1, figsize=(18,9))

            ax = self.fig.add_subplot(161)
            ax.set_title('predicted azimuth offset') 
            im = ax.imshow(azOff_re, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            self.fig.colorbar(im,fraction=0.07, orientation='horizontal',ticks=np.linspace(vmin10,vmax10,num=5),label='meter')
            #self.fig.colorbar(im,fraction=0.07, orientation='horizontal',label='meter')



            ax = self.fig.add_subplot(162)
            ax.set_title('raw azimuth offset') 
            im = ax.imshow(azOff_cmp, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            self.fig.colorbar(im,fraction=0.07, orientation='horizontal',ticks=np.linspace(vmin10,vmax10,num=5),label='meter')

            ax = self.fig.add_subplot(163)
            ax.set_title('filtered azimuth offset')
            im = ax.imshow(azOff, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            self.fig.colorbar(im,fraction=0.07, orientation='horizontal',ticks=np.linspace(vmin10,vmax10,num=5),label='meter')

            # Plot range offset.
            vmin = np.nanmin(rngOff_re)
            vmax = np.nanmax(rngOff_re)

            vmin10 = np.floor(vmin*10)/10-pad
            vmax10 = np.ceil(vmax*10)/10+pad

            ax = self.fig.add_subplot(164)
            ax.set_title('predicted range offset')
            im = ax.imshow(rngOff_re, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            self.fig.colorbar(im,fraction=0.07, orientation='horizontal',ticks=np.linspace(vmin10,vmax10,num=5),label='meter')

            ax = self.fig.add_subplot(165)
            ax.set_title('raw range offset') 
            im = ax.imshow(rngOff_cmp, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            self.fig.colorbar(im,fraction=0.07, orientation='horizontal',ticks=np.linspace(vmin10,vmax10,num=5),label='meter')

            ax = self.fig.add_subplot(166)
            ax.set_title('filtered range offset')          
            im = ax.imshow(rngOff, cmap=cm.jet, vmax=vmax10, vmin=vmin10)
            self.fig.colorbar(im,fraction=0.07, orientation='horizontal',ticks=np.linspace(vmin10,vmax10,num=5),label='meter')

            # Finalize.
            self.fig.suptitle(title)
            #plt.show()
            
            figdir = os.path.join(os.path.join(self.workdir,'figs', self.trackname,'radar'))
            if not os.path.exists(figdir):
                os.makedirs(figdir)

            self.fig.savefig(os.path.join(figdir,'offset_' + title + ".pdf"), format='pdf')
            self.fig.savefig(os.path.join(figdir,'offset_' + title + ".png"), format='png')

            return 0

        def _save_az_rng(self, azOffset_filtered, rngOffset_filtered):

            doc = self.doc

            redo = 1

            os.system('rm -r ' + doc.offset_folder + '/filt*' + str(doc.runid)+'*')

            azOffsetName = os.path.join(doc.offset_folder, 'filtAzimuth_' + str(doc.runid) +'.off')
            nbands = 1
            bandType = 6

            # Convert back to pixels in real interim days.
            azOffset_filtered_ori = azOffset_filtered / doc.azPixelSize * doc.interim
            rngOffset_filtered_ori = rngOffset_filtered / doc.rngPixelSize * doc.interim
 
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
           
            ################ 

            rngOffsetName = os.path.join(doc.offset_folder, 'filtRange_' + str(doc.runid) + '.off')
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


        def _offset_filter(self,offsetfile,snrfile,title):

            # Filter the offset fields.
            # The intermediate unit is meter/day.
            # The final unit is still in pixel.

            doc = self.doc

            ds = gdal.Open(offsetfile)
            azOffset = ds.GetRasterBand(1).ReadAsArray()
            rngOffset = ds.GetRasterBand(2).ReadAsArray()

            ds = gdal.Open(snrfile)
            data_snr = ds.GetRasterBand(1).ReadAsArray()

            # reference
            refer_azOffset = doc.grossDown
            refer_rngOffset = doc.grossAcross

            # convert reference from pixels to meter per day.
            refer_azOffset = refer_azOffset * doc.azPixelSize
            refer_rngOffset = refer_rngOffset * doc.rngPixelSize
 
            # save the reference (already saved as grossDown/Across)
            #doc.refer_azOffset = refer_azOffset
            #doc.refer_rngOffset = refer_rngOffset


            # Convert observation from pixels to meter per day.
            azOffset = azOffset * doc.azPixelSize / doc.interim
            rngOffset = rngOffset * doc.rngPixelSize / doc.interim

            # Generate first mask for stagnent regions, used for deriving mis-coreg.
            mask = self._get_refer_mask(azOffset, refer_azOffset, rngOffset, refer_rngOffset, thrs=1)

            # Get the mis-coregistration.
            az_mis, rng_mis = self._get_misCoreg(azOffset, refer_azOffset, rngOffset, refer_rngOffset, mask)
            #print(az_mis, rng_mis)

            # Mis-coregistration correction.
            azOffset = azOffset - az_mis
            rngOffset = rngOffset - rng_mis

            # Generate second mask, used for filtering.
            mask = self._get_refer_mask(azOffset, refer_azOffset, rngOffset, refer_rngOffset, thrs=1)

            # Filtering.
            azOffset_filtered = self._run_offset_filter(azOffset, data_snr, mask=mask, label='az',refer=refer_azOffset)
            rngOffset_filtered = self._run_offset_filter(rngOffset, data_snr, mask=mask, label='rng',refer=refer_rngOffset)

            # Cross nan validation.
            nan_ind = np.logical_or(np.isnan(azOffset_filtered), np.isnan(rngOffset_filtered))
            azOffset_filtered[nan_ind] = np.nan
            rngOffset_filtered[nan_ind] = np.nan

            # Display.
            self._display_az_rng(azOffset_filtered, rngOffset_filtered, refer_azOffset, refer_rngOffset, azOffset, rngOffset, title)

            # Convert the unit from meter per day to pixel.
            # Save the offset fields.
            self._save_az_rng(azOffset_filtered, rngOffset_filtered)

            
            return 0

        def offset_filter(self,s_date,e_date):

            doc = self.doc

            for offsetfield in self.offsetfields:
    
                date1str = offsetfield[0]
                date2str = offsetfield[1]

                date1 = offsetfield[2]
                date2 = offsetfield[3]

                # Skip the dates outside of the range.
                if date1 < s_date.date() or date2 > e_date.date():
                    continue

                doc.interim = (date2-date1).days

                doc.offset_folder = os.path.join(self.trackfolder,self.offsetFolder,date1str+'_'+date2str)

                offset_outprefix = os.path.join(doc.offset_folder,doc.outprefix)

                # Prepare the files.
                offsetfile = offset_outprefix + '_run_' + str(doc.runid) + '.bip'
                snrfile = offset_outprefix + '_run_' + str(doc.runid) + '_snr.bip'


                title = date1str+'_'+date2str

                # check if we should do it
                exist = False

                if exist == False:
                    print('Should work on', title)
                    if self.exe == True:
                        print('work on it!')
                        func = self._offset_filter
                        p = Process(target=func, args=(offsetfile,snrfile,title,))
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

        def _geocode(self,s_date=None,e_date=None):
            
            doc = self.doc

            latfile = doc.offsetLatFile
            lonfile = doc.offsetLonFile

            losfile = doc.offsetLosFile

            # Define grid, approximate 500 x 500 m (important)
            lon_step = 0.02
            lat_step = 0.005

            # Geocode LOS file.
            print('Geocoding LOS file ...')
            cmd1 = 'geocodeGdal.py -l ' + latfile + ' -L ' + lonfile + ' -x ' + str(lon_step) + ' -y ' + str(lat_step) + ' -f ' + losfile
            print(cmd1)

            filedir = os.path.dirname(losfile)
            losfilename = os.path.basename(losfile)
            gc_losfile = os.path.join(filedir,'gc_'+losfilename) 

            print('Convert LOS to observational vectors ...')
            cmd2 = 'los2enu.py -los {losfile}'.format(losfile = gc_losfile)
            print(cmd2)

            if self.exe:
                #os.system('rm ' + self.trackfolder + '/' + self.geometry + '/gc*')
                #os.system('rm ' + self.trackfolder + '/' + self.geometry + '/temp*')
                os.system(cmd1)
                os.system(cmd2)

            # generate the vectors in enu coordinates

            for offsetfield in self.offsetfields:
    
                date1str = offsetfield[0]
                date2str = offsetfield[1]

                date1 = offsetfield[2]
                date2 = offsetfield[3]

                # Skip the dates outside of the range.
                if date1 < s_date.date() or date2 > e_date.date():
                    continue


                doc.interim = (date2-date1).days

                doc.offset_folder = os.path.join(self.trackfolder,self.offsetFolder,date1str+'_'+date2str)

                # geocode offsetfields

                azOffset = os.path.join(doc.offset_folder, 'filtAzimuth_' + str(doc.runid) + '.off')
                rngOffset = os.path.join(doc.offset_folder, 'filtRange_' + str(doc.runid) + '.off')

                cmd1 = 'geocodeGdal.py -l ' + latfile + ' -L ' + lonfile + ' -x ' + str(lon_step) + ' -y ' + str(lat_step) + ' -f ' + azOffset
                cmd2 = 'geocodeGdal.py -l ' + latfile + ' -L ' + lonfile + ' -x ' + str(lon_step) + ' -y ' + str(lat_step) + ' -f ' + rngOffset

                # Parallel computing is problematic.

                ## check if we should do it
                #exist = False
                #if exist == False:
                #    print(cmd1)
                #    print(cmd2)

                #    if self.exe:
                #        os.system('rm ' + doc.offset_folder + '/gc*')
                #        os.system('rm ' + doc.offset_folder + '/*temp*')
                #    
                #        func =  globals()['run_denseoffset_bash']
                #        p = Process(target=func, args=(cmd1,self.exe))
                #        self.jobs.append(p)
                #        p.start()
 
                #        p = Process(target=func, args=(cmd2,self.exe))
                #        self.jobs.append(p)
                #        p.start()
                #        
                #        self.doc_count = self.doc_count + 2

                #        # nprocess controller
                #        if self.doc_count >= self.nproc:
                #            for ip in range(self.nproc):
                #                print(ip)
                #                self.jobs[ip].join() #wait here

                #            self.doc_count = 0
                #            self.jobs = []

                # Non-parallel computing.
                if self.exe:
                   #os.system('rm ' + doc.offset_folder + '/gc*')
                   #os.system('rm ' + doc.offset_folder + '/*temp*')
 
                    os.system(cmd1)
                    os.system(cmd2)
                else:
                    print(cmd1)
                    print(cmd2)


            return 0

        def geocode(self, s_date=None, e_date=None):
            self._geocode(s_date, e_date)

        def _plot_geocoded(self):

            from EvansPlot import EvansPlot
            from offsetField import offsetField

            doc = self.doc

            for offsetfield in self.offsetfields:

                date1str = offsetfield[0]
                date2str = offsetfield[1]

                date1 = offsetfield[2]
                date2 = offsetfield[3]

                title = date1str+'_'+date2str

                doc.interim = (date2-date1).days

                doc.offset_folder = os.path.join(self.trackfolder,self.offsetFolder,date1str+'_'+date2str)

                geo_azOffset = os.path.join(doc.offset_folder,'gc_filtAzimuth_' + str(doc.runid) + '.off')
                geo_rngOffset = os.path.join(doc.offset_folder,'gc_filtRange_' + str(doc.runid) + '.off')

                geo_los = os.path.join(self.trackfolder,self.geometry,'gc_los_offset_' + str(doc.runid) + '.rdr')


                # Set up the geocoded offsetfield
                myOff = offsetField()
                
                myOff.setOffsetField(geo_azOffset, geo_rngOffset)
                myOff.setLos(geo_los)
                
                myOff.setPixelSize(doc.azPixelSize, doc.rngPixelSize)
                myOff.setDays(doc.interim)

                # Plot it
                plot_Object = ['speed','azi','rng']
                plot_Object = ['speed']

                # Speed maps.
                if 'speed' in plot_Object:

                    # Choose one offsetfield to plot.
                    # Here, temporarily choose the first offsetfield.

                    if offsetfield == self.offsetfields[0]:

                        speedPlot = EvansPlot()
                        speedPlot.setup(name = 'Evans')
                        speedPlot.showGL()

                        figdir = os.path.join(os.path.join(self.workdir,'figs', self.trackname,'geocoded'))
                        if not os.path.exists(figdir):
                            os.makedirs(figdir)
                        
                        myOff.showSpeedMap(speedPlot,colorbar = True)
                    else:
                        myOff.showSpeedMap(speedPlot,colorbar = False)

                    figtitle = "Evans Ice Stream speed"
                    speedPlot.ax.set_title(figtitle,fontsize=15)

                    figname = os.path.join(figdir,'gc_speed_' + title)
 
                    speedPlot.fig.savefig(figname + ".png", format='png')

                    # undo
                    myOff.undoMap(speedPlot)


                # Azimuthal offsets.
                # Figure out the approximate value range from velocity model.
                # reference
                refer_azOffset = doc.grossDown
                refer_azOffset = refer_azOffset * doc.azPixelSize
                vmin = np.nanmin(refer_azOffset)
                vmax = np.nanmax(refer_azOffset)

                if 'azi' in plot_Object:

                    if offsetfield == self.offsetfields[0]:

                        aziPlot = EvansPlot()
                        aziPlot.setup(name = 'Evans')
                        aziPlot.showGL()

                        figdir = os.path.join(os.path.join(self.workdir,'figs', self.trackname,'geocoded'))
                        if not os.path.exists(figdir):
                            os.makedirs(figdir)
                        
                        myOff.showAzimuthMap(aziPlot,colorbar = True, value_range=[vmin,vmax])
                    else:
                        myOff.showAzimuthMap(aziPlot,colorbar = False, value_range=[vmin,vmax])

                    figtitle = "Evans Ice Stream azimuth offset map"
                    aziPlot.ax.set_title(figtitle,fontsize=15)

                    figname = os.path.join(figdir,'gc_azi_' + title)
 
                    aziPlot.fig.savefig(figname + ".png", format='png')

                    # undo
                    myOff.undoMap(aziPlot)


                # Range offsets.
                refer_rngOffset = doc.grossAcross
                refer_rngOffset = refer_rngOffset * doc.rngPixelSize
                vmin = np.nanmin(refer_rngOffset)
                vmax = np.nanmax(refer_rngOffset)

                if 'rng' in plot_Object:

                    if offsetfield == self.offsetfields[0]:

                        rngPlot = EvansPlot()
                        rngPlot.setup(name = 'Evans')
                        rngPlot.showGL()

                        figdir = os.path.join(os.path.join(self.workdir,'figs', self.trackname,'geocoded'))
                        if not os.path.exists(figdir):
                            os.makedirs(figdir)
                        
                        myOff.showRangeMap(rngPlot,colorbar = True, value_range=[vmin,vmax])
                    else:
                        myOff.showRangeMap(rngPlot,colorbar = False, value_range=[vmin,vmax])

                    figtitle = "Evans Ice Stream range offset map"
                    rngPlot.ax.set_title(figtitle,fontsize=15)

                    figname = os.path.join(figdir,'gc_rng_' + title)
 
                    rngPlot.fig.savefig(figname + ".png", format='png')

                    # undo
                    myOff.undoMap(rngPlot)
            return 0

        def _plot_geocoded_all(self,last=True):

            from EvansPlot import EvansPlot
            from offsetField import offsetField

            doc = self.doc

            # Pick the best offsetfield.
            minNanCount = None

            for offsetfield in self.offsetfields:

                self.doc_count = self.doc_count + 1
                print(self.doc_count)
    
                date1str = offsetfield[0]
                date2str = offsetfield[1]

                date1 = offsetfield[2]
                date2 = offsetfield[3]

                title = date1str+'_'+date2str

                doc.interim = (date2-date1).days

                doc.offset_folder = os.path.join(self.trackfolder,self.offsetFolder,date1str+'_'+date2str)

                geo_azOffsetFile = os.path.join(doc.offset_folder,'gc_filtAzimuth_' + str(doc.runid) + '.off')
                geo_rngOffsetFile = os.path.join(doc.offset_folder,'gc_filtRange_' + str(doc.runid) + '.off')

                geo_losFile = os.path.join(self.trackfolder,self.geometry,'gc_los_offset_' + str(doc.runid) + '.rdr')

                # read the data and count the number nan values
                ds = gdal.Open(geo_azOffsetFile)
                geo_azOffset = ds.GetRasterBand(1).ReadAsArray()
                
                #print(geo_azOffset)
                #print(geo_azOffset[19,39])
                #print(geo_azOffset.shape)
                #print(type(geo_azOffset))
                #print(geo_azOffset[0,0] ==0 )

                nanCount = np.isnan(geo_azOffset).sum()
                #print(geo_azOffsetFile)
                #print(nanCount)

                if minNanCount is None or  nanCount < minNanCount:
                    minNanCount = nanCount
                    best_azFile = geo_azOffsetFile
                    best_rngFile = geo_rngOffsetFile
                    interim_days = doc.interim

                # break


            # Plot the best offsetfield.
            print('the best file')
            print(best_azFile)
            print(best_rngFile)
            print(minNanCount)
            print(interim_days)

            # Set up the geocoded offsetfield.
            myOff = offsetField()
            
            myOff.setOffsetField(best_azFile, best_rngFile)
            myOff.setLos(geo_losFile)
            
            myOff.setPixelSize(doc.azPixelSize, doc.rngPixelSize)
            myOff.setDays(interim_days)
            print(interim_days)

            # Plot it.
            plot_Object = ['speed','azi','rng']
            plot_Object = ['speed']

            if 'speed' in plot_Object:
                try:
                    self.speedPlotAll
                except:
                    self.speedPlotAll = EvansPlot()
                    self.speedPlotAll.setup(name = 'Evans')
                    self.speedPlotAll.showGL()

                    self.figdir = os.path.join(os.path.join(self.workdir, 'figs'))
                    if not os.path.exists(self.figdir):
                        os.makedirs(self.figdir)
                    
                    myOff.showSpeedMap(self.speedPlotAll,colorbar = True)
                    figtitle = "Evans Ice Stream speed map (descending tracks)"
                    self.speedPlotAll.ax.set_title(figtitle,fontsize=15)

                else:
                    myOff.showSpeedMap(self.speedPlotAll,colorbar = False)

                if last:
                    figname = os.path.join(self.figdir,'gc_speed')
                    self.speedPlotAll.fig.savefig(figname + ".png", format='png')

            return 0
       
        def plot_geocoded(self,label=None,last=True):

            if (label is None) or (label == 'seperate'):
                self._plot_geocoded()
            elif label == 'all':
                self._plot_geocoded_all(last=last)

        def postprocess(self,s_date=None, e_date=None):

            # generate prediction from the velocity model
            # mask out bad values
            self.offset_filter(s_date,e_date)

            # mask and filter
            #self.run_MaskAndFilter_track()

            # remove the bad values

            # interpolation
             
            return 0

        def gc_lon_lat_axis(self):

            doc = self.doc

            gc_losfile = os.path.join(self.trackfolder,self.geometry,'gc_los_offset_' + str(doc.runid) + '.rdr')

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

        def point_index(self,point):

            doc = self.doc

            if not (hasattr(doc,'lon_list') and hasattr(doc,'lat_list')):
                self.gc_lon_lat_axis()

            lon_list = doc.lon_list
            lat_list = doc.lat_list

            lon, lat = point

            if len(np.where(lon_list == lon)[0])==1:
                ind_x = np.where(lon_list == lon)[0][0]
            else:
                ind_x = None

            if len(np.where(lat_list == lat)[0])==1:
                ind_y = np.where(lat_list == lat)[0][0]
            else:
                ind_y = None

            return (ind_x, ind_y)


        def extract_offset_set_series(self, point_set=None, dates=None):
            # Strategy:
                # Open each offsetfield once, find offsets for all points.
                # Offsetfield outer loop.
                # Points inner loop.
        
            doc = self.doc

            ind_set = {}
            for point in point_set:
                ind_set[point] = self.point_index(point)

            #print(ind_set)
            #print(len(ind_set))

            # Initialization.
            pairs_set = {}
            offsets_set = {}
            for point in point_set:
                pairs_set[point] = []
                offsets_set[point] = []

            # Look at the offsetfields
            # Load in the los file.
            geo_losFile = os.path.join(self.trackfolder,self.geometry,'gc_los_offset_' + str(doc.runid) + '.rdr')
            print(geo_losFile)
            ds = gdal.Open(geo_losFile)
            losField = ds.GetRasterBand(1).ReadAsArray()

            #print('complete offsetfields: ', self.offsetfields)
            for offsetfield in self.offsetfields:

                self.doc_count = self.doc_count + 1
    
                date1str = offsetfield[0]
                date2str = offsetfield[1]

                date1 = offsetfield[2]
                date2 = offsetfield[3]

                # Only use dates in allowed range.
                if not ((date1 in dates) and (date2 in dates)):
                    continue

                title = date1str+'_'+date2str

                doc.interim = (date2-date1).days

                doc.offset_folder = os.path.join(self.trackfolder,self.offsetFolder,date1str+'_'+date2str)

                geo_rngOffsetFile = os.path.join(doc.offset_folder,'gc_filtRange_' + str(doc.runid) + '.off')
                geo_azOffsetFile = os.path.join(doc.offset_folder,'gc_filtAzimuth_' + str(doc.runid) + '.off')
                #print(geo_rngOffsetFile)

                # Read in offsetfield.
                ds = gdal.Open(geo_rngOffsetFile)
                rngOffsetField = ds.GetRasterBand(1).ReadAsArray()
                
                ds = gdal.Open(geo_azOffsetFile)
                azOffsetField = ds.GetRasterBand(1).ReadAsArray()

                #fig = plt.figure(1,figsize=(5,5))
                #ax = fig.add_subplot(111)
                #ax.imshow(rngOffsetField)
                #plt.show()

                for point in point_set:
                    ind_x = ind_set[point][0]
                    ind_y = ind_set[point][1]

                    # Available in this track offsetfield
                    if ind_x is not None and ind_y is not None:

                        rngOffset = rngOffsetField[ind_y, ind_x]
                        azOffset = azOffsetField[ind_y, ind_x]
                        los = losField[ind_y, ind_x]


                        # The value is valid.
                        if (not np.isnan(rngOffset) and not np.isnan(rngOffset) and los > 5):

                            # Convert offset from pixel to meter.
                            rngOffset = rngOffset * doc.rngPixelSize
                            azOffset = azOffset * doc.azPixelSize

                            # Comply to the los vector definition.
                            # Pointed from ground to satellite.
                            rngOffset = -rngOffset

                            # Save them
                            pairs_set[point].append([date1,date2])
                            offsets_set[point].append([rngOffset, azOffset])

                # End of this offsetfield.

            #print(pairs_set[point_set[0]])
            #print(offsets_set[point_set[0]])

            return (pairs_set, offsets_set)


        def extract_offset_series(self, point, dates):

            doc = self.doc
            ind_x, ind_y = self.point_index(point)

            pairs = []
            offsets = []

            # Load in the los file.
            geo_losFile = os.path.join(self.trackfolder,self.geometry,'gc_los_offset_' + str(doc.runid) + '.rdr')
            print(geo_losFile)
            ds = gdal.Open(geo_losFile)
            losField = ds.GetRasterBand(1).ReadAsArray()

            for offsetfield in self.offsetfields:

                self.doc_count = self.doc_count + 1
    
                date1str = offsetfield[0]
                date2str = offsetfield[1]

                date1 = offsetfield[2]
                date2 = offsetfield[3]

                # Only use dates in allowed range.
                if not ((date1 in dates) and (date2 in dates)):
                    continue

                title = date1str+'_'+date2str

                doc.interim = (date2-date1).days

                doc.offset_folder = os.path.join(self.trackfolder,self.offsetFolder,date1str+'_'+date2str)

                geo_rngOffsetFile = os.path.join(doc.offset_folder,'gc_filtRange_' + str(doc.runid) + '.off')
                geo_azOffsetFile = os.path.join(doc.offset_folder,'gc_filtAzimuth_' + str(doc.runid) + '.off')
                #print(geo_rngOffsetFile)

                ds = gdal.Open(geo_rngOffsetFile)
                rngOffsetField = ds.GetRasterBand(1).ReadAsArray()
              
                ds = gdal.Open(geo_azOffsetFile)
                azOffsetField = ds.GetRasterBand(1).ReadAsArray()


                if ind_x is not None and ind_y is not None:
                    
                    azOffset = azOffsetField[ind_y, ind_x]
                    rngOffset = rngOffsetField[ind_y, ind_x]
                    los = losField[ind_y, ind_x]

                    if (not np.isnan(rngOffset) and not np.isnan(rngOffset) and los > 5):

                        # Convert offset from pixel to meter.
                        rngOffset = rngOffset * doc.rngPixelSize
                        azOffset = azOffset * doc.azPixelSize

                        # Comply to the los vector definition.
                        # Pointed from ground to satellite.
                        rngOffset = -rngOffset

                        if not np.isnan(rngOffset) and not np.isnan(azOffset) and los>5:
                            #print(date1str+'_'+date2str, ': ', rngOffset, azOffset)
                            pairs.append([date1,date2])
                            offsets.append([rngOffset, azOffset])

            return (pairs, offsets)

