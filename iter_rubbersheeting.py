#!/usr/bin/env python3
# Author: Minyan Zhong

import os
import datetime
import numpy as np
from dense_offset import dense_offset
import gdal
import pickle
import pathlib

# ISCE
import isce
import isceobj
from isceobj.Util.decorators import use_api
from isceobj.Util.Poly2D import Poly2D
import stdproc
from stdproc.stdproc import crossmul 

fmt='%Y%m%d'

class rubber_sheeting():

    def __init__(self, proj, track_number, date1, date2, mode):

        self.proj = proj
        self.track_number = track_number
        self.date1 = date1
        self.date2 = date2

        self.date1str = date1.strftime(fmt)
        self.date2str = date2.strftime(fmt)

        # case mode
        # pair: pair based rb, one resampling
        # track: stack based rb, two sampling
        self.mode = mode

        if self.proj == "CSK-Rutford":
            self.stack = "stripmap"
            self.projfolder = "/net/kraken/nobak/mzzhong/CSK-Rutford"

            self.trackname = "track_" + str(self.track_number).zfill(3) + '_0'
            self.pairname = self.date1str + '_' + self.date2str

            self.trackfolder = os.path.join(self.projfolder, "track_"+str(self.track_number).zfill(3)+"_0")
            self.pairfolder = os.path.join(self.projfolder, 
                                    "track_"+str(self.track_number).zfill(3)+"_0", 
                                    'pairs', self.pairname)

            if self.mode=="pair":
                self.runfolder = self.pairfolder
            elif self.mode=="track":
                self.runfolder = self.trackfolder

            # Interferogram folder
            self.ifgfolder = os.path.join(self.runfolder, "Igrams", self.pairname)
            pathlib.Path(self.ifgfolder).mkdir(parents=True, exist_ok=True)

        # Locate the master and slave SLC
        # symbolic link for pair mode
        if self.mode=="pair":
            self.master_SLC = os.path.join(self.runfolder, 'raw', self.date1str, self.date1str+ ".raw.slc")
            self.master_SLC_data_folder = os.path.join(self.runfolder, 'raw', self.date1str)
    
            self.slave_SLC = os.path.join(self.runfolder, 'raw', self.date2str, self.date2str+ ".raw.slc")
            self.slave_coregSLC = os.path.join(self.runfolder, 'coregSLC','Coarse',self.date2str, self.date2str+ ".slc")

            self.slave_SLC_data_folder = os.path.join(self.runfolder, 'raw', self.date2str)

        elif self.mode == "track":

            self.master_SLC = os.path.join(self.runfolder, 'coregSLC','Coarse',self.date1str, self.date1str+ ".slc")
            self.master_SLC_data_folder = os.path.join(self.runfolder, 'raw', self.date1str)
    
            self.slave_SLC = os.path.join(self.runfolder, 'coregSLC','Coarse',self.date2str, self.date2str+ ".slc")
            self.slave_SLC_data_folder = os.path.join(self.runfolder, 'raw', self.date2str)

        print(self.master_SLC)
        print(self.slave_SLC)

        # Set up offset field parameters
        self.runid = 20190921

    def getShape(self, gdalfile):
    
        dataset = gdal.Open(gdalfile, gdal.GA_ReadOnly)
        return dataset.RasterYSize, dataset.RasterXSize

    def prepare(self):

        readme_file = os.path.join(self.runfolder, "readme")
        rawfolder = self.runfolder+'/raw'
        redo=1
        if not os.path.exists(rawfolder) or redo==1:
            pathlib.Path(rawfolder).mkdir(exist_ok=True)
            date1raw = rawfolder + '/' + self.date1str
            date2raw = rawfolder + '/' + self.date2str

            cmd=" ".join(("ln -s", self.trackfolder + '/raw/' + self.date1str, date1raw))
            #print(cmd)
            if not os.path.exists(date1raw):
                os.system(cmd)
            cmd=" ".join(("ln -s", self.trackfolder + '/raw/' + self.date2str, date2raw))
            #print(cmd)
            if not os.path.exists(date2raw):
                os.system(cmd)

            # Modify master VRT file

        if not os.path.exists(self.runfolder+'/run_files'):
            cmd = 'stackStripMap.py -s {workdir}/raw --bbox "-78.3 -73.5 -85 -65" -d /net/kraken/nobak/mzzhong/CSK-Rutford/DEM/Rutford_surface.dem -w {workdir} -m {master} -W slc --useGPU'.format(workdir=self.runfolder, master=self.date1str)
            f = open(readme_file,"w")
            f.write(cmd)
            f.close()
            os.system("cat " + readme_file + "|bash")
        else:
            print("skip readme")

        if not os.path.exists(self.runfolder + "/raw_crop"):
            os.system("ln -s " + self.runfolder+"/raw" + " " + self.runfolder + "/raw_crop")

        # Run topo and geo2rdr'
        if not os.path.exists(self.slave_coregSLC):
            os.system("bash " + self.runfolder + "/run_files/run_2_master")
            os.system("bash " + self.runfolder + "/run_files/run_4_geo2rdr_coarseResamp")
        else:
            print("skip topo and geo2rdr")

    def run(self, iter_num, exe=False):

        self.iter_num = iter_num

        exe = True

        redo_offset = 0
        redo_filter = 1
        redo_rb = 0
        redo_resamp = 0
        redo_ifg = 0
        redo_geo = 0

        # Preparation
        workdir = self.projfolder
        offset = dense_offset(stack=self.stack, workdir = workdir, nproc=1, runid=self.runid, exe=exe)

        if self.proj == "CSK-Rutford":
            trackname = self.trackname
            pairname = self.pairname
        else:
            raise Exception("Undefined")

        # 0. Initialize offset object
        offset.initiate_track_pair(trackname=trackname, pairname = pairname, mode=self.mode)

        rb_suffix = "rb_" + str(self.iter_num)

        self.rb_suffix = rb_suffix
        offset.initiate_rb_pair(self.iter_num, rb_suffix)

        # 1. Run offset field
        offset.run_offset_pair_iter(iter_num=self.iter_num, exe=exe)

        # 2. Run offset field filter
        offset.offset_filter_pair_iter(iter_num=self.iter_num, redo=redo_filter, exe=exe)
        self.offset_rng_mis = offset.rng_mis

        # 3. Run rubber sheeting
        azOffset = offset.azOffsetName
        rngOffset = offset.rngOffsetName

        if self.mode == "pair":

            # Rubber sheeting

            if self.iter_num==0:
                azOffsetBefore = os.path.join(self.runfolder, "offsets", self.date2str, "azimuth.off")
                rngOffsetBefore = os.path.join(self.runfolder, "offsets", self.date2str, "range.off")
            else:
                azOffsetBefore = os.path.join(self.ifgfolder, "azimuth" + "_" + "rb_" + str(self.iter_num-1) +".off")
                rngOffsetBefore = os.path.join(self.ifgfolder, "range" + "_" + "rb_" + str(self.iter_num-1) +".off")
    
            if not os.path.exists(azOffsetBefore) or not os.path.exists(rngOffsetBefore):
                raise Exception("geometry files are missing: ",azOffsetBefore, rngOffsetBefore)

            azOffsetResamp = os.path.join(self.ifgfolder, "azimuth_off_resamp"+ "_" + rb_suffix + ".off")
            rngOffsetResamp = os.path.join(self.ifgfolder, "range_off_resamp" + "_" + rb_suffix + ".off")
    
            azOffsetNow = os.path.join(self.ifgfolder, "azimuth" + "_" + rb_suffix +".off")
            rngOffsetNow = os.path.join(self.ifgfolder, "range" + "_" + rb_suffix +".off")
    
            print(azOffsetBefore)
            print(rngOffsetBefore)
            print(azOffsetResamp)
            print(rngOffsetResamp)
            print(azOffsetNow)
            print(rngOffsetNow)
     
            if not os.path.exists(azOffsetNow) or redo_rb == 1:
                self._rubbersheet_offsetfields(azOffset, azOffsetBefore, azOffsetResamp, azOffsetNow)
            if not os.path.exists(rngOffsetNow) or redo_rb ==1:
                self._rubbersheet_offsetfields(rngOffset, rngOffsetBefore, rngOffsetResamp, rngOffsetNow)

        elif self.mode == "track":
            # Simply resample the obtained offset field
            azOffsetNow = os.path.join(self.ifgfolder, "azimuth" + "_" + rb_suffix +".off")
            rngOffsetNow = os.path.join(self.ifgfolder, "range" + "_" + rb_suffix +".off")

            if not os.path.exists(azOffsetNow) or redo_rb == 1:
                self._resample_offsetfields(azOffset, azOffsetNow)
            if not os.path.exists(rngOffsetNow) or redo_rb ==1:
                self._resample_offsetfields(rngOffset, rngOffsetNow)

        # 4. Resampling of SLC
        master_SLC = self.master_SLC
        slave_SLC = self.slave_SLC
        #resampled_slave_SLC = os.path.join(self.runfolder, "coregSLC", "Coarse", self.date2str, self.date2str + '_' + rb_suffix + '.slc')
        resampled_slave_SLC = os.path.join(self.ifgfolder, self.date2str + '_' + rb_suffix + '.slc')

        if not os.path.exists(resampled_slave_SLC) or redo_resamp==1:
            self._resample_SLC(master_SLC, slave_SLC, azOffsetNow, rngOffsetNow, resampled_slave_SLC)
        else:
            print("Existing resampled slave slc: ", resampled_slave_SLC)

        # 5. Form interferform
        for i in [0]:

            master_SLC = self.master_SLC
            if i==0:
                slave_SLC = self.slave_coregSLC
                ifgname = os.path.join(self.ifgfolder, self.pairname)
                ifg = ifgname + '.int'    
                filt_ifg = os.path.join(self.ifgfolder, "filt_" + self.pairname + ".int")
                filt_cor = os.path.join(self.ifgfolder, "filt_" + self.pairname + ".cor")
                unwname = os.path.join(self.ifgfolder, "filt_" + self.pairname)
                masterdir = os.path.join(self.runfolder, "raw_crop", self.date1str, "data")
 
            else:
                slave_SLC = resampled_slave_SLC
                ifgname = os.path.join(self.ifgfolder, self.pairname + '_' + rb_suffix)
                ifg = ifgname + '.int'    
                filt_ifg = os.path.join(self.ifgfolder, "filt_" + self.pairname + "_" + rb_suffix + ".int")
                filt_cor = os.path.join(self.ifgfolder, "filt_" + self.pairname + "_" + rb_suffix + ".cor")
                unwname = os.path.join(self.ifgfolder, "filt_" + self.pairname + "_" + rb_suffix)
                masterdir = os.path.join(self.runfolder, "raw_crop", self.date1str, "data")
  
            print("ifg master: ", master_SLC)
            print("ifg slave: ", slave_SLC)
 
            # Set unw file   
            unwfile = unwname + "_snaphu.unw"
            print(unwfile)
    
            program1 = "/net/jokull/state/partition1/home/mzzhong/dev/stackProcs/20180618/stripmap_stack/crossmul.py"
#            slave_SLC2 = "/net/kraken/nobak/mzzhong/CSK-Rutford/track_201_0/pairs/20130901_20130902/coregSLC/Coarse/20130902/20130902.slc"
#            cmdd = " ".join([program1, "-m", master_SLC, "-s", slave_SLC2, "-o", ifgname+'_2', "-a", "10", "-r", "10"])
#            print(cmdd)
#            os.system(cmdd)
#            #print(slave_SLC)
#            #print(slave_SLC2)
#            print(stop)
    
            program2 = "/net/jokull/state/partition1/home/mzzhong/dev/stackProcs/20180618/stripmap_stack/FilterAndCoherence.py"
            program3 = "/net/jokull/state/partition1/home/mzzhong/dev/stackProcs/20180618/stripmap_stack/unwrap.py"
    
            if not os.path.exists(unwfile) or redo_ifg == 1:
                cmd1 = " ".join([program1, "-m", master_SLC, "-s", slave_SLC, "-o", ifgname, "-a", "10", "-r", "10"])
                print(cmd1)
                os.system(cmd1)
    
                cmd2 = " ".join([program2, "-i", ifg, "-f", filt_ifg, "-c", filt_cor, "-s", "0.8"])
                print(cmd2)
                os.system(cmd2)
    
                cmd3 = " ".join([program3, "-i", filt_ifg, "-c", filt_cor, "-u", unwname, "-s", masterdir, "-d", "100", "-m", "snaphu"])
                print(cmd3)
                os.system(cmd3)
            else:
                print("Unwrapped file exists")

        # 6. Convert interferograms radian to meters
        # Get wavelength

        # Load master object
        with open(os.path.join(self.master_SLC_data_folder, 'data.pkl'), 'rb') as f:
            master = pickle.load(f)

        wvl = master.getInstrument().getRadarWavelength()
        geofile = unwname + "_snaphu_m.unw"
        print(geofile)

        if not os.path.exists(geofile) or redo_geo == 1:
            cmd4 = '''imageMath.py --eval="a_0;a_1/2.0/PI*{wvl}/2.0" --a={unw} -t float -s BIL -o {geofile}'''.format(wvl = wvl,unw = unwfile, geofile=geofile)
            print(cmd4)
            os.system(cmd4)
        else:
            print("Geofile exists")

        # 7. Get offset field
        print(dir(offset.doc))
        print("pixelSize: ", offset.doc.rngPixelSize)
        raw_offset = rngOffset
        raw_offset_coregis = os.path.join(self.ifgfolder, "raw_range_offset_m" + '.off')

        # Deal with raw offset
        # Flip the sign because the definition of los vector
        cmd5 = '''imageMath.py --eval="-(a_0-{rng_mis})*{pixelSize}" --a={raw_offset} -t float -s BIL -o {raw_offset_coregis}'''.format(rng_mis=self.offset_rng_mis, \
                                                                                                                raw_offset=raw_offset, \
                                                                                                                raw_offset_coregis=raw_offset_coregis,\
                                                                                                                  pixelSize=offset.doc.rngPixelSize)
        # Do it
        print(cmd5)
        os.system(cmd5)

        # Deal with fine offset
        fine_offset = geofile
        dataset = gdal.Open(fine_offset, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(2).ReadAsArray()
        fine_offset_coregis = os.path.join(self.ifgfolder, "fine_range_offset_m" + '.off')

        if self.track_number==201:
            miscoregis = np.nanmedian(data[30:1300,800:1800])
        elif self.track_number==128:
            miscoregis = np.nanmedian(data[20:200,50:250])
        else:
            raise Exception("Undefined")

        print(miscoregis)

        cmd6 = '''imageMath.py --eval="a_1-{miscoregis}" --a={fine_offset} -t float -s BIL -o {fine_offset_coregis}'''.format(miscoregis=miscoregis, \
                                                                                                                fine_offset=fine_offset, \
                                                                                                                fine_offset_coregis=fine_offset_coregis)
        ## Do it
        print(cmd6)
        os.system(cmd6)

        return 0

    def _resample_offsetfields(self, Offset, reference, resampledOffset):
    
        length, width = self.getShape(reference)
        print('oversampling the filtered and masked offsets to the width and length:', width, ' ', length )
        cmd = 'gdal_translate -of ENVI  -outsize  ' + str(width) + ' ' + str(length) + ' ' + Offset + ' ' + resampledOffset
        os.system(cmd)
    
        img = isceobj.createImage()
        img.setFilename(resampledOffset)
        img.setWidth(width)
        img.setLength(length)
        img.setAccessMode('READ')
        img.bands = 1
        img.dataType = 'FLOAT'
        img.scheme = 'BIP'
        #img.createImage()
        img.renderHdr()
        img.renderVRT()
        #img.finalizeImage()

        return 0

    def _rubbersheet_offsetfields(self, Offset, OffsetBefore, resampledOffset, outName):
    
        length, width = self.getShape(OffsetBefore)
        print('oversampling the filtered and masked offsets to the width and length:', width, ' ', length )
        cmd = 'gdal_translate -of ENVI  -outsize  ' + str(width) + ' ' + str(length) + ' ' + Offset + ' ' + resampledOffset
        os.system(cmd)
    
        img = isceobj.createImage()
        img.setFilename(resampledOffset)
        img.setWidth(width)
        img.setLength(length)
        img.setAccessMode('READ')
        img.bands = 1
        img.dataType = 'FLOAT'
        img.scheme = 'BIP'
        #img.createImage()
        img.renderHdr()
        img.renderVRT()
        #img.finalizeImage()
    
        print ('Adding the dense offsets to the offsets Before. Output: ', outName)
        
        # option 1
        cmd = "imageMath.py -e='a+b' -o " + outName + " -t float  --a=" + OffsetBefore + " --b=" + resampledOffset + " -s BIP -t float"
        print(cmd)
        os.system(cmd)

        # option 2
        # Using gdal_calc seems faster
        # gdal_calc.py -A 20091219_20100321_usingCoarseCoreg/filtRngOff_resampled.bil -B ../offsets/20100321/range.off --outfile=range_sum.off --calc="A+B" --format=ENVI --type=Float32
        #cmd = "gdal_calc.py -A " + OffsetBefore + " -B " + resampledOffset + " --outfile=" + outName + ' --calc="A+B" --format=ENVI --type=Float32'
        #print (cmd)
        #os.system(cmd)
    
        return 0

    def _resample_SLC(self, master_slc, slave_slc, azOffset, rngOffset, resampled_slave_slc):

        # Old way
#        dirname = os.path.dirname(slave_slc)
#        slaveShelveDir = os.path.join(dirname, "slaveShelve")
#        masterShelveDir = os.path.join(dirname, "masterShelve")
#
#        with open(os.path.join(slaveShelveDir, 'data.pkl'), 'rb') as f:
#            slave = pickle.load(f)
#            doppler = slave._dopplerVsPixel
#    
#        azpoly = None
#        rgpoly = None
    
#        with open(os.path.join(masterShelveDir, 'data.pkl'), 'rb') as f:
#            master = pickle.load(f)

        with open(os.path.join(self.slave_SLC_data_folder, 'data.pkl'), 'rb') as f:
            slave = pickle.load(f)
            doppler = slave._dopplerVsPixel

        azpoly = None
        rgpoly = None

        # Load master object
        with open(os.path.join(self.master_SLC_data_folder, 'data.pkl'), 'rb') as f:
            master = pickle.load(f)


#########################################################
        # Below is adapted from Heresh's resampleSlc
        #runResamp(slave, inps.offsets, outfile, doppler, azpoly,rgpoly, 
        #    flatten=(not inps.noflat), zero=inps.zero, dims = inps.dims, master = master)
        
#        burst = slave
#        offdir = inps.offsets
#        outname = outfile
#        doppler = doppler
#        azpoly = azpoly
#        rgpoly = rgpoly
#        master = master
#        flatten = True
#        zero = False
#        dims = None 

        rngImg = isceobj.createImage()
        rngImg.load(rngOffset + '.xml')
        rngImg.setAccessMode('READ')

        aziImg = isceobj.createImage()
        aziImg.load(azOffset + '.xml')
        aziImg.setAccessMode('READ')

        width = rngImg.getWidth()
        length = rngImg.getLength()
       
        inimg = isceobj.createSlcImage()
        inimg.load(slave.getImage().filename + '.xml')
        inimg.setAccessMode('READ')

        prf = slave.getInstrument().getPulseRepetitionFrequency()

        # Zero is False, factor = 1
        zero = False
        if zero:
            factor = 0.0
        else:
            factor = 1.0

        try:
            print('Polynomial doppler provided')
            coeffs = [factor * 2*np.pi*val/prf for val in doppler._coeffs]
        except:
            print('List of coefficients provided')
            coeffs = [factor * 2*np.pi*val/prf for val in doppler]
    
        zcoeffs = [0. for val in coeffs]
    
        print("coeffs: ", coeffs)
        print("zcoeffs: ", zcoeffs)

        # Create doppler object
        dpoly = Poly2D()
        dpoly.initPoly(rangeOrder=len(coeffs)-1, azimuthOrder=0, coeffs=[coeffs])

        # Create resampling object
        rObj = stdproc.createResamp_slc()
        rObj.slantRangePixelSpacing = slave.getInstrument().getRangePixelSize()
        rObj.radarWavelength = slave.getInstrument().getRadarWavelength()
        rObj.dopplerPoly = dpoly
        rObj.azimuthOffsetsPoly = None
        rObj.rangeOffsetsPoly = None

        rObj.imageIn = inimg
        rObj.flatten = True
        rObj.outputWidth = width
        rObj.outputLines = length
        rObj.residualRangeImage = rngImg
        rObj.residualAzimuthImage = aziImg

        rObj.startingRange = slave.startingRange
        rObj.referenceStartingRange = master.startingRange
        rObj.referenceSlantRangePixelSpacing = master.getInstrument().getRangePixelSize()
        rObj.referenceWavelength = master.getInstrument().getRadarWavelength()


        imgOut = isceobj.createSlcImage()
        imgOut.setWidth(width)
        imgOut.filename = resampled_slave_slc
        imgOut.setAccessMode('write')

        rObj.resamp_slc(imageOut=imgOut)
        imgOut.renderHdr()

        return 0

    def geocode(self, iter_num):

        ## Geocode ###
        ## Generate lat lon files
        lonfile = os.path.join(self.runfolder, "merged", "geom_master","lon.rdr")
        ml_lonfile = os.path.join(self.runfolder, "merged", "geom_master","lon_a10_r10.rdr")

        latfile = os.path.join(self.runfolder, "merged", "geom_master","lat.rdr")
        ml_latfile = os.path.join(self.runfolder, "merged", "geom_master","lat_a10_r10.rdr")

        losfile = os.path.join(self.runfolder, "merged", "geom_master","los.rdr")
        ml_losfile = os.path.join(self.runfolder, "merged", "geom_master","los_a10_r10.rdr")

        if not os.path.exists(ml_lonfile) or not os.path.exists(ml_latfile):
            ## Multi-look
            cmd1 = "looks.py -i {full} -r 10 -a 10 -o {downsample}".format(full=lonfile, downsample=ml_lonfile)
            cmd2 = "looks.py -i {full} -r 10 -a 10 -o {downsample}".format(full=latfile, downsample=ml_latfile)
            #cmd3 = "looks.py -i {full} -r 10 -a 10 -o {downsample}".format(full=losfile, downsample=ml_losfile)
            os.system(cmd1)
            os.system(cmd2)

        offset_data = os.path.join(self.ifgfolder, "fine_range_offset_m" + '.off')
        gc_offset_data = os.path.join(self.ifgfolder, "gc_fine_range_offset_m" + '.off')
        lon_step = 0.005
        lat_step = 0.001

        if self.track_number==201:
            bbox = '"' + ' '.join(["-78.88","-78.3","-85.5","-79.7"]) + '"'
        elif self.track_number==128:
            bbox = '"' + ' '.join(["-78.82","-78.32","-84.7","-80.4"]) + '"'
        else:
            raise Exception("Undefined")


        gc_cmd = "geocodeGdal.py -l {latfile} -L {lonfile} -x {lon_step} -y {lat_step} -f {infile} -b {bbox_SNWE}".format(\
        latfile = ml_latfile,\
        lonfile = ml_lonfile,\
        lon_step = str(lon_step),\
        lat_step = str(lat_step),\
        infile = offset_data,\
        bbox_SNWE = bbox)

        print(offset_data)
        print(gc_cmd)
        os.system(gc_cmd)

def main():

    proj = "CSK-Rutford"
    track_number = 128

    if track_number == 201:
        gs = []
        gs.append((datetime.date(2013,9,1),datetime.date(2013,9,2)))
        gs.append((datetime.date(2013,10,19),datetime.date(2013,10,20)))
        gs.append((datetime.date(2013,11,20),datetime.date(2013,11,21))) # bad
        gs.append((datetime.date(2013,12,6),datetime.date(2013,12,7)))
        gs.append((datetime.date(2013,12,22),datetime.date(2013,12,23)))
        gs.append((datetime.date(2014,1,23),datetime.date(2014,1,24)))
        gs.append((datetime.date(2014,8,3),datetime.date(2014,8,4)))

    if track_number==128:
        gs=[]
        gs.append((datetime.date(2014,7,29),datetime.date(2014,7,30)))
        gs.append((datetime.date(2014,8,14),datetime.date(2014,8,15)))
 

    for g in gs:
        date1, date2 = g
        mode = "pair"
        rs = rubber_sheeting(proj=proj, track_number = track_number, date1=date1, date2=date2, mode=mode)
    
        iter_num = 0
        rs.prepare()
        rs.run(iter_num)
        rs.geocode(iter_num)


if __name__=='__main__':
    main()
