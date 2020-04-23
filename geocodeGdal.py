#!/usr/bin/env python3
########################
#Author: Heresh Fattahi
#Copyright 2016
# Modified by Minyan Zhong
######################
import argparse
import isce
import isceobj
import os
import gdal
import xml.etree.ElementTree as ET

import numpy as np

import subprocess
import time

def createParser():
    '''
    Create command line parser.
    '''

    parser = argparse.ArgumentParser( description='Create DEM simulation for merged images')

    parser.add_argument('-l','--lat', dest='latFile', type=str, required=True,
            help = 'latitude file in radar coordinate')

    parser.add_argument('-L','--lon', dest='lonFile', type=str, required=True,
            help = 'longitude file in radar coordinate')

    parser.add_argument('-f', '--filelist', dest='prodlist', type=str, required=True,
            help='Input file to be geocoded')

    parser.add_argument('-b', '--bbox', dest='bbox', type=str, default="",
            help='Bounding box (SNWE)')

    parser.add_argument('-x', '--lon_step', dest='lonStep', type=str, default=0.001,
            help='output pixel size (longitude) in degrees. Default 0.001')

    parser.add_argument('-y', '--lat_step', dest='latStep', type=str, default=0.001,
            help='output pixel size (latitude) in degrees. Default 0.001')

    parser.add_argument('-o', '--xoff', dest='xOff', type=int, default=0,
            help='Offset from the begining of geometry files in x direction. Default 0.0')

    parser.add_argument('-p', '--yoff', dest='yOff', type=int, default=0,
            help='Offset from the begining of geometry files in y direction. Default 0.0')

    parser.add_argument('-r', '--resampling_method', dest='resamplingMethod', type=str, default='near',
            help='Resampling method (gdalwarp resamplin methods)')

    parser.add_argument('--srcnodata', dest='srcnodata', type=str, default="0",
            help='invalid data in src file')

    parser.add_argument('--dstnodata', dest='dstnodata', type=str, default="0",
            help='invalid data in dst file')

    parser.add_argument('--prefix', dest='prefix', type=str, default="gc",
            help='prefix of geocoded file')

    parser.add_argument('--suffix', dest='suffix', type=str, default="",
            help='suffix of geocoded file')

    return parser

def cmdLineParse(iargs = None):
    '''
    Command line parser.
    '''

    parser = createParser()
    inps =  parser.parse_args(args = iargs)

    if inps.bbox != "":
        inps.bbox = [val for val in inps.bbox.split()]
        if len(inps.bbox) != 4:
            raise Exception('Bbox should contain 4 floating point values')


    inps.prodlist = inps.prodlist.split()
    return inps

def prepare_lat_lon(inps):

    latFile = os.path.abspath(inps.latFile)
    lonFile = os.path.abspath(inps.lonFile)
    cmd = 'isce2gis.py vrt -i ' + latFile
    #os.system(cmd)
    subprocess.call(cmd, shell=True)
    cmd = 'isce2gis.py vrt -i ' + lonFile
    #os.system(cmd)
    subprocess.call(cmd, shell=True)

    success = False
    while not success:
        try:   
            width, length =  getSize(latFile)
            widthFile , lengthFile = getSize(inps.prodlist[0])
        
            xOff = inps.xOff
            yOff = inps.yOff
    
            tempLat = os.path.join(os.path.dirname(inps.prodlist[0]), 'tempLAT.vrt')
            tempLon = os.path.join(os.path.dirname(inps.prodlist[0]), 'tempLON.vrt')
    
            cmd = 'gdal_translate -of VRT -srcwin ' + str(xOff) + ' ' + str(yOff) \
                    +' '+ str(width - xOff) +' '+ str(length - yOff) +' -outsize ' + str(widthFile) + \
                    ' '+ str(lengthFile)  + ' -a_nodata 0 ' + latFile +'.vrt ' +  tempLat
    
            subprocess.call(cmd, shell=True)
    
            cmd = 'gdal_translate -of VRT -srcwin ' + str(xOff) + ' ' + str(yOff) \
                +' '+ str(int(width-xOff)) +' '+ str(int(length-yOff)) +' -outsize ' + str(widthFile) +\
                    ' '+ str(lengthFile)  + ' -a_nodata 0 ' + lonFile +'.vrt ' +  tempLon
    
            subprocess.call(cmd, shell=True)
    
            if os.path.exists(tempLat) and os.path.exists(tempLon):
                success = True
        except:
            time.sleep(1)

    return tempLat, tempLon

    # gdal_translate -of VRT -srcwin  384 384 64889 12785 -outsize 1013 199 ../../COMBINED/GEOM_MASTER/LAT.rdr LAT_off.vrt
    

def writeVRT(infile, latFile, lonFile):
#This function is modified from isce2gis.py
            latFile = os.path.abspath(latFile)
            lonFile = os.path.abspath(lonFile)
            infile = os.path.abspath(infile)
            cmd = 'isce2gis.py vrt -i ' + infile
            #os.system(cmd)
            subprocess.call(cmd, shell=True)

            tree = ET.parse(infile + '.vrt')
            root = tree.getroot()

            meta = ET.SubElement(root, 'metadata')
            meta.attrib['domain'] = "GEOLOCATION"
            meta.tail = '\n'
            meta.text = '\n    '


            #rdict = { 'Y_DATASET' : os.path.relpath(latFile , os.path.dirname(infile)),
            #          'X_DATASET' :  os.path.relpath(lonFile , os.path.dirname(infile)),

            rdict = { 'Y_DATASET' : latFile ,
                      'X_DATASET' :  lonFile ,
                      'X_BAND' : "1",
                      'Y_BAND' : "1",
                      'PIXEL_OFFSET': "0",
                      'LINE_OFFSET' : "0",
                      'LINE_STEP' : "1",
                      'PIXEL_STEP' : "1" }

            for key, val in rdict.items():
                data = ET.SubElement(meta, 'mdi')
                data.text = val
                data.attrib['key'] = key
                data.tail = '\n    '

            data.tail = '\n'
            tree.write(infile + '.vrt')

            #print(stop)


def runGeo(inps):

    success = False

    while not success:
        try:
            for rfile in inps.prodlist:
                cmd = 'isce2gis.py envi -i ' + rfile
                #os.system(cmd)
                subprocess.call(cmd, shell=True)
        
            latFile, lonFile = prepare_lat_lon(inps)
        
            if inps.bbox:
                print("BBOX set by user:")
                SNWE = str(inps.bbox[0]) + ', ' + str(inps.bbox[1]) + ', ' + str(inps.bbox[2]) + ', ' + str(inps.bbox[3])
                print("SNWE: ", SNWE)
                WSEN = str(inps.bbox[2]) + ' ' + str(inps.bbox[0]) + ' ' + str(inps.bbox[3]) + ' ' + str(inps.bbox[1])
                print("WSEN (input to gdalwarp): ", WSEN)
            else:
                print("Inferring BBOX from LAT LON file:")
                minlat, maxlat, minlon, maxlon = getBound(latFile,lonFile)
                # Round to grid points
                # Option 1
                #coe_lon = np.round(1/float(inps.lonStep))
                #coe_lat = np.round(1/float(inps.latStep))
            
                # Option 2
                # round to 0.1 degree
                coe_lon = 10
                coe_lat = 10
            
                minlat = np.ceil(minlat * coe_lat)/coe_lat
                minlon = np.ceil(minlon * coe_lon)/coe_lon
            
                maxlat = np.floor(maxlat * coe_lat)/coe_lat
                maxlon = np.floor(maxlon * coe_lon)/coe_lon
                SNWE = str(minlat) + ", " + str(maxlat) + ", " + str(minlon) + ", " +str(maxlon)
                print("SNWE: ", SNWE)
                WSEN = str(minlon) + ' ' + str(minlat) + ' ' + str(maxlon) + ' ' +str(maxlat)
                print("WSEN (input to gdalwarp): ", WSEN)
        
            for rfile in inps.prodlist:
               rfile = os.path.abspath(rfile)
               print ('geocoding ' + rfile)
               #cmd = 'isce2gis.py vrt -i '+ rfile + ' --lon ' + lonFile + ' --lat '+ latFile
               #os.system(cmd)
        
               writeVRT(rfile, latFile, lonFile)
        
               folder = os.path.dirname(rfile)
               filename = os.path.basename(rfile) 
        
               # Number of bands
               ds = gdal.Open(rfile)
               bands = ds.RasterCount

               filebasename, ending = filename.split(".")

               outfilename = inps.prefix + "_" + filebasename
               if inps.suffix!="":
                    outfilename = outfilename + "_" + inps.suffix

               outfilename = outfilename + "." + ending

               outname = os.path.join(folder, outfilename)

               if outname.split('.')[-1] in ['bil','bip','BIL','BIP'] and bands>1:
                  outname = ''.join(outname.split('.')[0:-1])+'.bsq'
               print("outname: ", outname)
               #return
        
               if bands == 1:
                  srcnodata = '"' + inps.srcnodata + '"'
                  dstnodata = '"' + inps.dstnodata + '"'
                  cmd = 'gdalwarp -of ENVI -geoloc  -te '+ WSEN + ' -tr ' + str(inps.lonStep) + ' ' + str(inps.latStep) + ' -srcnodata '+ srcnodata+ ' -dstnodata ' + dstnodata + ' -r ' +inps.resamplingMethod +' ' + rfile +'.vrt '+ outname
               elif bands == 2:
        
                  srcnodata = '"' + inps.srcnodata + ' '+ inps.srcnodata + '"'
                  dstnodata = '"' + inps.dstnodata + ' '+ inps.dstnodata + '"'
        
                  cmd = 'gdalwarp -of ENVI -geoloc  -te '+ WSEN + ' -tr ' + str(inps.lonStep) + ' ' + str(inps.latStep) + ' -srcnodata '+ srcnodata+ ' -dstnodata ' + dstnodata + ' -r ' +inps.resamplingMethod +' ' + rfile +'.vrt '+ outname
        
               elif bands == 3:
        
                  srcnodata = '"' + inps.srcnodata + ' '+ inps.srcnodata + ' ' + inps.srcnodata + '"'
                  dstnodata = '"' + inps.dstnodata + ' '+ inps.dstnodata + ' ' + inps.dstnodata + '"'
        
                  cmd = 'gdalwarp -of ENVI -geoloc  -te '+ WSEN + ' -tr ' + str(inps.lonStep) + ' ' + str(inps.latStep) + ' -srcnodata '+ srcnodata+ ' -dstnodata ' + dstnodata + ' -r ' +inps.resamplingMethod +' ' + rfile +'.vrt '+ outname
        
               print(cmd)
               #os.system(cmd)
               subprocess.call(cmd, shell=True)
               #time.sleep(2)
               write_xml(outname, bands)

            if os.path.exists(outname) and os.path.exists(outname + '.xml') and os.path.exists(outname + '.vrt'):
                success = True
        except:
            time.sleep(1)

def getSize(f):

    ds=gdal.Open(f, gdal.GA_ReadOnly)
    b=ds.GetRasterBand(1)
    width = b.XSize
    length = b.YSize
    ds = None
    return width, length

def getBound(latfile,lonfile): #added by Minyan Zhong
    
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

    return minlat,maxlat,minlon,maxlon

def get_lat_lon(f):

    ds=gdal.Open(f, gdal.GA_ReadOnly)
    b=ds.GetRasterBand(1)
    width  = b.XSize
    length = b.YSize
    minLon = ds.GetGeoTransform()[0]
    deltaLon = ds.GetGeoTransform()[1]
    maxLat = ds.GetGeoTransform()[3]
    deltaLat = ds.GetGeoTransform()[5]
    minLat = maxLat + (b.YSize)*deltaLat
    ds = None
    return maxLat, deltaLat, minLon, deltaLon, width, length

def write_xml(outFile,bands): 

    maxLat, deltaLat, minLon, deltaLon, width, length = get_lat_lon(outFile)

    unwImage = isceobj.Image.createImage()
    unwImage.setFilename(outFile)
    unwImage.setWidth(width)
    unwImage.setLength(length)
    unwImage.bands = bands

    # only work for los and offsetfields, need to change !
    if bands == 1:
        unwImage.scheme = 'BIP'
    elif bands == 2:
        unwImage.scheme = 'BSQ'
    elif bands == 3:
        unwImage.scheme = 'BSQ'
 
    unwImage.dataType = 'FLOAT'
    unwImage.setAccessMode('read')
    
    unwImage.coord2.coordDescription = 'Latitude'
    unwImage.coord2.coordUnits = 'degree'
    unwImage.coord2.coordStart = maxLat 
    unwImage.coord2.coordDelta = deltaLat 
    unwImage.coord1.coordDescription = 'Longitude'
    unwImage.coord1.coordUnits = 'degree'
    unwImage.coord1.coordStart = minLon 
    unwImage.coord1.coordDelta = deltaLon 

   # unwImage.createImage()
    unwImage.renderHdr()
    unwImage.renderVRT()

def main(iargs=None):
    '''
    Main driver.
    '''
    inps = cmdLineParse(iargs)
    runGeo(inps)
 
   
if __name__ == '__main__':
    main()


