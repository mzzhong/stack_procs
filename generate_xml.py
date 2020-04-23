#!/usr/bin/env python3
#
# Author: Heresh Fattahi
# Copyright 2016
#
import argparse
import isce
import isceobj
import os   
import gdal
import sys
            
def get_lat_lon(file):

    ds=gdal.Open(file)
    b=ds.GetRasterBand(1)
    width=b.XSize
    minLon = ds.GetGeoTransform()[0]
    deltaLon = ds.GetGeoTransform()[1]
    maxLat = ds.GetGeoTransform()[3]
    deltaLat = ds.GetGeoTransform()[5]
    minLat = maxLat + (b.YSize)*deltaLat

    return maxLat, deltaLat, minLon, deltaLon, width

def write_xml(outFile): 
    
    maxLat, deltaLat, minLon, deltaLon, width = get_lat_lon(outFile)
    
    unwImage = isceobj.Image.createImage()
    unwImage.setFilename(outFile)
    unwImage.setWidth(width)
    unwImage.bands = 1
    unwImage.scheme = 'BSQ'
    unwImage.dataType = 'SHORT'
    unwImage.setAccessMode('read')

    unwImage.coord2.coordDescription = 'Latitude'
    unwImage.coord2.coordUnits = 'degree'
    unwImage.coord2.coordStart = maxLat
    unwImage.coord2.coordDelta = deltaLat
    unwImage.coord1.coordDescription = 'Longitude'
    unwImage.coord1.coordUnits = 'degree'
    unwImage.coord1.coordStart = minLon
    unwImage.coord1.coordDelta = deltaLon
    
    unwImage.createImage()
    unwImage.renderHdr()


write_xml(sys.argv[1])

