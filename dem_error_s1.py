#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pickle

import isce
import isceobj

import gdal

import xml.etree.ElementTree as ET

def get_master(workdir, name):

    option = 1
    if option == 1:
        f = open(os.path.join(workdir, name, "readme"))
        line = f.readlines()[0]
        f.close()
        elements = line.split()
        for i, element in enumerate(elements):
            if element == "-m" or element == "--master":
                master = elements[i+1]
                break
    else:
        raise Exception()
    
    return master

#tracklist = [37, 52, 169, 65, 7, 50, 64]
tracklist = [52, 169, 65, 7, 50, 64]

proj = "S1-Evans-v2"
workdir="/net/kraken/nobak/mzzhong/S1-Evans-v2"

# Loop through all tracks
for track in tracklist:
    print(track)
    track_name = "track_" + str(track)
    # find the master
    master = get_master(workdir, track_name)
    print(master)

    # find the starting range of master
    master_data_xml = os.path.join(workdir, track_name, "master", "IW1.xml")

    print(master_data_xml)

    f = open(master_data_xml)
    lines = f.readlines()
    read_rangepixelsize = 0
    read_startingrange = 0
    for line in lines:
        if len(line.split())>1 and line.split()[1].startswith('name="rangepixelsize"'):
            read_rangepixelsize = 1
        elif len(line.split())>1 and line.split()[1].startswith('name="startingrange"'):
            read_startingrange = 1
        elif read_rangepixelsize:
            rangepixelsize = line.split('>')[1].split('<')[0]
            read_rangepixelsize = 0
            print(rangepixelsize)
        elif read_startingrange:
            startingrange = line.split('>')[1].split('<')[0]
            read_startingrange = 0
            print(startingrange)

    rangePixelSize = float(rangepixelsize)
    startingRange = float(startingrange)

    print("range pixel size (m): ", rangePixelSize)
    print("starting range (m): ", startingRange)

    # find the range size of master
    master_slc = os.path.join(workdir, track_name, 'merged/SLC', master, master + '.slc.full')
    dataset = gdal.Open(master_slc)
    print("done")
    xsize = dataset.RasterXSize
    ysize = dataset.RasterYSize
    print("Xsize: ", xsize)
    print("Ysize: ", ysize)

    rangeLength = xsize
    print("range length (num of pixels): ", rangeLength)

    # Form range arr
    range_arr = startingRange + np.arange(rangeLength) * rangePixelSize
    range_arr = range_arr.reshape(1, rangeLength)
    print("range_arr: ", range_arr)
    
    # Read in local incidence angle
    los_file = os.path.join(workdir, track_name, "merged","geom_master","los.rdr.full")
    print(los_file)
    dataset = gdal.Open(los_file)
    inc = dataset.GetRasterBand(1).ReadAsArray()
    print(inc.shape)

    demfactor = np.zeros(shape = inc.shape)

    for y in range(ysize):
        print(y, ysize)
        inv_demfactor = range_arr * np.sin(inc[y, :] / 180 * np.pi)
        inv_demfactor[np.isnan(inv_demfactor)] = np.nan
        demfactor[y, :] = 1 / inv_demfactor

    demfactor_2 = demfactor[::100,::100]
    plot_factor = 1
    if plot_factor:
        fig = plt.figure(1, figsize=(10,10))
        fig.clf()
        ax = fig.add_subplot(111)
        cs = ax.imshow(demfactor_2, cmap='jet')
        fig.colorbar(cs, shrink=0.7)
        fig.savefig("demfactor_" + track_name + ".png")

    # Save the factor
    demfactor_file = os.path.join(workdir, track_name, "merged","geom_master","demfactor.rdr.full")

    driver = gdal.GetDriverByName( 'ENVI' )
    nbands = 1
    bandType = 7

    dst_ds = driver.Create(demfactor_file, demfactor.shape[1], demfactor.shape[0], nbands, bandType)
    dst_ds.GetRasterBand(1).WriteArray( demfactor, 0 ,0 )
    dst_ds = None

    outImg = isceobj.createImage()
    outImg.setDataType('DOUBLE')
    outImg.setFilename(demfactor_file)
    outImg.setBands(1)
    outImg.scheme = 'BIP'
    outImg.setWidth(demfactor.shape[1])
    outImg.setLength(demfactor.shape[0])
    outImg.setAccessMode('read')
    outImg.renderHdr()
    outImg.renderVRT()

    del inc
    del demfactor
