#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pickle

import isce
import isceobj

import gdal

from CSK_Rutford_Utils import CSK_Rutford_Utils

proj = "CSK-Evans"

def get_master(workdir, name):
    # Find master
    if proj == "CSK-Rutford":
        option = 1
    elif proj == "CSK-Evans":
        option = 2
    else:
        raise Exception()

    # choose the same one as the old setup
    if option == 1:
        f = open(os.path.join(workdir, name, "readme"))
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
    
    return master

if proj == "CSK-Rutford":
    tracklist = [8,10,23,25,40,52,55,67,69,82,97,99,114,126,128,129,141,143,156,158,171,172,173,186,188,201,203,215,218,230,231,232]
elif proj == "CSK-Evans":
    tracklist = range(22)
else:
    raise ValueError()

if proj == "CSK-Rutford":
    old_workdir="/net/kraken/nobak/mzzhong/CSK-Rutford"
    workdir="/net/kraken/nobak/mzzhong/CSK-Rutford-v2"
elif proj == "CSK-Evans":
    workdir = "/marmot-nobak/mzzhong/CSK-Evans-v3" 

# look through all tracks
for track in tracklist:
    print(track)
    track_name = "track_" + str(track).zfill(3) + '_0'

    if proj == "CSK-Rutford":
        # find the master
        master = get_master(old_workdir, track_name)
        print(master)
    elif proj == "CSK-Evans":
        master = get_master(workdir, track_name)
        print(master)

    # find the starting range of master
    master_data_pkl = os.path.join(workdir, track_name, "raw", master,"data.pkl")
    
    with open(master_data_pkl, 'rb') as f:
        master_data = pickle.load(f)

    #print(dir(master_data))
    startingRange = master_data.getStartingRange()
    print("starting range (m): ", startingRange)

    # find the range size of master
    master_slc = os.path.join(workdir, track_name, "raw", master, master + '.raw.slc.vrt')
    dataset = gdal.Open(master_slc)
    #print(dir(dataset))
    xsize = dataset.RasterXSize
    ysize = dataset.RasterYSize
    print(xsize, ysize)

    rangeLength = xsize
    print("range length (num of pixels): ", rangeLength)

    # find the range pixel size
    param_file = os.path.join(workdir, track_name, "raw", master, "params.txt")
    f = open(param_file,'r')
    rangePixelSize = float(f.readlines()[1].split(":")[1])
    f.close()

    print("range pixel size: ", rangePixelSize)

    # form range arr
    range_arr = startingRange + np.arange(rangeLength) * rangePixelSize
    range_arr = range_arr.reshape(1,rangeLength)
    print(range_arr)

    # read in local incidence angle
    los_file = os.path.join(workdir, track_name, "merged","geom_master","los.rdr.vrt")
    dataset = gdal.Open(los_file)
    inc = dataset.GetRasterBand(1).ReadAsArray()
    #print(inc)
    #print(inc.shape)

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

    # save the factor
    demfactor_file = os.path.join(workdir, track_name, "merged","geom_master","demfactor.rdr")

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
