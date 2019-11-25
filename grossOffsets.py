#!/usr/bin/env python3
# Generate grossOffsets (pixel) based on velocity field.
# Author: Minyan Zhong


import numpy as np
import pyproj
import subprocess
import isce
import isceobj
from iscesys.Component.ProductManager import ProductManager as PM
import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import gdal

from scipy.interpolate import interp2d, griddata
from scipy import signal
from scipy import ndimage

import matplotlib.pyplot as plt

import os

class grossOffsets:

    def __init__(self, figpath=None, runid=None):
        
        self.nfig = 1
        self.figsize = (10,10)

        ## Antarctica Velocity File
        self.vel_file = '/net/jokull/nobak/mzzhong/Ant_Data/velocity_models/antarctica_ice_velocity_900m.nc'
        self.vel_file_v2 = '/home/mzzhong/jokull-nobak-net/Ant_Data/velocity_models/antarctica_ice_velocity_450m_v2.nc'

        self.vProj = pyproj.Proj('+init=EPSG:3031')

        self.figpath = figpath
        self.runid = runid

    def setMode(self,mode):

        if mode == 'interior' or mode == 'exterior':
            self.mode = mode
        else:
            raise Exception('Wrong gross offset mode')

    def setLatFile(self,val):
        self.latfile = val

    def setLonFile(self,val):
        self.lonfile = val

    def setLosFile(self,val):
        self.losfile = val

    def setXSize(self,val):
        self.XSize = val

    def setYize(self,val):
        self.YSize = val

    def setMargin(self,val):
        self.margin = val

    def setWinSizeHgt(self,val):
        self.winSizeHgt = val

    def setWinSizeWidth(self,val):
        self.winSizeWidth = val

    def setSearchSizeHgt(self,val):
        self.searchSizeHgt = val

    def setSearchSizeWidth(self,val):
        self.searchSizeWidth = val

    def setSkipSizeHgt(self,val):
        self.skipSizeHgt = val

    def setSkipSizeWidth(self,val):
        self.skipSizeWidth = val

    # exterior mode
    def setOffsetLat(self,lat):
        self.lat = lat

    def setOffsetLon(self,lon):
        self.lon = lon

    def setOffsetInc(self,inc):
        self.inc = inc

    def setOffsetAzi(self,azi):
        self.azi = azi

    def setNumWinDown(self,numWinDown):
        self.numWinDown = numWinDown

    def setNumWinAcross(self,numWinAcross):
        self.numWinAcross = numWinAcross

    def setbTemp(self,val):
        self.bTemp = val

    def setPixelSize(self,azPixelSize,rngPixelSize):
        
        self.azPixelSize = azPixelSize
        self.rngPixelSize = rngPixelSize

    def get_veloData(self):
        
        print("getting velocity data...")
        fh=Dataset(self.vel_file,mode='r')
        self.vx = fh.variables['vx'][:]
        self.vy = fh.variables['vy'][:]
        self.vx = np.flipud(self.vx)
        self.vy = np.flipud(self.vy) 
        self.v = np.sqrt(np.multiply(self.vx,self.vx)+np.multiply(self.vy,self.vy))
        print(self.v.shape)

        # Set up the coordinates.
        self.x0 = np.arange(-2800000,2800000,step=900)
        self.y0 = np.arange(-2800000,2800000,step=900)+200

        #x,y = np.meshgrid(self.x0,self.y0)

        #self.AntVeloDataMap = Basemap(width=5600000,height=5600000,\
        #                        resolution='l',projection='stere',\
        #                        lat_ts=-71,lat_0=-90,lon_0=0)

        #self.vel_lon, self.vel_lat= self.vProj(x,y,inverse="true")


    def get_veloData_v2(self):
        
        print("getting velocity data...")
        fh=Dataset(self.vel_file_v2,mode='r')
        
        #self.vx = fh.variables['VX'][:]
        #self.vy = fh.variables['VY'][:]
        
        #self.vx = np.flipud(self.vx)
        #self.vy = np.flipud(self.vy) 
        
        #self.v = np.sqrt(np.multiply(self.vx,self.vx)+np.multiply(self.vy,self.vy))
        #print(self.v.shape)

        # Set up the coordinates.
        #self.x0 = fh.variables['x'][:]
        #self.y0 = fh.variables['y'][:]

        #print(self.x0)
        #print(self.y0)
        #print(stop)


        dataset = gdal.Open("/net/jokull/nobak/mzzhong/Ant_Data/data-visualization/Ant_Velo_VX.tif", gdal.GA_ReadOnly)
        band = dataset.GetRasterBand(1)
        self.vx = band.ReadAsArray()

        dataset = gdal.Open("/net/jokull/nobak/mzzhong/Ant_Data/data-visualization/Ant_Velo_VY.tif", gdal.GA_ReadOnly)
        band = dataset.GetRasterBand(1)
        self.vy = band.ReadAsArray()

        self.vx = np.flipud(self.vx)
        self.vy = np.flipud(self.vy) 

        self.v = np.sqrt(np.multiply(self.vx,self.vx)+np.multiply(self.vy,self.vy))

        self.x0 = np.arange(-2800000,2800000,step=450)
        self.y0 = np.arange(-2800000,2800000,step=450)+200


    def runGrossOffsets(self):

        ###Pieces of information needed
        
        ###These pieces of information come from the output of "topo" module from ISCE
        ### llh - size(3) - lat,lon,hgt of pixel under consideration
        ### los - size(2) - inc, azi LOS angles
        
        
        ###These pieces of information come from an external velocity product, e.g from NSIDC
        
        ###  vx - scalar  - Velocity in x direction at pixel under consideration
        ###  vy - scalar  - Velocity in y direction at pixel under consideration
        ###  vproj - string - Projection system of the velocity field
        ###        - EPSG:3031 for Antarctica
        ###        - EPSG:3413 for Greenland
        
        
        #### The equations below describe the operations needed for a single pixel
        #### I will use Greenland as an example. Easy to change for Antarctica by changing the coordinate system.
        
        ### Step 0: Set up projection transformers for ease of use
        self.llhProj = pyproj.Proj('+init=EPSG:4326')       ##Standard lat,lon, hgt
        self.xyzProj = pyproj.Proj('+init=EPSG:4978')    ##Standard xyz (ECEF)

        # From xy to lat lon.
        refPt = self.vProj(0.0, 0.0, inverse=True)

        print(refPt)
        
        ### Step 1: Set up radar image information
        azPixelSize = self.azPixelSize
        rngPixelSize = self.rngPixelSize
        
        ### Step 2: Cut the data
        version = 'v2' 
        print('Obtain the velocity data...')
        
        if version == 'v1':
            self.get_veloData()
        elif version == 'v2':
            self.get_veloData_v2()

        print('Extract the data to this radar scene...')

        # The following code is to be consistent with "get_offset_geometry" in dense_offset.py
        if self.mode == 'interior':

            numWinDown = (self.YSize - self.margin*2 - self.searchSizeHgt*2 - self.winSizeHgt) // self.skipSizeHgt
            numWinAcross = (self.XSize - self.margin*2 - self.searchSizeWidth*2 - self.winSizeWidth) // self.skipSizeWidth

            lat = np.zeros(shape=(numWinDown,numWinAcross)) 
            lon = np.zeros(shape=(numWinDown,numWinAcross))
            inc = np.zeros(shape=(numWinDown,numWinAcross))
            azi = np.zeros(shape=(numWinDown,numWinAcross))

            self.centerOffsetHgt = self.winSizeHgt//2-1
            self.centerOffsetWidth = self.winSizeWidth//2-1

        elif self.mode == 'exterior':

            numWinDown = self.numWinDown 
            numWinAcross = self.numWinAcross

            lat = self.lat
            lon = self.lon
            inc = self.inc
            azi = self.azi

        print("Number of winows in down direction, Number of window in across direction")
        print(numWinDown, numWinAcross)

        cut_vx = np.zeros(shape=(numWinDown,numWinAcross))
        cut_vy = np.zeros(shape=(numWinDown,numWinAcross))
        cut_v = np.zeros(shape=(numWinDown,numWinAcross))
        pixel = np.zeros(shape=(numWinDown,numWinAcross))
        line = np.zeros(shape=(numWinDown,numWinAcross))
 
        for iwin in range(numWinDown):

            # Need to calculate lat lon in the interior mode.
            if self.mode == 'interior':

                print('Processing line: ',iwin)
                
                down = self.margin + self.skipSizeHgt * iwin  + self.centerOffsetHgt
                
                off = down*self.XSize

                # Warning: depend on the ENVI format. This is for BSQ
                off2 = self.YSize * self.XSize + down*self.XSize
    
                across_indices = self.margin + np.arange(numWinAcross)*skipSizeWidth + self.centerOffsetWidth
    
                latline = np.memmap(filename=self.latfile,dtype='float64',offset=8*off,shape=(self.XSize))
                lonline = np.memmap(filename=self.lonfile,dtype='float64',offset=8*off,shape=(self.XSize))
                incline = np.memmap(filename=self.losfile,dtype='float32',offset=4*off,shape=(self.XSize))
                aziline = np.memmap(filename=self.losfile,dtype='float32',offset=4*off2,shape=(self.XSize))
    
                lat[iwin,:] = latline[across_indices]
                lon[iwin,:] = lonline[across_indices]
                inc[iwin,:] = incline[across_indices]
                azi[iwin,:] = aziline[across_indices]

            #print(iwin,': ',lon[iwin,:])
            #print(iwin,': ',lat[iwin,:])
            #print(iwin,': ',inc[iwin,:])
            #print(iwin,': ',azi[iwin,:])

            #### Look up in MEaSUREs InSAR-Based Antarctica Ice Velocity Map

            # Convert lat lon to grid coordinates in polar stereographic projection.
            xyMap = pyproj.transform(self.llhProj, self.vProj, lon[iwin,:], lat[iwin,:])
            #xyMap = self.vProj(lon[iwin,:],lat[iwin,:])
 
            # Extract the values in the velocity model.
            if version == 'v1':
                interval = 900
            elif version == 'v2':
                interval = 450

            pixel[iwin,:] = np.clip((xyMap[0]-self.x0[0])/interval, 0, self.vx.shape[1]-1)
            line[iwin,:] = np.clip((xyMap[1]-self.y0[0])/interval, 0, self.vx.shape[0]-1)

            pixel_int = pixel[iwin,:].astype(int)
            line_int = line[iwin,:].astype(int)

            # For Debug
            #print(iwin,': ', 'location: ', xyMap[0],xyMap[1]) 
            #print(iwin,': ', 'location: ', lon[iwin,:],lat[iwin,:])
            
            cut_vx[iwin,:] = self.vx[line_int,pixel_int]
            cut_vy[iwin,:] = self.vy[line_int,pixel_int]
            #cut_v[iwin,:] = np.sqrt(np.multiply(cut_vx[iwin,:],cut_vx[iwin,:]),np.multiply(cut_vy[iwin,:],cut_vy[iwin,:]))


        #print('cut_vx: ',cut_vx)
        #plt.figure()
        #plt.imshow(np.abs(self.vx),vmin=0, vmax=1000, cmap='coolwarm')
        #plt.show()
        #print('cut_vy: ',cut_vy)

        # Mask out invalid value based on the value of lat (or lon) (only work for polar region)
        # Mask out zero velocity

        # Removing bad data        
        #cut_vx = signal.medfilt(cut_vx, kernel_size = (5,5))
        #cut_vy = signal.medfilt(cut_vy, kernel_size = (5,5))

        #cut_vx = ndimage.median_filter(input=cut_vx,size=5,mode='nearest')
        #cut_vy = ndimage.median_filter(input=cut_vy,size=5,mode='nearest')

        cut_v = np.sqrt(np.multiply(cut_vx,cut_vx),np.multiply(cut_vy,cut_vy))
        
        valid = np.logical_and(inc!=0, cut_v!=0)


        ### Mask out invalid values ###
        # 1. Mask out invalid values at margin.
        cut_vx[inc==0] = np.nan
        cut_vy[inc==0] = np.nan
 
        # New logic: value 0 means unknown
        # 2. Mask out zero-speed values
        #cut_vx[cut_v==0] = np.nan
        #cut_vx[cut_v==0] = np.nan

        # 3. For version 1: There is a hole of biased values in velocity field, so mask it out here. (ad hoc)
        # Interpolation method
        # method 1: Linear interpolation
        # method 2: No interpolation
        if version == 'v1':
            
            hole_inds = np.logical_and (np.logical_and ( np.logical_and(lon> -75, lon< -74.5), lat <-77), lat > -77.2)
            cut_vx[hole_inds] = 0
            cut_vy[hole_inds] = 0

        if version == 'v1':
            interpolate_method = 2

        elif version == 'v2':
            interpolate_method = 2

        if interpolate_method == 1:

            print('Interpolating velocity field...')
            # Gerenate the new field for interpolation.

            # Linear interpolation.

            x0=np.arange(numWinAcross)
            y0=np.arange(numWinDown)
            xx,yy=np.meshgrid(x0,y0)
    
            grid_x,grid_y=[grid.ravel() for grid in np.meshgrid(x0,y0)]
            
            points = np.column_stack((xx[valid],yy[valid]))
            print(points.shape)
    
            # Interpolate VX.
            in_dat = cut_vx[valid]
            cut_vx_new = griddata(points, in_dat, (grid_x, grid_y), method='linear')
    
            # Interpolate VY.
            in_dat = cut_vy[valid]
            cut_vy_new = griddata(points, in_dat, (grid_x, grid_y), method='linear')
    
            # Get shaped.        
            cut_vx_new = cut_vx_new.reshape(numWinDown,numWinAcross)
            cut_vy_new = cut_vy_new.reshape(numWinDown,numWinAcross)

        elif interpolate_method == 2:
            # No interpolation.
            cut_vx_new = np.copy(cut_vx)
            cut_vy_new = np.copy(cut_vy)

       
        # Get Interpolated speed.
        cut_v_new = np.sqrt(np.multiply(cut_vx_new,cut_vx_new),np.multiply(cut_vy_new,cut_vy_new))

        print("The speed matrix") 
        print(cut_v_new)
        print("The shape of speed matrix")
        print(cut_v_new.shape)

        ##########
        #fig=plt.figure(10,figsize=(10,10))
        #ax = fig.add_subplot(111)
        #print(cut_v.shape)
        #ax.imshow(np.clip(cut_v_new,0,1000),cmap=plt.cm.viridis)
        #fig.savefig('10.png',format='png')

        ### Step 3: Convert XY velocity to EN velocity (clockwise rotation)

        print('Coverting XY to EN...')
        lonr = np.radians(lon - refPt[0])

        cut_ve = np.multiply(cut_vx_new, np.cos(lonr)) - np.multiply(cut_vy_new, np.sin(lonr))
        cut_vn = np.multiply(cut_vy_new, np.cos(lonr)) + np.multiply(cut_vx_new, np.sin(lonr))

        print('Polar stereographic velocity: ', [cut_vx_new, cut_vy_new])
        print('Local ENU velocity: ', [cut_ve, cut_vn])
 
        ####Step 4: Convert EN velocity to rng and azimuth
        #Local los and azi vector in ENU coordinate
        print(' Coverting EN to rdr...')
        incr = np.radians(inc)
        azir = np.radians(azi)
        losr = np.radians(azi-90.0)

        losenu=[ np.multiply(np.sin(incr),np.cos(losr)),
                 np.multiply(np.sin(incr),np.sin(losr)),
                 -np.cos(incr) ]
        
        azienu=[ np.cos(azir),
                 np.sin(azir),
                 0.0 ]

        # unit: pixel per day
        grossRangeOffset = (self.bTemp/365.25) * (cut_ve * losenu[0] + cut_vn * losenu[1])/ rngPixelSize
        grossAzimuthOffset = (self.bTemp/365.25) * (cut_ve * azienu[0] + cut_vn * azienu[1]) / azPixelSize

        # Add: 2019.02.12 Remove the invalid values at corner.
        # Mask out invalid values at margin.
        grossRangeOffset[inc==0] = np.nan
        grossAzimuthOffset[inc==0] = np.nan
 
        print('Gross azimuth offset: ', grossAzimuthOffset)
        print('Gross range offset: ', grossRangeOffset)

        ############ VISUALIZATION ##################
        ################### Float #############################
        fig=plt.figure(21,figsize=(16,9))
        #ax = fig.add_subplot(121)
        ax = fig.add_axes([0.05,0.05,0.4,0.9])
        ax.set_title('gross azimuth offset',fontsize=15)
        print(grossAzimuthOffset.shape)
        cax = ax.imshow(grossAzimuthOffset,cmap=plt.cm.coolwarm)
        cbar = fig.colorbar(cax,fraction=0.035,pad=0.04,ticks=np.arange(np.rint(np.nanmin(grossAzimuthOffset)),np.rint(np.nanmax(grossAzimuthOffset))+0.001))
        cbar.set_label("pixel",fontsize=15)

        #ax = fig.add_subplot(122)
        ax = fig.add_axes([0.55,0.05,0.4,0.9])
        ax.set_title('gross range offset',fontsize=15)
        print(grossRangeOffset.shape)
        cax = ax.imshow(grossRangeOffset,cmap=plt.cm.coolwarm)

        cbar = fig.colorbar(cax,fraction=0.035,pad=0.04,ticks=np.arange(np.rint(np.nanmin(grossRangeOffset)),np.rint(np.nanmax(grossRangeOffset))+0.001))
        cbar.set_label("pixel",fontsize=15)

        #fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        #fig.tight_layout()

        if self.figpath is not None:
            figname = os.path.join(self.figpath,'21.png')
        else:
            figname = '21.png'

        fig.savefig(figname,format='png')

        plt.close()

        # Save grossRangeOffset and grossAzimuthOffset as ISCE supported images.
        # Range
        if self.figpath is not None:
            rangeFileName = os.path.join(self.figpath, 'grossRange' + str(self.runid) +'.off')
        else:
            rangeFileName = os.path.join('grossRange.off')

        driver = gdal.GetDriverByName('ENVI')
        dst_ds = driver.Create(rangeFileName, xsize=grossRangeOffset.shape[1], ysize=grossRangeOffset.shape[0], bands=1, eType=gdal.GDT_Float32)
        dst_ds.GetRasterBand(1).WriteArray(grossRangeOffset,0,0)
        dst_ds = None

        outImage = isceobj.createImage()
        outImage.setDataType('FLOAT')
        outImage.setFilename(rangeFileName)
        outImage.setBands(1)
        outImage.scheme='BIL'
        outImage.setLength(grossRangeOffset.shape[0])
        outImage.setWidth(grossRangeOffset.shape[1])
        outImage.setAccessMode('read')
        outImage.renderHdr()

        # Azimuth
        if self.figpath is not None:
            azimuthFileName = os.path.join(self.figpath, 'grossAzimuth' + str(self.runid) + '.off')
        else:
            azimuthFileName = os.path.join('grossAzimuth.off')

        driver = gdal.GetDriverByName('ENVI')
        dst_ds = driver.Create(azimuthFileName, xsize=grossAzimuthOffset.shape[1], ysize=grossAzimuthOffset.shape[0], bands=1, eType=gdal.GDT_Float32)
        dst_ds.GetRasterBand(1).WriteArray(grossAzimuthOffset,0,0)
        dst_ds = None

        outImage = isceobj.createImage()
        outImage.setDataType('FLOAT')
        outImage.setFilename(azimuthFileName)
        outImage.setBands(1)
        outImage.scheme='BIL'
        outImage.setLength(grossAzimuthOffset.shape[0])
        outImage.setWidth(grossAzimuthOffset.shape[1])
        outImage.setAccessMode('read')
        outImage.renderHdr()


        ######### Round to integer ###########################
        fig=plt.figure(22,figsize=(16,9))
        #ax = fig.add_subplot(121)
        ax = fig.add_axes([0.05,0.05,0.4,0.9])
        ax.set_title('gross azimuth offset',fontsize=15)
        print(grossRangeOffset.shape)
        cax = ax.imshow(np.rint(grossAzimuthOffset),cmap=plt.cm.coolwarm)
        cbar = fig.colorbar(cax,fraction=0.035,pad=0.04,ticks=np.arange(np.rint(np.nanmin(grossAzimuthOffset)),np.rint(np.nanmax(grossAzimuthOffset))+0.001))
        cbar.set_label("pixel",fontsize=15)

        #ax = fig.add_subplot(122)
        ax = fig.add_axes([0.55,0.05,0.4,0.9])
        ax.set_title('gross range offset',fontsize=15)
        print(grossRangeOffset.shape)
        cax = ax.imshow(np.rint(grossRangeOffset),cmap=plt.cm.coolwarm)

        print("Before plotting the gross offsets (min and max): ", np.rint(np.nanmin(grossAzimuthOffset)),np.rint(np.nanmax(grossAzimuthOffset)))

        cbar = fig.colorbar(cax,fraction=0.035,pad=0.04,ticks=np.arange(np.rint(np.nanmin(grossRangeOffset)),np.rint(np.nanmax(grossRangeOffset))+0.001))
        cbar.set_label("pixel",fontsize=15)

        #fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        #fig.tight_layout()

        if self.figpath is not None:
            figname = os.path.join(self.figpath,'22.png')
        else:
            figname = '22.png'

        fig.savefig(figname,format='png')

        plt.close()

        return grossAzimuthOffset, grossRangeOffset

def main():

    grossObj = grossOffsets()
    grossObj.runGrossOffsets()

if __name__=='__main__':
    main()

