#!/usr/bin/env python3
import numpy as np
import numpy.ma as ma

import xml.etree.ElementTree as ET
import argparse
import os
import sys

from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import gdal

import matplotlib.pyplot as plt

import pickle
import os

import pyproj


class Ant_data():

    def get_veloData_v0(self):
    
        # Velocity data.
        print("Getting velocity data...")

        velo_dir = '/net/jokull/nobak/mzzhong/Ant_Plot/Data/velocity_models'

        use_nc = 1
    
        if os.path.isfile(velo_dir + '/'+'AntVelo.npz') and use_nc == 0:
    
            npzfile = np.load(velo_dir + '/AntVelo.npz')
    
            v = npzfile['v']
            vel_lon = npzfile['vel_lon']
            vel_lat = npzfile['vel_lat']
        else:
            vel_file = velo_dir + '/antarctica_ice_velocity_900m.nc'
            fh=Dataset(vel_file,mode='r')
            vx = fh.variables['vx'][:]
            vy = fh.variables['vy'][:]
            v = np.sqrt(np.multiply(vx,vx)+np.multiply(vy,vy))
    
            x0 = np.arange(-2800000,2800000,step=900)+2800000
            y0 = np.arange(-2800000,2800000,step=900)+200+2800000
            x,y = np.meshgrid(x0,y0)
    
            AntVeloDataMap = Basemap(width=5600000,height=5600000,\
                                resolution='l',projection='stere',\
                                lat_ts=-71,lat_0=-90,lon_0=0)
    
            vel_lon, vel_lat= AntVeloDataMap(x,y,inverse="true")
            v = np.flipud(v)
    
            v = v/365

            #np.savez(velo_dir + '/AntVelo.npz',v=v,vel_lon=vel_lon,vel_lat=vel_lat)
    
        return v, vel_lon,vel_lat
 

    def get_veloData_v1_wrong(self):

        # Velocity data.
        print("Getting velocity data...")

        velo_dir = '/net/jokull/nobak/mzzhong/Ant_Plot/Data/velocity_models'
        file_basename = 'antarctica_ice_velocity_900m.nc'
        vel_file = os.path.join(velo_dir, file_basename)


        npz_filebasename = 'AntVelo_v1.npz'
        redo = 0
        if os.path.isfile(velo_dir + '/' + npz_filebasename) and redo == 0:
    
            npzfile = np.load(velo_dir + '/' + npz_filebasename)
    
            self.vel_lon = npzfile['vel_lon']
            self.vel_lat = npzfile['vel_lat']
            self.ve = npzfile['ve']
            self.vn = npzfile['vn']
            self.v_comb = npzfile['v_comb']

            self.vy = npzfile['vx']
            self.vx = npzfile['vy']
            self.v = npzfile['v']



        else:
            fh=Dataset(vel_file,mode='r')
    
            for key in fh.variables.keys():
                print(fh.variables[key])
    
            x0 = np.arange(-2800000,2800000,step=900)
            y0 = np.arange(-2800000,2800000,step=900)+200
            x,y = np.meshgrid(x0,y0)
            vProj = pyproj.Proj('+init=EPSG:3031')
            llhProj = pyproj.Proj('+init=EPSG:4326')
            vel_lon, vel_lat = pyproj.transform(vProj, llhProj, x,y)
    
            # Read in vx, vy.
            vx = fh.variables['vx'][:]
            vy = fh.variables['vy'][:]

            # Per day.
            vx = vx/365
            vy = vy/365

            # Angle.
            # Find the projection after one day.
            x_1 = x + vx
            y_1 = y + vy
            vel_lon_1, vel_lat_1 = pyproj.transform(vProj, llhProj, x_1, y_1)

            vel_lon_shift = vel_lon_1 - vel_lon
            vel_lat_shift = vel_lat_1 - vel_lat

            # Radius.
            radii = 6371*1000
            # East component. Rad times radius. 
            ve = np.multiply(np.deg2rad(vel_lon_shift), radii * np.cos(np.deg2rad(np.abs(vel_lat))))
            # North component.
            vn = np.deg2rad(vel_lat_shift) * radii

            # Speed.
            v = np.sqrt(np.multiply(vx,vx)+np.multiply(vy,vy))

            # Speed for checking.
            v_comb = np.sqrt(np.multiply(ve,ve)+np.multiply(vn,vn))
            
            # Flip all velocity matrix before saving.
            vx = np.flipud(vx)
            vy = np.flipud(vy)
            v = np.flipud(v)

            vn = np.flipud(vn)
            ve = np.flipud(ve)
            v_comb = np.flipud(v_comb)            


            # Save to the class.
            self.vel_lon = vel_lon
            self.vel_lat = vel_lat

            self.vx = vx
            self.vy = vy
            self.vn = vn
            self.ve = ve

            self.v = v
            self.v_comb = v_comb
           
            # Clip the model.
            # Must be kept.
            #ind_lon = np.logical_and(vel_lon > -90, vel_lon < -60)
            #ind_lat = np.logical_and(vel_lat > -79, vel_lat < -72)
            #
            #ind_kept = np.logical_or(ind_lon,ind_lat)
    
            #left = np.where(np.sum(ind_kept,axis=0)>0)[0][0]
            #right = np.where(np.sum(ind_kept,axis=0)>0)[0][-1]
    
            #up = np.where(np.sum(ind_kept,axis=1)>0)[0][0]
            #down = np.where(np.sum(ind_kept,axis=1)>0)[0][-1]
    
            #print(up,down,left,right)
    
            #v = v[up : down,  left:right]
            #vel_lon = vel_lon[up:down, left:right]
            #vel_lat = vel_lat[up:down, left:right]

            np.savez(velo_dir + '/'+ npz_filebasename, v_comb=v_comb, ve=ve, vn=vn, vel_lon=vel_lon, vel_lat=vel_lat, v=v, vx=vx, vy=vy)


    def get_veloData_v1(self):

        # Velocity data.
        print("Getting velocity data...")

        velo_dir = '/net/jokull/nobak/mzzhong/Ant_Plot/Data/velocity_models'
        file_basename = 'antarctica_ice_velocity_900m.nc'
        vel_file = os.path.join(velo_dir, file_basename)


        npz_filebasename = 'AntVelo_v1.npz'
        redo = 0
        if os.path.isfile(velo_dir + '/' + npz_filebasename) and redo == 0:
    
            npzfile = np.load(velo_dir + '/' + npz_filebasename)
    
            self.vel_lon = npzfile['vel_lon']
            self.vel_lat = npzfile['vel_lat']
            self.ve = npzfile['ve']
            self.vn = npzfile['vn']
            self.v_comb = npzfile['v_comb']

            self.vy = npzfile['vx']
            self.vx = npzfile['vy']
            self.v = npzfile['v']

        else:
            fh=Dataset(vel_file,mode='r')
    
            for key in fh.variables.keys():
                print(fh.variables[key])
    
            x0 = np.arange(-2800000,2800000,step=900)
            y0 = np.arange(-2800000,2800000,step=900)+200
            x,y = np.meshgrid(x0,y0)
            vProj = pyproj.Proj('+init=EPSG:3031')
            llhProj = pyproj.Proj('+init=EPSG:4326')
            vel_lon, vel_lat = pyproj.transform(vProj, llhProj, x,y)
    
            # Read in vx, vy.
            vx = fh.variables['vx'][:]
            vy = fh.variables['vy'][:]
            vx = np.flipud(vx)
            vy = np.flipud(vy)

            # Velocity meter per day.
            vx = vx/365
            vy = vy/365

            # Convert from XY to EN
            # refPt.
            refPt = vProj(0.0, 0.0, inverse=True)
            # refPt: lon 0, lat -90

            print('Coverting XY to EN...')
            lonr = np.radians(vel_lon - refPt[0])

            ve = np.multiply(vx, np.cos(lonr)) - np.multiply(vy, np.sin(lonr))
            vn = np.multiply(vy, np.cos(lonr)) + np.multiply(vx, np.sin(lonr))

            # Speed.
            v = np.sqrt(np.multiply(vx,vx)+np.multiply(vy,vy))

            # Speed for checking.
            v_comb = np.sqrt(np.multiply(ve,ve)+np.multiply(vn,vn))
            
            # Save to the class.
            self.vel_lon = vel_lon
            self.vel_lat = vel_lat

            self.vx = vx
            self.vy = vy
            self.vn = vn
            self.ve = ve

            self.v = v
            self.v_comb = v_comb
           
            # Clip the model.
            # Must be kept.
            #ind_lon = np.logical_and(vel_lon > -90, vel_lon < -60)
            #ind_lat = np.logical_and(vel_lat > -79, vel_lat < -72)
            #
            #ind_kept = np.logical_or(ind_lon,ind_lat)
    
            #left = np.where(np.sum(ind_kept,axis=0)>0)[0][0]
            #right = np.where(np.sum(ind_kept,axis=0)>0)[0][-1]
    
            #up = np.where(np.sum(ind_kept,axis=1)>0)[0][0]
            #down = np.where(np.sum(ind_kept,axis=1)>0)[0][-1]
    
            #print(up,down,left,right)
    
            #v = v[up : down,  left:right]
            #vel_lon = vel_lon[up:down, left:right]
            #vel_lat = vel_lat[up:down, left:right]

            np.savez(velo_dir + '/'+ npz_filebasename, v_comb=v_comb, ve=ve, vn=vn, vel_lon=vel_lon, vel_lat=vel_lat, v=v, vx=vx, vy=vy)
 
    def get_veloData_v2(self):
    
        # Velocity data.
        print("Getting velocity data...")

        velo_dir = '/net/jokull/nobak/mzzhong/Ant_Plot/Data/velocity_models'
        file_basename = 'Antarctica_ice_velocity_2016_2017_1km_v01.nc'
        vel_file = os.path.join(velo_dir, file_basename)

        fh=Dataset(vel_file,mode='r')

        # Available keys:
        # 'coord_system', 'x', 'y', 'lat', 'lon', 'VX', 'VY', 
        # 'STDX', 'STDY', 'ERRX', 'ERRY', 'CNT'

        for key in fh.variables.keys():
            print(fh.variables[key])

        # Coordinates transformation.
        vel_lon = fh.variables['lon'][:]
        vel_lat = fh.variables['lat'][:]


        #x0 = np.arange(-2800000,2800000,step=900)
        #y0 = np.arange(-2800000,2800000,step=900)+200
        #x,y = np.meshgrid(x0,y0)
        #vProj = pyproj.Proj('+init=EPSG:3031')
        #llhProj = pyproj.Proj('+init=EPSG:4326')
        #vel_lon, vel_lat = pyproj.transform(vProj, llhProj, x,y)

        # Velocity.
        vx = fh.variables['VX'][:]
        vy = fh.variables['VY'][:]
        v = np.sqrt(np.multiply(vx,vx)+np.multiply(vy,vy))
        v = np.flipud(v)
        v = v/365

        fig = plt.figure(1,figsize=(10,10))
        ax = fig.add_subplot(111)
        im = ax.imshow(v)
        fig.colorbar(im)
        plt.show()


        print(v)
        print(stop)

        return v, vel_lon,vel_lat


    def get_glData(self):
        
        # Groundline data.
        print("Getting grounding line data...")
    
        #if os.path.isfile('/net/jokull/nobak/mzzhong/Ant_Plot/Data/GL_Points.npz'):
        #    npzfile = np.load('/net/jokull/nobak/mzzhong/Ant_Plot/Data/GL_Points.npz')
        #    gldata = npzfile['gldata']
        #else:
        #    gl_file = '/net/jokull/nobak/mzzhong/Ant_Plot/Data/GL_Points.txt'
        #    gldata = np.genfromtxt(gl_file,dtype=None)
        #    gldata[gldata[:,0]>180,0] = gldata[gldata[:,0]>180,0]-360
    
        #    np.savez('/net/jokull/nobak/mzzhong/Ant_Plot/Data/GL_Points.npz',gldata=gldata)

        gl_file = '/net/jokull/nobak/mzzhong/Ant_Plot/Data/GL_Points_Evans_min.txt'
        gldata = np.genfromtxt(gl_file,dtype=None)
        gldata[gldata[:,0]>180,0] = gldata[gldata[:,0]>180,0]-360
    
        return gldata


    def get_CryoSat_DEM(self):

        # Velocity data.
        print("Getting cryosat dem data...")

        dem_dir = '/net/jokull/nobak/mzzhong/Ant_Plot/Data/DEM-CryoSat'
        file_basename = 'Antarctica_Cryosat2_1km_DEMv1.0.nc'
        dem_file = os.path.join(dem_dir, file_basename)


        fh=Dataset(dem_file,mode='r')

        #for key in fh.variables.keys():
        #    print(fh.variables[key])

        #x0 = np.arange(-2819500,2820500,step=1000)
        #y0 = np.arange(-2419500,2420500,step=1000)

        print(fh.variables)
        x0 = fh.variables['x'][:]
        y0 = fh.variables['y'][:]
        
        dem = fh.variables['z'][:]

        x,y = np.meshgrid(x0,y0)
        vProj = pyproj.Proj('+init=EPSG:3031')
        llhProj = pyproj.Proj('+init=EPSG:4326')
        dem_lon, dem_lat = pyproj.transform(vProj, llhProj, x,y)


        self.dem_lon = dem_lon
        self.dem_lat = dem_lat
        self.dem = dem


        write_to_xyz = True

        if write_to_xyz:
            r_dem_lon = np.round(dem_lon * 50)/50
            r_dem_lat = np.round(dem_lat * 200)/200
            f = open(os.path.join(dem_dir,'CryoSat_DEM.xyz'),'w')
            keys = []

            for row in range(dem_lon.shape[0]):
                for col in range(dem_lon.shape[1]):
                    lon = r_dem_lon[row, col]
                    lat = r_dem_lat[row, col]
                    z = dem[row,col]

                    if lon>-100 and lon<-55 and lat>-85 and lat<-70 and not np.isnan(z):
                        f.write(str(lon)+' '+str(lat)+' '+str(z)+'\n')
                        keys.append((lon,lat))

            f.close()
        
        return 0

if __name__=='__main__':
    ant = Ant_data() 
    ant.get_CryoSat_DEM()
