#!/usr/bin/env python3

# Read, process, and convert Antarctica data.

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

        velo_dir = '/net/jokull/nobak/mzzhong/Ant_Data/velocity_models'

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
 
    def get_veloData_v1(self):

        # Velocity data.
        print("Getting velocity data...")

        velo_dir = '/net/jokull/nobak/mzzhong/Ant_Data/velocity_models'
        file_basename = 'antarctica_ice_velocity_900m.nc'
        vel_file = os.path.join(velo_dir, file_basename)


        npz_filebasename = 'AntVelo_v1.npz'
        npz_filebasename_error = 'AntVelo_v1_error.npz'

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

            #####################
            npzfile = np.load(velo_dir + '/' + npz_filebasename_error)
    
            self.ve_err = npzfile['ve_err']
            self.vn_err = npzfile['vn_err']
            self.err = npzfile['err']
            self.err_comb = npzfile['err_comb']

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
            err = fh.variables['err'][:]

            vx = np.flipud(vx)
            vy = np.flipud(vy)
            err = np.flipud(err)

            # Velocity meter per day.
            vx = vx/365
            vy = vy/365
            err = err/365

            vx_err = np.multiply(err, np.divide(vx, np.sqrt(vx**2 + vy**2)+0.0001))
            vy_err = np.multiply(err, np.divide(vy, np.sqrt(vx**2 + vy**2)+0.0001))

            
            # Convert from XY to EN
            # refPt.
            refPt = vProj(0.0, 0.0, inverse=True)
            # refPt: lon 0, lat -90

            print('Coverting XY to EN...')
            lonr = np.radians(vel_lon - refPt[0])

            ve = np.multiply(vx, np.cos(lonr)) - np.multiply(vy, np.sin(lonr))
            vn = np.multiply(vy, np.cos(lonr)) + np.multiply(vx, np.sin(lonr))

            ve_err = np.abs(np.multiply(vx_err, np.cos(lonr)) - np.multiply(vy_err, np.sin(lonr)))
            vn_err = np.abs(np.multiply(vy_err, np.cos(lonr)) + np.multiply(vx_err, np.sin(lonr)))

            # Speed.
            v = np.sqrt(np.multiply(vx,vx)+np.multiply(vy,vy))

            # Speed for checking.
            v_comb = np.sqrt(np.multiply(ve,ve)+np.multiply(vn,vn))

            # err for checking.
            err_comb = np.sqrt(np.multiply(ve_err,ve_err)+np.multiply(vn_err,vn_err))
 
            
            # Save to the class.
            self.vel_lon = vel_lon
            self.vel_lat = vel_lat

            self.vx = vx
            self.vy = vy

            self.ve = ve
            self.vn = vn

            self.v = v
            self.v_comb = v_comb


            self.vx_err = vx_err
            self.vy_err = vy_err
            self.ve_err = ve_err
            self.vn_err = vn_err
            self.err = err
            self.err_comb = err_comb

            # Save the results. 
            np.savez(velo_dir + '/'+ npz_filebasename, v_comb=v_comb, ve=ve, vn=vn, vel_lon=vel_lon, vel_lat=vel_lat, v=v, vx=vx, vy=vy)
            np.savez(velo_dir + '/'+ npz_filebasename_error, err_comb=err_comb, ve_err=ve_err, vn_err = vn_err, err=err)


        write_to_xyz = True
        if write_to_xyz:
            print('Writing...')

            vel_lon = np.round(self.vel_lon * 50)/50
            vel_lat = np.round(self.vel_lat * 200)/200
            vx = self.vx
            vy = self.vy
            v = self.v

            f = open(os.path.join(velo_dir,'FRIS_v1.xyz'),'w')
            keys = []

            for row in range(vel_lon.shape[0]):
                for col in range(vel_lon.shape[1]):

                    lon = vel_lon[row, col]
                    if lon > 180:
                        lon = lon-360

                    lat = vel_lat[row, col]
                    
                    if lon>-95 and lon<-30 and lat>-85 and lat<-72:
                        #vxx = vx[row, col]
                        #vyy = vy[row, col]
                        #z = np.sqrt(vxx**2 + vyy**2)

                        z = v[row, col]

                        if z > 0:
                            #print(str(lon)+' '+str(lat)+' '+str(z)+'\n')
                            f.write(str(lon)+' '+str(lat)+' '+str(z)+'\n')
                            #keys.append((lon,lat))

            f.close()
 

    def get_veloData_v2(self):

        # Available keys:
        # 'coord_system', 'x', 'y', 'lat', 'lon', 'VX', 'VY', 
        # 'STDX', 'STDY', 'ERRX', 'ERRY', 'CNT'

        # Velocity data.
        print("Getting velocity data...")

        velo_dir = '/net/jokull/nobak/mzzhong/Ant_Data/velocity_models'
        file_basename = 'another/antarctica_ice_velocity_450m_v2.nc'
        vel_file = os.path.join(velo_dir, file_basename)

        npz_filebasename = 'AntVelo_v2.npz'
        #npz_filebasename_error = 'AntVelo_v2_error.npz'

        redo = 0
        if os.path.isfile(velo_dir + '/' + npz_filebasename) and redo == 0:

            print('Loading v2...') 
            npzfile = np.load(velo_dir + '/' + npz_filebasename)
    
            self.vel_lon = npzfile['vel_lon']
            self.vel_lat = npzfile['vel_lat']
            self.ve = npzfile['ve']
            self.vn = npzfile['vn']
            self.v_comb = npzfile['v_comb']

            self.vy = npzfile['vx']
            self.vx = npzfile['vy']
            self.v = npzfile['v']

            #####################
            #npzfile = np.load(velo_dir + '/' + npz_filebasename_error)
    
            #self.ve_err = npzfile['ve_err']
            #self.vn_err = npzfile['vn_err']
            #self.err = npzfile['err']
            #self.err_comb = npzfile['err_comb']


        else:
            fh=Dataset(vel_file,mode='r')
    
            for key in fh.variables.keys():
                print(fh.variables[key])

            # Coordinates transformation.
            vel_lon = fh.variables['lon'][:]
            vel_lat = fh.variables['lat'][:]
    
            vx = fh.variables['VX'][:].data
            vy = fh.variables['VY'][:].data

            print(vel_lon)

            plt.figure(1, figsize=(10,10))
            plt.imshow(np.abs(vel_lon), vmin=0, vmax=1000, cmap='coolwarm')
            plt.savefig('1.png')

            vx = vx/365
            vy = vy/365

            # Downsample:
            #vel_lon = vel_lon[0::2,0::2]
            #vel_lat = vel_lat[0::2,0::2]
            #vx = vx[0::2,0::2]
            #vy = vy[0::2,0::2] 

            # Convert from XY to EN
            # refPt.
            vProj = pyproj.Proj('+init=EPSG:3031')
            refPt = vProj(0.0, 0.0, inverse=True)
            # refPt: lon 0, lat -90

            print('Converting XY to EN...')
            lonr = np.radians(vel_lon - refPt[0])

            ve = np.multiply(vx, np.cos(lonr)) - np.multiply(vy, np.sin(lonr))
            vn = np.multiply(vy, np.cos(lonr)) + np.multiply(vx, np.sin(lonr))

            # Speed.
            v = np.sqrt(vx**2 + vy**2)

            # Speed for checking.
            v_comb = np.sqrt(np.multiply(ve,ve)+np.multiply(vn,vn))

            # Save to the class.
            self.vel_lon = vel_lon
            self.vel_lat = vel_lat

            self.vx = vx
            self.vy = vy

            self.ve = ve
            self.vn = vn

            self.v = v
            self.v_comb = v_comb

            print('Saving...')

            # Save the results. 
            np.savez(velo_dir + '/'+ npz_filebasename, v_comb=v_comb, ve=ve, vn=vn, vel_lon=vel_lon, vel_lat=vel_lat, v=v, vx=vx, vy=vy)
            
            #np.savez(velo_dir + '/'+ npz_filebasename_error, err_comb=err_comb, ve_err=ve_err, vn_err = vn_err, err=err)

        write_to_xyz = True

        if write_to_xyz:
            print('Writing...')

            vel_lon = np.round(self.vel_lon * 50)/50
            vel_lat = np.round(self.vel_lat * 200)/200
            vx = self.vx
            vy = self.vy
            v = self.v

            f = open(os.path.join(velo_dir,'Evans_v2.xyz'),'w')
            keys = []

            for row in range(vel_lon.shape[0]):
                for col in range(vel_lon.shape[1]):

                    lon = vel_lon[row, col]
                    if lon > 180:
                        lon = lon-360

                    lat = vel_lat[row, col]
                    
                    if lon>-95 and lon<-60 and lat>-82 and lat<-72:
                        #vxx = vx[row, col]
                        #vyy = vy[row, col]
                        #z = np.sqrt(vxx**2 + vyy**2)

                        z = v[row, col]

                        # Only non-zero values are output.
                        if z > 0:
                            #print(str(lon)+' '+str(lat)+' '+str(z)+'\n')
                            f.write(str(lon)+' '+str(lat)+' '+str(z)+'\n')
                        #keys.append((lon,lat))

            f.close()
        
        return 0

    def get_glData(self):
        
        # Groundline data.
        print("Getting grounding line data...")
    
        #if os.path.isfile('/net/jokull/nobak/mzzhong/Ant_Data/GL_Points.npz'):
        #    npzfile = np.load('/net/jokull/nobak/mzzhong/Ant_Data/GL_Points.npz')
        #    gldata = npzfile['gldata']
        #else:
        #    gl_file = '/net/jokull/nobak/mzzhong/Ant_Data/GL_Points.txt'
        #    gldata = np.genfromtxt(gl_file,dtype=None)
        #    gldata[gldata[:,0]>180,0] = gldata[gldata[:,0]>180,0]-360
    
        #    np.savez('/net/jokull/nobak/mzzhong/Ant_Data/GL_Points.npz',gldata=gldata)

        gl_file = '/net/jokull/nobak/mzzhong/Ant_Data/GL/GL_Points_Evans_min.txt'
        gldata = np.genfromtxt(gl_file,dtype=None)
        gldata[gldata[:,0]>180,0] = gldata[gldata[:,0]>180,0]-360
    
        return gldata


    def get_CryoSat_DEM(self):

        # Velocity data.
        print("Getting cryosat dem data...")

        dem_dir = '/net/jokull/nobak/mzzhong/Ant_Data/DEM-CryoSat'
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
            f = open(os.path.join(dem_dir,'FRIS_DEM.xyz'),'w')
            keys = []

            for row in range(dem_lon.shape[0]):
                for col in range(dem_lon.shape[1]):
                    lon = r_dem_lon[row, col]

                    lat = r_dem_lat[row, col]
                    z = dem[row,col]

                    if lon>-95 and lon<-30 and lat>-85 and lat<-72 and not np.isnan(z):
                        f.write(str(lon)+' '+str(lat)+' '+str(z)+'\n')
                        keys.append((lon,lat))

            f.close()
        
        return 0

if __name__=='__main__':
    ant = Ant_data()
    #ant.get_veloData_v1()
    ant.get_veloData_v2()
    #ant.get_CryoSat_DEM()
