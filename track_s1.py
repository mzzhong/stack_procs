#!/usr/bin/env python3

# read S1 zip file and plot the tracks

import os
import xml.etree.ElementTree as ET

#import cartopy.crs as ccrs

import time
import datetime
from datetime import date

def get_lat_lon_v2(safefile):

    from xml.etree import ElementTree as ET
    import zipfile
    lats = []
    lons = []
    
    zf = zipfile.ZipFile(safefile,'r')
    fname = os.path.join(os.path.basename(safefile).replace('zip','SAFE'), 'preview/map-overlay.kml')
    xmlstr = zf.read(fname)
    xmlstr=xmlstr.decode('utf-8')
    start = '<coordinates>'
    end = '</coordinates>'
    pnts = xmlstr[xmlstr.find(start)+len(start):xmlstr.find(end)].split()
    
    for pnt in pnts:
        lons.append(float(pnt.split(',')[0]))
        lats.append(float(pnt.split(',')[1]))

    lons.append(lons[0])
    lats.append(lats[0])
    return lons, lats

class track_s1():
    def __init__(self,track_num = None):
        self.path = '/net/kraken/nobak/mzzhong/S1-Evans/footprints' + '/' + 'track_' + str(track_num)
        self.color = 'k'
        self.track_num = track_num

    def showTracks(self,plot=None):

        count = 0

        for product in os.listdir(self.path):
            if product.endswith('.zip'):
                print(product)
                safefile = self.path + '/' + product
                lon,lat = get_lat_lon_v2(safefile)
                print(lon)
                print(lat)
                #plot.ax.plot(lon,lat, self.color, transform=ccrs.PlateCarree())

                # Output the footprint                
                count = count + 1
                f = open('/home/mzzhong/insarRoutines/track_footprints/s1/'+str(self.track_num)+'_'+str(count)+'.txt','w')
                for ii in range(len(lon)):
                     f.write(str(lon[ii]) + ' ' + str(lat[ii]) + '\n')
                f.close()



#                # time
#                #ax.text(lon[0],lat[0],root[1][1].text.split(' ')[1].split('.')[0],transform=ccrs.PlateCarree())
#                #ax.text(lon[0],lat[0],dirname,transform=ccrs.PlateCarree())
#                count = count + 1

def main():

    mytrack = track_s1(64)
    mytrack.showTracks()

    return 0

if __name__=="__main__":
    main()
