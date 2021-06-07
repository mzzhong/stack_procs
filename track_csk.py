#!/usr/bin/env python3

# Minyan Zhong # 

# read CSK xml files, and plot the tracks
# generate data xml files for roiApp.py

import os
import xml.etree.ElementTree as ET

import cartopy.crs as ccrs

import datetime
from datetime import date

import numpy as np

import pandas as pd 

from CSK_Utils import CSK_Utils

class track_csk():
    def __init__(self, track_num = None, date = None, satellite = None, direction = None):
        if track_num is not None:
            self.track_num = track_num
 
        if date:
            self.year = date.year
            self.month = date.month
            self.day = date.day
            self.date = date

        if satellite:
            self.satellite = satellite

        if direction:
            self.direction = direction
    
    def timestr1(self, t1):
        date, time=t1.split()
        year, month, day=date.split('-')
        hour, minute, second = time.split(':')
        second = int(round(float(second)))
        second = str(second)
        second = second.zfill(2)
        return year + month + day + hour + minute + second

    def timestr2(self,t1):
        dt = pd.to_datetime(t1)
        dt = dt.round('1s')
        return dt.strftime('%Y%m%d%H%M%S')

    def ariafilename(self,sate,direction,t1,t2):
       #CSKS3_RAW_B_HI_01_HH_RA_SF_20171119140722_20171119140729
        date, time=t1.split()
        year, month, day=date.split('-')
        P0 = '/'.join(["/net/kraken/nobak/mzzhong/CSKData/csk-raw_b",year,month,day]) + '/'
        P1 = sate
        P2 = 'RAW_B_HI_01_HH'
        if direction == 'asc':
            P3 = 'RA'
        else:
            P3 = 'RD'
        P4 = 'SF'
        P5 = self.timestr2(t1)
        P6 = self.timestr2(t2)
        name = '_'.join([P1,P2,P3,P4,P5,P6])
 
        h5filename = P0 + name + '/' + name + '.h5'

        return(h5filename)

    def show_track_footprint(self):

        # Information provided.
        # Date
        # Satellite
        # Direction.
        # These uniquely determines the track number.
        
        datestr = self.date.strftime("%Y%m%d")
        sate = self.satellite
    
        count = 0
        #xmlfileDir = "/net/jokull/nobak/mzzhong/CSKData/xmlOnly/88.86.187.38"
        xmlfileDir = "/net/kraken/bak/mzzhong/CSKData-FTP/xmlOnly/88.86.187.38"
        csk_footprint_dir = '/home/mzzhong/insarRoutines/track_footprints/csk-evans-2021'

        # create the folder
        os.makedirs(csk_footprint_dir, exist_ok=True)

        ## Find the track_coordinates.
        track_coordinates = []
        # sate is used
        for dirname in os.listdir(xmlfileDir + '/' + sate):
            for filename in os.listdir(xmlfileDir+'/'+ sate+'/'+dirname):
                # date is used
                if filename.endswith('.h5.xml') and datestr in filename:
                    print(filename)
                    tree = ET.parse(xmlfileDir+'/'+sate+'/'+dirname+'/'+filename)
                    root = tree.getroot()
                    lat=[]
                    lon=[]
                    for i in range(4):
                        #print(root[1][i+3].tag)
                        coord=root[1][i+3].text.split(' ')
                        lat.append(float(coord[0]))
                        lon.append(float(coord[1]))

                    # Find the start time the end of the scence
                    fmt="%Y-%m-%d %H:%M:%S"
                    start_time_str=root[1][1].text.split('.')[0]
                    start_time = datetime.datetime.strptime(start_time_str, fmt)
    
                    # Figure out this frame is ascending and descending.
                    if lon[1]>lon[0]:
                        current_direction='asc'
                    else:
                        current_direction='desc'
    
                    # Figure out which track the frame belongs to.
                    #   Date used
                    derived_track_num = CSK_Utils().date_sate_direc_2track(day=self.date,
                                                        sate = self.satellite,
                                                        direc = current_direction)
                      
                    # Find the start time and end time.
                    # data[1] start time; data[2] end time; data[9] satelillte number
                    # Now direction used.
                    if current_direction == self.direction:
                        # Track number also matches
                        assert derived_track_num==self.track_num, "track number doesn't pass check"

                        t1 = root[1][1].text
                        t2 = root[1][2].text
    
                        # Form the frame polygon.
                        lat[3],lat[2] = lat[2],lat[3]
                        lon[3],lon[2] = lon[2],lon[3]
                        lat.append(lat[0])
                        lon.append(lon[0])
    
                        # Save this frame.
                        track_coordinates.append((lon,lat,start_time))

        # Sort the frames according to start time
        track_coordinates = sorted(track_coordinates, key=lambda x:x[2])
        print(track_coordinates)

        # Work on this track, if it exists.
        if len(track_coordinates)>0:
            if not os.path.exists(csk_footprint_dir+'/'+'T'+str(self.track_num)):
                os.mkdir(csk_footprint_dir+'/'+'T'+str(self.track_num))

            # Plot the footprint/tracks.
            # Plot each side separately.
            asc_sides = []
            desc_sides = []

            for frame_coordinates in track_coordinates:
                lon, lat, start_time = frame_coordinates
                for ii in range(4):
                    
                    side = [lon[ii],lon[ii+1]], [lat[ii], lat[ii+1]]

                    k = (lat[ii]-lat[ii+1]) / (lon[ii]-lon[ii+1])
                    
                    if k >0:
                        desc_sides.append(side)
                    else:
                        asc_sides.append(side)

                # count the number of frames                    
                count = count + 1

            track_color = (218/255,165/255,32/255)

            if self.direction == 'asc':

                # The two ends.
                minlon = 1000
                minside = []

                maxlon = -1000
                maxside = []
                for side in desc_sides:
                    lon, lat = side
                    if min(lon) < minlon:
                        minlon = min(lon)
                        minside = side
                    if max(lon) > maxlon:
                        maxlon = max(lon)
                        maxside = side
                    
                asc_sides.append(maxside)
                asc_sides.append(minside)

                # Output
                for i, side in enumerate(asc_sides):
                    f = open(csk_footprint_dir+'/'+'T'+str(self.track_num)+'/'+'asc_track_'+str(i)+'.txt','w')
                    f.write(str(side[0][0]) + ' ' + str(side[1][0]) + '\n')
                    f.write(str(side[0][1]) + ' ' + str(side[1][1]) + '\n')
                    f.close()

            elif self.direction == 'desc':

                # The two ends.
                minlon = 1000
                minside = []

                maxlon = -1000
                maxside = []
                for side in asc_sides:
                    lon, lat = side
                    if min(lon) < minlon:
                        minlon = min(lon)
                        minside = side
                    if max(lon) > maxlon:
                        maxlon = max(lon)
                        maxside = side
                    
                desc_sides.append(maxside)
                desc_sides.append(minside)

                # Output
                for i, side in enumerate(desc_sides):
                    f = open(csk_footprint_dir+'/'+'T'+str(self.track_num)+'/'+'desc_track_'+str(i)+'.txt','w')
                    f.write(str(side[0][0]) + ' ' + str(side[1][0]) + '\n')
                    f.write(str(side[0][1]) + ' ' + str(side[1][1]) + '\n')
                    f.close()
 
 
            # Output the frames
            filefolder = csk_footprint_dir  + '/' + 'T' + str(self.track_num)

            if not os.path.exists(filefolder):
                os.mkdir(filefolder)

            for i, frame_coordinates in enumerate(track_coordinates):
                lon, lat, start_time = frame_coordinates
                filebasename = '_'.join(['frame_' + str(i) + '.txt'])
                filename = os.path.join(filefolder, filebasename)
                print(filename)
                f = open(filename,'w')
                for ii in range(5):
                    f.write(str(lon[ii])+' '+str(lat[ii])+'\n')
                f.close()

            ## Mark the time.
            ##ax.text(lon[0],lat[0],root[1][1].text.split(' ')[1].split('.')[0],transform=ccrs.PlateCarree())
            ##ax.text(lon[0],lat[0],dirname,transform=ccrs.PlateCarree())
            ##if count == 0:
            #    #plot.ax.text(lon[0],lat[0],daystr[4:8],transform=ccrs.PlateCarree())

            #if generate_xml: 
            #    h5filename = self.ariafilename(self.satellite,current_direction,t1,t2)
            #    h5filelist.write(h5filename + "\n")
            #    if count>1:
            #        h5filestr = h5filestr + ', '
            #    h5filestr = h5filestr + h5filename 
      
        print("total number of frames: ", count)

        return 0

    def generateXml(self,daystr, h5filestr):
        root = ET.Element('component', name='Sensor')
        ET.SubElement(root, 'property', name='HDF5').text = h5filestr
        ET.SubElement(root, 'property', name='OUTPUT').text = 'raw/'+daystr+'.raw'
        
        tree = ET.ElementTree(root)
        tree.write(daystr + '.xml')

def check_availablity(track_num, sate, date, dirc):

    #xmlfileDir = "/net/jokull/nobak/mzzhong/CSKData/xmlOnly/88.86.187.38"
    xmlfileDir = "/net/kraken/bak/mzzhong/CSKData-FTP/xmlOnly/88.86.187.38"

    # There are many frames of this satellite from different tracks.
    count_frame = 0
    for dirname in os.listdir(xmlfileDir + '/' + sate):
        if os.path.isdir(xmlfileDir+'/'+ sate +'/'+ dirname):
            # look through the files for this product
            for filename in os.listdir(xmlfileDir+'/'+ sate +'/'+ dirname):
                if filename.endswith('h5.xml') and date.strftime('%Y%m%d') in filename:
                    # check if the direction is correct
                    tree = ET.parse(xmlfileDir+'/'+sate+'/'+dirname+'/'+filename)
                    root = tree.getroot()
                    lat=[]
                    lon=[]
                    for i in range(4):
                        #print(root[1][i+3].tag)
                        coord=root[1][i+3].text.split(' ')
                        lat.append(float(coord[0]))
                        lon.append(float(coord[1]))

                    # Figure out this frame is ascending and descending.
                    if lon[1]>lon[0]:
                        current_direction='asc'
                    else:
                        current_direction='desc'

                    if current_direction == dirc:
                        count_frame += 1

    if count_frame == CSK_Utils().numOfFrames[track_num]:
        return True
    else:
        return False
 
def find_available_data(track_num, dirc):
    # search the first month data
    first_date = date(2017,11,16)
    last_date = date(2017,12,20)
 
    dates = CSK_Utils().track2date(track = track_num, first = first_date, last = last_date)

    #print(dates)

    for sate in sorted(dates.keys()):
        #print(sate)
        sate_dates = dates[sate]
        for today in sate_dates:
            #print(sate, today, check_availablity(sate, today, dirc))
            if check_availablity(track_num, sate, today, dirc):
                return (today, sate)

def main():
    track_numbers = range(0,22)
    for track_num in track_numbers:
        if track_num<11:
            dirc = "asc"
        else:
            dirc = "desc"

        date, sate = find_available_data(track_num, dirc)
        print(track_num, date, sate)
        
        track_csk_obj = track_csk(track_num = track_num, date=date, satellite = sate, direction=dirc)

        track_csk_obj.show_track_footprint()

if __name__ == '__main__':
    main()
