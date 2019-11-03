#!/usr/bin/env python3

import datetime
from datetime import date

import pickle

class CSK_Rutford_Utils():
    def __init__(self):

        # Initialization

        self.repo = "/net/kraken/nobak/mzzhong/CSK-Rutford/raw_data/aria-csk-dav.jpl.nasa.gov/repository/products/csk_rawb/v0.6"

        self.satelist = ['CSKS1','CSKS2','CSKS3','CSKS4']

        self.tracklist = [8,10,23,25,40,52,55,67,69,82,97,99,114,126,128,129,141,143,156,158,171,172,173,186,188,201,203,215,218,230,231,232]
        #print(len(self.tracklist))
 
        self.ref_day = {}
        # If CSKS2 is the reference
        self.ref_day['CSKS2']=0
        self.ref_day['CSKS3']=1
        self.ref_day['CSKS4']=4
        self.ref_day['CSKS1']=8

        #self.numOfFrames = [3,8,6,9,9,6,6,7,6,6,5,4,4,4,8,8,8,8,9,9,5,5]

        # From date to tracks (CSKS2)
        ref_track={}  
        ref_track[date(2013,8,19)] = [8, 10]
        ref_track[date(2013,8,20)] = [23, 25]
        ref_track[date(2013,8,21)] = [40]
        ref_track[date(2013,8,22)] = [52, 55]
        ref_track[date(2013,8,23)] = [67, 69]
        ref_track[date(2013,8,24)] = [82]
        ref_track[date(2013,8,9)] =  [97, 99]
        ref_track[date(2013,8,10)] = [114]
        ref_track[date(2013,8,11)] = [126, 128, 129]
        ref_track[date(2013,8,12)] = [141, 143]
        ref_track[date(2013,8,13)] = [156, 158]
        ref_track[date(2013,8,14)] = [171, 172, 173]
        ref_track[date(2013,8,15)] = [186, 188]
        ref_track[date(2013,8,16)] = [201, 203]
        ref_track[date(2013,8,17)] = [215, 218]
        ref_track[date(2013,8,18)] = [230, 231, 232]

        self.ref_track = ref_track

        # From track to date
        ref_date={}
        for track in self.tracklist:
            for day, crsp_tracks in ref_track.items():
                if track in crsp_tracks:
                    ref_date[track] = day

        self.ref_date = ref_date

        # New startdate and endate for CSK-RIS mission
        self.startdate = date(2013,1,1)
        self.stopdate = date(2014,12,31)

        self.maxframenum = 100

        # Load the acquistion time of each track
        aq_time_filename="/net/kraken/nobak/mzzhong/CSK-Rutford/raw_data/logistics/CSK_RIS_aq_time_of_tracks.pkl"
        with open(aq_time_filename,'rb') as f:
            self.aq_time = pickle.load(f)
        #print(self.aq_time)


    def date2track(self,day,sate=None,direction=None):
        
        tracks={}
        satelist = ['CSKS1','CSKS2','CSKS3','CSKS4']
        for satename in satelist:
            tracks[satename]=[]
            
            if sate != None and satename != sate:
                continue

            # For CSKS2 
            csks2date = day - datetime.timedelta(days=self.ref_day[satename])
            for ref_date,crsp_tracks in self.ref_track.items():
                
                if (csks2date - ref_date).days % 16 == 0 :
                    tracks[satename] = crsp_tracks

        return tracks

    def date_sate_direc_2track(self,day,sate,direc):
        
        tracks={}
        satelist = ['CSKS1','CSKS2','CSKS3','CSKS4']
        for satename in satelist:
            tracks[satename]=[]
            
            if sate != None and satename != sate:
                continue
            
            csks2date = day - datetime.timedelta(days=self.ref_day[satename])
            for ref_date,crsp_tracks in self.ref_track.items():
                
                if (csks2date - ref_date).days % 16 == 0 :
                    tracks[satename] = crsp_tracks

        track_num = None
        for this_track in tracks[sate]:
            if direc == 'asc' and this_track<=10:
                track_num = this_track
            if direc == 'desc' and this_track>=11:
                track_num = this_track
        
        return track_num


    def date_track2sate(self,day,track):
        
        satelist = ['CSKS1','CSKS2','CSKS3','CSKS4']
        for satename in satelist:

            csk2date = self.ref_date[track]
            cskdate = csk2date + datetime.timedelta(days=self.ref_day[satename])

            if (cskdate-day).days % 16 == 0:
                return satename
            
    def track2date(self,track,sate=None,first=None,last=None):
        if first == None:
            first = date(2017,11,16)
        if last == None:
            last = date(2020,11,16)

        satelist = ['CSKS1','CSKS2','CSKS3','CSKS4']

        dates={}
        for satename in satelist:
            dates[satename]=[]
            
            if sate != None and satename != sate:
                continue

            today = first 
            while today < last:
                
                if (today - (self.ref_date[track] + datetime.timedelta(days=self.ref_day[satename]))).days % 16 == 0:
                    break

                today = today + datetime.timedelta(days=1)

            while today < last:
                dates[satename].append(today)
                today = today + datetime.timedelta(days=16)

        return(dates)

if __name__=="__main__":
    CSK_Rutford_Utils()
    #tracks = CSK_Utils().date2track(day=date(2017,11,21),sate='CSKS3')
    
    #dates = CSK_Rutford_Utils().track2date(track=1,sate=None,first=date(2013,1,1),last=date(2016,1,1))
    #print(dates)
