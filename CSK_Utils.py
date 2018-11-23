#!/usr/bin/env python3

import datetime
from datetime import date

class CSK_Utils():
    def __init__(self):

        self.satelist = ['CSKS1','CSKS2','CSKS3','CSKS4']
 
        self.ref_day={}
        self.ref_day['CSKS2']=0
        self.ref_day['CSKS3']=1
        self.ref_day['CSKS4']=4
        self.ref_day['CSKS1']=8


        self.numOfFrames = [3,8,6,9,9,6,6,7,6,6,5,4,4,4,8,8,8,8,9,9,5,5]

        ref_track={}  
        ref_track[date(2018,3,3)] = [9,20]
        ref_track[date(2018,3,4)] = [6]
        ref_track[date(2018,3,5)] = [3]
        ref_track[date(2018,3,6)] = [0,13]
        ref_track[date(2018,3,7)] = [16]
        ref_track[date(2018,3,8)] = [10,19]
        ref_track[date(2018,3,9)] = [7]
        ref_track[date(2018,3,10)] = [4]
        ref_track[date(2018,3,11)] = [1,12]
        ref_track[date(2018,3,12)] = [15]
        ref_track[date(2018,3,13)] = [18]
        ref_track[date(2018,3,14)] = [8,21]
        ref_track[date(2018,3,15)] = [5]
        ref_track[date(2018,3,16)] = [2,11]
        ref_track[date(2018,3,17)] = [14]
        ref_track[date(2018,3,18)] = [17]

        self.ref_track = ref_track

        ref_date={}
        for track in range(22):
            for day, crsp_tracks in ref_track.items():
                if track in crsp_tracks:
                    ref_date[track] = day


        self.ref_date = ref_date

        self.startdate = date(2017,11,16)
        self.enddate = date(2020,11,16)

        self.maxframenum = 100

    def date2track(self,day,sate=None,direction=None):
        
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

        return tracks

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
    CSK_Utils()
    #tracks = CSK_Utils().date2track(day=date(2017,11,21),sate='CSKS3')
    
    dates = CSK_Utils().track2date(track=1,sate=None,first=date(2017,11,16),last=date(2018,12,20))
    print(dates)
