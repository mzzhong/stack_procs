#!/usr/bin/env python3

import datetime
from datetime import date

class S1_Utils():
    def __init__(self):

        tracks = [37, 52, 169, 65]
        self.numberOfFrames = {}
        self.numberOfFrames[37] = 5
        self.numberOfFrames[52] = 5
        self.numberOfFrames[169] = 4
        self.numberOfFrames[65] = 4

        # from date to track
        ref_track={}  
        ref_track[date(2018,2,13)] = [37]
        ref_track[date(2018,2,20)] = [52]
        ref_track[date(2018,2,22)] = [169]
        ref_track[date(2018,2,21)] = [65]

        self.ref_track = ref_track

        # from track to date
        ref_date={}
        for track in tracks:
            for day, crsp_tracks in ref_track.items():
                if track in crsp_tracks:
                    ref_date[track] = day
        self.ref_date = ref_date

        self.startdate = date(2014,5,1)
        self.enddate = date(2024,5,1)

    def track2date(self,track,first=None,last=None):
        
        if first == None:
            first = self.startdate
        if last == None:
            last = self.enddate

        today = first 
        while today < last:
            if (today - self.ref_date[track]).days % 6 == 0:
                break
            today = today + datetime.timedelta(days=1)

        dates = [today]

        today = today + datetime.timedelta(days=6)

        while today <= last:
            dates.append(today)
            today = today + datetime.timedelta(days=6)
 
        return(dates)

if __name__=="__main__":
    S1_Utils()
