#!/usr/bin/env python3
import datetime
from datetime import date

class CSK_Utils():
    def __init__(self):

        self.satelist = ['CSKS1','CSKS2','CSKS3','CSKS4']

        # number of frames
        self.numOfFrames = [3,8,6,9,9,6,6,7,6,6,5,4,4,4,8,8,8,8,9,9,5,5]

        # reduced_track_frames
        # This is for ordering data focused on the Evans Ice Shelf
        # order time: 2020 and 2021
        # 2020 order: only descending
        # 2021 order: ascending 1,2,3,4 and fix descending 11 and 12
        self.reduced_track_frames = {}
        self.reduced_track_frames[0] = []
        self.reduced_track_frames[1] = [0,1,2,3]
        self.reduced_track_frames[2] = [0,1,2,3,4]
        self.reduced_track_frames[3] = [0,1,2,3,4,5]
        self.reduced_track_frames[4] = [0,1,2,3,4,5]

        # 2020 version
        #self.reduced_track_frames[5] = [1,2,3]
        # 2021 version
        self.reduced_track_frames[5] = [1,2,3]

        self.reduced_track_frames[6] = []
        self.reduced_track_frames[7] = []
        self.reduced_track_frames[8] = []
        self.reduced_track_frames[9] = []
        self.reduced_track_frames[10] = []
        self.reduced_track_frames[11] = [0,1,2,3]
        self.reduced_track_frames[12] = [0,1,2,3]
        self.reduced_track_frames[13] = [0,1,2,3]
        self.reduced_track_frames[14] = [5,6,7]
        self.reduced_track_frames[15] = [5,6,7]
        self.reduced_track_frames[16] = [4,5,6,7]
        self.reduced_track_frames[17] = [3,4,5,6,7]
        self.reduced_track_frames[18] = [3,4,5,6]
        self.reduced_track_frames[19] = []
        self.reduced_track_frames[20] = []
        self.reduced_track_frames[21] = []

        # ref_day
        self.ref_day={}
        self.ref_day['stage1'] = {}
        self.ref_day['stage2'] = {}

        self.ref_day['stage1']['CSKS2']=0
        self.ref_day['stage1']['CSKS3']=1
        self.ref_day['stage1']['CSKS4']=4
        self.ref_day['stage1']['CSKS1']=8

        self.ref_day['stage2']['CSKS3']=-4
        self.ref_day['stage2']['CSKS4']=-1
        self.ref_day['stage2']['CSKS2']=0
        self.ref_day['stage2']['CSKS1']=8

        # Use CSK2 as reference
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

        # Critertion to remove bad acquisition based on coverage number at Evans
        min_cov = {}
        for it in range(22):
            min_cov[it] = 0
        min_cov[3] = 6
        min_cov[4] = 6
        min_cov[5] = 5
        min_cov[7] = 6
        min_cov[9] = 5.1
        min_cov[11] = 3
        min_cov[13] = 3.3
        min_cov[14] = 7
        min_cov[16] = 7
        min_cov[17] = 7
        min_cov[18] = 8

        self.min_cov = min_cov

    def get_ref_day(self, satename, day):
        stage = 'stage1'
        if satename == 'CSKS3' and day > date(2019,6,9):
            stage = 'stage2'
        if satename == 'CSKS4' and day > date(2019,5,11):
            stage = 'stage2'

        #stage='stage1'
        return self.ref_day[stage][satename]

    def date2track(self,day,sate=None,direction=None):
        '''
        day: datetime.date
        sate: satename
        direction: 'A' or 'B'
        '''
        tracks={}
        satelist = ['CSKS1','CSKS2','CSKS3','CSKS4']

        # Satellite may not be provided. Try every satellite
        for satename in satelist:
            tracks[satename]=[]
            # If satellite the given and it fails to match the current one
            if sate != None and satename != sate:
                continue
            
            # Convert the date to the corresponding date of S2
            csks2date = day - datetime.timedelta(days=self.get_ref_day(satename, day))

            for ref_date, crsp_tracks in self.ref_track.items():
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
            
            csks2date = day - datetime.timedelta(days=self.get_ref_day(satename,day))
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
            cskdate = csk2date + datetime.timedelta(days=self.get_ref_day(satename, day))

            if (cskdate-day).days % 16 == 0:
                return satename
            
    def track2date(self,track, sate=None, first=None, last=None):
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

            # Just to try every date
            today = first 
            while today < last:
                diff_day = (today - self.ref_date[track] - datetime.timedelta(days=self.get_ref_day(satename, today))).days
                if diff_day % 16 == 0: 
                    dates[satename].append(today)
                today = today + datetime.timedelta(days=1)
        
        return(dates)

    def decompose_xml(self, xmlfile):
        f = open(xmlfile)
        line = f.readlines()[1]
        f.close()
    
        # Get product name
        h5_name = line.split('<')[4].split('>')[1]
    
        # Satellite name
        sate = h5_name[0:5]
    
        # Date & time
        acq_datefmt = h5_name.split('_')[8]
        acq_time = datetime.datetime(int(acq_datefmt[:4]),int(acq_datefmt[4:6]),int(acq_datefmt[6:8]), int(acq_datefmt[8:10]),int(acq_datefmt[10:12]),int(acq_datefmt[12:14]))
   
        # find the track candidates 
        tracks = self.date2track(day = acq_time.date(), sate=sate)[sate]
   
        # direction
        direction = h5_name.split('_')[6][1]
        assert direction in ['A','D'], print("Fail to find the direction")

        # remove the track candidate whose direction is wrong
        if direction == 'A':
            track = [ i for i in tracks if i<=10 ]
        else:
            track = [ i for i in tracks if i>=11 ]
        assert len(track)==1, print("Fail to derive the track number")

        return (track[0], sate, acq_time)

if __name__=="__main__":
    CSK_Utils()
    #tracks = CSK_Utils().date2track(day=date(2017,11,21),sate='CSKS3')
    dates = CSK_Utils().track2date(track=0,sate='CSKS1',first=date(2017,11,16),last=date(2020,12,31))
    print(dates)
