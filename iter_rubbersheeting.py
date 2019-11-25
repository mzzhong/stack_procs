#!/usr/bin/env python3
# Author: Minyan Zhong

import os
import datetime

from dense_offset import dense_offset

fmt='%Y%m%d'

class rubber_sheeting():

    def __init__(self, proj, track_number, date1, date2):

        self.proj = proj
        self.track_number = track_number
        self.date1 = date1
        self.date2 = date2

        self.date1str = date1.strftime(fmt)
        self.date2str = date2.strftime(fmt)

        if self.proj == "CSK-Rutford":
            self.stack = "stripmap"
            self.proj_folder = "/net/kraken/nobak/mzzhong/CSK-Rutford"
            self.pair_folder = os.path.join(self.proj_folder, 
                                    "track_"+str(self.track_number).zfill(3)+"_0", 
                                    'pairs', self.date1str + '_' + self.date2str)
            print(self.pair_folder)

        # Locate the master and slave SLC
        self.master_SLC = os.path.join(self.pair_folder, 'raw', self.date1str, self.date1str+ ".raw.slc")
        self.slave_SLC = os.path.join(self.pair_folder, 'raw', self.date2str, self.date2str+ ".raw.slc")
        print(self.master_SLC)
        print(self.slave_SLC)

        # Set up offset field parameters
        self.runid = 20190921

    def _set_resampled_slave_SLC(self):

        if self.iter_num==1:
            self.resampled_slave_SLC = os.path.join(self.pair_folder, 'coregSLC','Coarse',self.date2str, self.date2str+ ".slc")
            print(self.resampled_slave_SLC)

        return 0

    def _run(self, exe=False):

        exe = True

        # Preparation
        workdir = self.proj_folder
        offset = dense_offset(stack=self.stack, workdir = workdir, nproc=1, runid=self.runid, exe=exe)

        if self.proj == "CSK-Rutford":
            trackname = "track_" + str(self.track_number).zfill(3) + '_0'
            pairname = self.date1str + '_' + self.date2str
        else:
            raise Exception("undefined")

        offset.initiate_track_pair(trackname=trackname, pairname = pairname)

        offset.initiate_rubbersheeting(self.iter_num)

        # Run offset field
        offset.run_offset_pair_iter(iter_num=self.iter_num, exe=exe)
        
        # Run filter
        offset.offset_filter_pair_iter(iter_num=self.iter_num, exe=exe)
        
        return 0

    def _resample_offsetfield(self):
        
        return 0

    def _resample_SLC(self):
        
        return 0

    def run(self, iter_num):

        self.iter_num = iter_num

        self._run()

        return 0

def main():

    iter_times = 1

    proj = "CSK-Rutford"
    track_number = 201
    date1 = datetime.date(2013,9,1)
    date2 = datetime.date(2013,9,2)

    rs = rubber_sheeting(proj=proj, track_number = track_number, date1=date1, date2=date2) 

    for iter_num in range(iter_times):
        rs.run(iter_num)

if __name__=='__main__':
    main()
