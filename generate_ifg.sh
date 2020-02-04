#!/bin/bash

#usage: crossmul.py [-h] -m MASTER -s SLAVE [-o PREFIX] [-a AZLOOKS]
#                   [-r RGLOOKS]
#
#Generate offset field between two Sentinel swaths
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -m MASTER, --master MASTER
#                        Master image
#  -s SLAVE, --slave SLAVE
#                        Slave image
#  -o PREFIX, --outdir PREFIX
#                        Prefix of output int and amp files
#  -a AZLOOKS, --alks AZLOOKS
#                        Azimuth looks
#  -r RGLOOKS, --rlks RGLOOKS
#                        Range looks

date1=20170228
date2=20170312
master=/net/jokull/nobak/mzzhong/S1-Evans/track_7/merged/SLC/20170228/20170228.slc.full
slave=/net/jokull/nobak/mzzhong/S1-Evans/track_7/merged/SLC/20170312/20170312.slc.full
azimuth_misreg=/net/jokull/nobak/mzzhong/S1-Evans/track_7/cuDenseOffsets/grossAzimuthfull_12d.off
range_misreg=/net/jokull/nobak/mzzhong/S1-Evans/track_7/cuDenseOffsets/grossRangefull_12d.off

fixImageXml.py -f -i $master
fixImageXml.py -f -i $slave
fixImageXml.py -f -i $azimuth_misreg
fixImageXml.py -f -i $range_misreg

outprefix='/net/jokull/nobak/mzzhong/S1-Evans/track_7/cuDenseOffsets/20170228_20170312/rs0_gpu_20170312'
runid=20180703
outsuffix='.slc.full'
interval=12
resampleMergedSlc.py --master $master --slave $slave -o ${outprefix}${outsuffix} -x $interval -a ${azimuth_misreg} -r ${range_misreg} -useGPU

crossmul -m $master -s $slave -o $prefix -a 10 -r 3
