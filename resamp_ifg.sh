#!/bin/bash
date1=20170626
date2=20170702

workdir=/net/jokull/nobak/mzzhong/S1-Evans/track_65_pairs/${date1}_${date2}

# model(0) or estimated(1)
#id='mod'
id='est'

master=${workdir}/merged/SLC/${date1}/${date1}.slc.full
slave=${workdir}/merged/SLC/${date2}/${date2}_level0.slc.full

coarse_range_off=${workdir}/merged/offset/range.off.full
coarse_azimuth_off=${workdir}/merged/offset/azimuth.off.full

interval=6

if [ "$id" = "mod" ]; then
    echo 0
    azimuth_misreg=${workdir}/cuDenseOffsets/grossAzimuthfull_${interval}d.off
    range_misreg=${workdir}/cuDenseOffsets/grossRangefull_${interval}d.off
else
    echo 1
    azimuth_misreg=${workdir}/cuDenseOffsets/${date1}_${date2}/filtAzimuth_20180703_v9_full.off
    range_misreg=${workdir}/cuDenseOffsets/${date1}_${date2}/filtRange_20180703_v9_full.off
fi

fixImageXml.py -f -i $master
fixImageXml.py -f -i $slave
fixImageXml.py -f -i $azimuth_misreg
fixImageXml.py -f -i $range_misreg


# Rubber sheeting
rs_azimuth_misreg=${workdir}/merged/offset/rs_azimuth_${id}.off.full
rs_range_misreg=${workdir}/merged/offset/rs_range.off_${id}.full


imageMath.py -e='a+b' -o $rs_azimuth_misreg -t float  --a=$coarse_azimuth_off --b=$azimuth_misreg
imageMath.py -e='a+b' -o $rs_range_misreg -t float  --a=$coarse_range_off --b=$range_misreg

fixImageXml.py -f -i $rs_azimuth_misreg
fixImageXml.py -f -i $rs_range_misreg


# Resample

#kind='_gpu'
kind=''

outprefix="${workdir}/cuDenseOffsets/${date1}_${date2}/rs_${id}_${date2}"
outsuffix='.slc.full'

resampled=${outprefix}${outsuffix}

#if [ ! -f ${outprefix}${outsuffix} ]; then
    echo 'run'
    if [ "$kind" = "_gpu" ]; then
        echo 'GPU resample'
        #resampleMergedSlc.py --master $master --slave $slave -o $resampled -a ${azimuth_misreg} -r ${range_misreg} -useGPU
        exit
    else
        echo 'CPU resample'
        resampleMergedSlc.py --master $master --slave $slave -o $resampled -a ${rs_azimuth_misreg} -r ${rs_range_misreg}
    fi
#else
#    echo 'skip'
#fi

# Interferogram
echo $resampled
echo $azimuth_misreg
echo $range_misreg
prefix="${workdir}/cuDenseOffsets/${date1}_${date2}/rs_${id}"
generateMergedIgram.py -m $master -s $resampled -o $prefix



#usage: generateIgram.py [-h] -m MASTER -s SLAVE [-f] [-i INTERFEROGRAM]
#                        [-p INTPREFIX] [-v]
#
#Use polynomial offsets and create burst by burst interferograms
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -m MASTER, --master MASTER
#                        Directory with master acquisition
#  -s SLAVE, --slave SLAVE
#                        Directory with slave acquisition
#  -f, --flatten         Flatten the interferograms with offsets if needed
#  -i INTERFEROGRAM, --interferogram INTERFEROGRAM
#                        Path for the interferogram
#  -p INTPREFIX, --interferogram_prefix INTPREFIX
#                        Prefix for the interferogram
#  -v, --overlap         Flatten the interferograms with offsets if needed


#generateIgram : 
#master : /net/jokull/nobak/mzzhong/S1-Evans/track_7/coreg_slaves/20170228
#slave : /net/jokull/nobak/mzzhong/S1-Evans/track_7/coreg_slaves/20170312
#interferogram : /net/jokull/nobak/mzzhong/S1-Evans/track_7/interferograms/20170228_20170312
#flatten : False
#interferogram_prefix : fine
#overlap : False
