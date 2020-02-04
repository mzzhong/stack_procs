workdir=/net/jokull/nobak/mzzhong/S1-Evans/track_65_pairs/20170626_20170702

azoff=${workdir}/cuDenseOffsets/20170626_20170702/filtAzimuth_20180703_v9.off
rgoff=${workdir}/cuDenseOffsets/20170626_20170702/filtRange_20180703_v9.off
target=${workdir}/merged/SLC/20170626/20170626.slc.full

azout=${workdir}/cuDenseOffsets/20170626_20170702/filtAzimuth_20180703_v9_full.off
rgout=${workdir}/cuDenseOffsets/20170626_20170702/filtRange_20180703_v9_full.off

resampleOffsets_modi.py -i ${azoff} -t ${target} -o ${azout}

resampleOffsets_modi.py -i ${rgoff} -t ${target} -o ${rgout}

#  -h, --help            show this help message and exit
#  -i INPUT, --input INPUT
#                        input file
#  -t TARGETFILE, --target_file TARGETFILE
#                        the reference file that the input will be interpolated
#                        to its size and added to it
#  -o OUTPUT, --output OUTPUT
#                        output file
#  -a ADD, --add ADD     Whether or not add the target file to the resampled
#                        offsets (defaut True)
