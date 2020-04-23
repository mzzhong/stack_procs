#!/bin/bash
fixImageXml.py -f -i lat.rdr
fixImageXml.py -f -i lon.rdr
fixImageXml.py -f -i los.rdr
looks.py -i lat.rdr -o lat_r4_a4.rdr -r 4 -a 4 &
id1=$!
looks.py -i lon.rdr -o lon_r4_a4.rdr -r 4 -a 4 &
id2=$!
looks.py -i los.rdr -o los_r4_a4.rdr -r 4 -a 4 &
id3=$!

wait $id1
wait $id2
wait $id3

fixImageXml.py -f -i lat_r4_a4.rdr
fixImageXml.py -f -i lon_r4_a4.rdr
fixImageXml.py -f -i los_r4_a4.rdr

isce2gis.py envi -i lat.rdr
isce2gis.py envi -i lon.rdr
isce2gis.py envi -i los.rdr
isce2gis.py envi -i lat_r4_a4.rdr
isce2gis.py envi -i lon_r4_a4.rdr
isce2gis.py envi -i los_r4_a4.rdr

