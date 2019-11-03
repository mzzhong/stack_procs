#!/usr/bin/env python3

import os
import glob

# Write all the data names to this file
all_data_link = "csk_rutford_data_link.txt"
ff = open(all_data_link,'w')

# Extract coordinates
def run(filename=None, outdir=""):

    import xml.etree.ElementTree as ET

    if filename == None:
        filename = "DFDN_CSKS2_RAW_B_HI_03_HH_RA_SF_20130101030712_20130101030719.h5.xml"

    tree = ET.parse(filename)
    root = tree.getroot()

    #for child in root:
    #    print(child.tag, child.attrib)

    #for i in range(len(root)):
    #    for j in range(len(root[i])):
    #        print(i,j,root[i][j].tag, root[i][j].text)

    # Coordinates:
    bl = root[1][3].text.split()[:2]
    br = root[1][4].text.split()[:2]
    tl = root[1][5].text.split()[:2]
    tr = root[1][6].text.split()[:2]

    #elements = os.path.basename(filename).split('_')
    #coord_txt_name = outdir + '/' + elements[5] + '/' + '_'.join([elements[1], elements[5], elements[7], elements[9], elements[10].split('.')[0]]) + '.txt'

    coord_txt_name = os.path.dirname(filename) + '/' + 'footprint.txt'

    # Check if the data is on Rutford
    lons =  [ float(x) for x in [bl[1], br[1], tl[1], tr[1]] ]
    lats =  [ float(x) for x in [bl[0], br[0], tl[0], tr[0]] ]

    if min(lons)>-88 and max(lons)<-79 and min(lats)>-79.5 and max(lats)<-76:
        
        print('In Rutford')
        print(coord_txt_name)
        f = open(coord_txt_name,'w')
        f.write(bl[1] + ' ' + bl[0]+'\n')
        f.write(br[1] + ' ' + br[0]+'\n')
        f.write(tr[1] + ' ' + tr[0]+'\n')
        f.write(tl[1] + ' ' + tl[0]+'\n')
        f.write(bl[1] + ' ' + bl[0]+'\n')
        
        f.close()

        dirname = os.path.dirname(filename)
        idname = dirname.split('/')[-1]
        data_link ="https://" +  '/'.join(dirname.split('/')[7:]) + '/' + idname + '.tar.gz'
        data_name = dirname +'/' + idname + '.tar.gz'
        print(data_link)
        print(data_name)

        # Save the data name for download later
        ff.write(data_link + '\n')

    return 0

# Loop through the directory
def main():

    repo = "/net/kraken/nobak/mzzhong/CSK-Rutford/raw_data/aria-csk-dav.jpl.nasa.gov/repository/products/csk_rawb/v0.6"

    years = [2013, 2014]
    years = [2014]
    months = range(1,13)

    for y in years:
        for m in months:
            year = str(y)
            month = str(m).zfill(2)
        
            #outdir = "/net/kraken/nobak/mzzhong/CSK-Rutford/raw_data/footprints/" + year + "_" + month
        
            month_repo = repo + '/' + year + '/' + month
            #print(month_repo)
        
            for iday in range(1,32):
                day = str(iday).zfill(2)
                day_repo = month_repo +'/'+ day

                if os.path.exists(day_repo):
                    for product in os.listdir(day_repo):
                        product_repo = day_repo + '/' + product
                        id_name = os.listdir(product_repo)[0]
                        product_repo2 = product_repo +'/'+ id_name

                        if product_repo2 == "/net/kraken/nobak/mzzhong/CSK-Rutford/raw_data/aria-csk-dav.jpl.nasa.gov/repository/products/csk_rawb/v0.6/2014/03/15/CSKS2_RAW_HI_10_HH_RA_20140315154907_20140315154914/EL20140315_771310_3420217.6.2":
                            print('good')
                            #print(stop)

                        #print(product_repo2)
                        try:
                            DFDN_xml = glob.glob(product_repo2 +'/' +'DFDN*')[0]
                        except:
                            print("xml file doesn't exist")
                            print(product_repo2)
                            continue
                        run(filename=DFDN_xml)

if __name__=="__main__":
    main()
    ff.close()
