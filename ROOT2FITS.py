import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def go(args):
    print " -------------"
    #print os.environ['PATH']
    #print os.environ['DYLD_LIBRARY_PATH']
    # print os.environ['PYTHONPATH']
    maptree = args.mtfile
    response = args.rsfile
    DIR ="mkdir TEMP "
    
    os.system(DIR)
    command='skymaps-maptree2fits  -i ' + str(maptree) + ' -o  ./TEMP/TEMP   --nBins 9 --first 1 '
    print command
    os.system(command)


    command1 = "aerie-apps-HealpixSigFluxMap -i  TEMP/TEMP_bin1.fits.gz TEMP/TEMP_bin2.fits.gz TEMP/TEMP_bin3.fits.gz  TEMP/TEMP_bin4.fits.gz TEMP/TEMP_bin5.fits.gz TEMP/TEMP_bin6.fits.gz TEMP/TEMP_bin7.fits.gz TEMP/TEMP_bin8.fits.gz TEMP/TEMP_bin9.fits.gz -b 1 2 3 4 5 6 7 8 9 -d " + str(response) + " -o " + str(args.out)  +".fits  --index 2.7 --negSignif --negFlux --nside  256 --pivot 7"
    print command1
    os.system(command1)
    DIR1 ="rm -rf TEMP "
    
    command2 = "python ~/software_base/aerie/src/map-maker/scripts/plotMercator.py  --tevcat-labels --coords G --milagro --interpolation  --contourmapcoord C --xyrange 5 87 -5 5 " + args.out + ".fits -o " + args.out + ".png"
    
    print command2
    os.system(command2)
    
    
    command3 = "rm -rf " + args.out + ".fits  "
    
    print command3
    os.system(command3)
    
    print "------------------------------------------------------"
    print "if there was a problem (os.system), just copy the following into bash"
    print "---"
    print DIR
    print command
    print command1
    print command2
    print DIR1
    print command3
    print "-------"

    #os.system(DIR1)



if __name__=="__main__":

    import argparse

    p = argparse.ArgumentParser(description="Example spectral fit with LiFF")
    p.add_argument("-m", "--maptreefile", dest="mtfile", help="LiFF MapTree ROOT file", required=True)
    p.add_argument("-o", "--output", dest="out", help="output name", required=True)
    p.add_argument("-r", "--responsefile", dest="rsfile", help="LiFF detector response ROOT file", required=True)

  
   

    args = p.parse_args()

    go(args)




