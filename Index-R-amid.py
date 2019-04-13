import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc("font", family="serif", size=14)

import warnings

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# This disable momentarily the printing of warnings, which you might get
# if you don't have the Fermi ST or other software installed
with warnings.catch_warnings():

    warnings.simplefilter("ignore")

    from threeML import *

# Make sure that the HAWC plugin is available

assert is_plugin_available("HAWCLike"), "HAWCLike is not available. Check your configuration"




mtfile = "PI0_residual.root"
rsfile = "response.root"
bmin = -2.
bmax = 2.
Out_text = "PI0.txt"


startBin = 1
stopBin = 9
ntransit = 1

for l in range (0,85):
    lmin = l - 0.5
    lmax = l + 0.5
    
    spectrum = Powerlaw()
    shape    = Rectangle_gal()


    source   = ExtendedSource("Extended",
                              spatial_shape=shape,
                              spectral_shape=spectrum)

    
    

    
    fluxUnit = 7.0 / (u.TeV * u.cm**2 * u.s)


    # Set spectral parameters (do it after the source definition to make sure
    # the units are handled correctly)

    spectrum.piv = 7. * u.TeV  # Pivot energy
    spectrum.piv.fix = True

    #normalization
    spectrum.K = 2.e-18
#    spectrum.K.bounds = (1e-17 * fluxUnit, 1e-10*fluxUnit)
    spectrum.K.fix = False

    #index
    spectrum.index = -2.6
    spectrum.index.bounds = (-4., -1.)
    spectrum.index.fix = False

    # Set spatial parameters (do it after the source definition to make sure
    # the units are handled correctly)
    shape.l_min = (lmin) * u.degree
    shape.l_min.fix = True
    #shape.lon0.bounds = ()
    
    shape.l_max = (lmax) * u.degree
    shape.l_max.fix = True
    #shape.lat0.bounds = ()


    shape.b_min =  (bmin) * u.degree
    shape.b_min.fix = True

    shape.b_max = (bmax) * u.degree
    shape.b_max.fix = True

    shape.K = 1.0
    shape.K.fix = True


    # Set up a likelihood model using the source
    lm = Model(source)
    llh = HAWCLike("Extended", mtfile, rsfile,fullsky=True) #if fullsky you need to define a ROI
    llh.set_active_measurements(startBin, stopBin)
    # llh.set_ROI(ra, dec, ROI ,True) #the whole map
    llh.set_strip_ROI(8. , 21 , -6, 6, fixed_ROI=False, galactic=True)

    # Double check the free parameters
    print("Likelihood model:\n")
    print "----------------------------------------"
    lm.display(complete=True)

    # Set up the likelihood and run the fit
    print("Performing likelihood fit...\n")
    datalist = DataList(llh)
    jl = JointLikelihood(lm, datalist, verbose=True)
    

    try:
        jl.set_minimizer("ROOT")
        result = jl.fit()

    except AttributeError:
        jl.set_minimizer("minuit")
        result = jl.fit()
    
    print result
    # Print the TS, significance, and fit parameters, and then plot stuff
    print("\nTest statistic:")
    TS = llh.calc_TS()
    sigmaa = np.sqrt(TS)

    print("Test statistic: %g" % TS)
    print("Significance:   %g\n" % sigmaa)
    print "-----------------------"
#parsing and cleaning for the same output format as get_errors()
    print "parsing and cleaning for the same output format as get_errors"
    item=[]
    print result
    word=result[0].to_string().split()
    print 'word', word
    # item.append(word[5])
    #for i in range (0, 20):
#   print "for i=",i, "===>",  word[i]

    G=open( Out_text , 'a')
    G.write("\n%s\t" % str(l))
    G.write("%s\t" % str(lmin))
    G.write("%s\t" % str(lmax))
    G.write("%s\t" % str(bmin))
    G.write("%s\t" % str(bmax))
    G.write("%s\t" % str(float(word[6])))
    G.write("%s\t" % str(float(word[7])))
    G.write("%s\t" % str(float(word[8])))
    G.write("%s\t" % str(word[16]))
    G.write("%s\t" % str(word[17]))
    G.write("%s\t" % str(word[18]))
    G.close()

G=open( Out_text , 'a')
G.write("\n%s\t" % "l")
G.write("%s\t" % "lmin")
G.write("%s\t" % "lmax")
G.write("%s\t" % "bmin")
G.write("%s\t" % "bmax")
G.write("%s\t" % "Norm")
G.write("%s\t" % "Neg_Norm_Err")
G.write("%s\t" % "Pos_Norm_Err")
G.write("%s\t" % "Index")
G.write("%s\t" % "Neg_Index_Err")
G.write("%s\t" % "Pos_Index_Err")
G.write("%s\t" % "Finished")
G.close()
