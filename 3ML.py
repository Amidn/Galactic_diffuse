import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import datetime


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


def go(args):
    startt = datetime.datetime.now()
   
    fluxUnit = 1. / (u.TeV * u.cm**2 * u.s)
    
    mtfile = './maptree_256.root'
    rsfile = './response.root'
    output = "source"
    NF = 1
    
    
    spectrum1 = Powerlaw()
    source1  = PointSource("A_1809",
                              272.4,
                              -19.0,
    spectral_shape=spectrum1)
          
    source1.position.ra.fix = True
    source1.position.ra.bounds = (271.5,273.5)
                              
    source1.position.dec.fix = True
    source1.position.dec.bounds = (-20., -18.5)
                                  
    spectrum1.K = 4.3e-14 * fluxUnit
    spectrum1.K.fix =True
                                              
    spectrum1.piv = 7.0 * u.TeV
    spectrum1.piv.fix =True
                                                      
    spectrum1.index = -2.51
    spectrum1.index.fix = True
    #---------------
    
  
    # --------------------
    spectrum2 = Powerlaw()
    source2  = PointSource("A_F1813_173",
                            273.47 ,
                            -17.35,
                            spectral_shape=spectrum2)
    source2.position.ra.fix = True
    source2.position.ra.bounds = (273.0,273.8)
    source2.position.dec.fix = True
    source2.position.dec.bounds = (-17.7, -17.0)

    spectrum2.K = 5.0e-14 * fluxUnit
    spectrum2.K.fix = True
    
    spectrum2.piv = 7.0 * u.TeV
    spectrum2.piv.fix = True
    
    spectrum2.index = -2.46
    spectrum2.index.fix = True
    
    

    
    
    spectrum3  = Powerlaw()
    source3    = PointSource("A_j1819",
                                274.83,
                                -15.13,
                                spectral_shape=spectrum3)

    spectrum3.K = 1.5e-14 * fluxUnit
    spectrum3.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit)
    spectrum3.K.fix = True
    
    spectrum3.index = -2.64
    spectrum3.index.bounds = (-4., 0.0)
    spectrum3.index.fix = True
    
    spectrum3.piv = 7.0 * u.TeV
    spectrum3.piv.fix = True
    
    
    # -------------
    
    spectrum4 = Powerlaw()
    source4   = PointSource("A_J1813_125",
                              273.38,
                              -12.56,
                              spectral_shape=spectrum4)
    spectrum4.K = 1.3e-14 * fluxUnit
    spectrum4.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit)
    spectrum4.K.fix = True
                              
    spectrum4.index = -2.610
    spectrum4.index.bounds = (-4., 0.0)
    spectrum4.index.fix = True
                              
    spectrum4.piv = 7.0 * u.TeV
    spectrum4.piv.fix = True
                              
                              #-------------------------
    spectrum5 = Powerlaw()
    shape5 = Gaussian_on_sphere()
    source5 = ExtendedSource("A_J1825", spatial_shape=shape5  ,spectral_shape=spectrum5)
    shape5.lon0 = 276.46 * u.degree
    shape5.lon0.fix = True
    shape5.lat0 = -13.4 * u.degree
    shape5.lat0.fix = True
                              
    shape5.sigma = 0.314 * u.degree
    shape5.sigma.fix = True
    shape5.sigma.max_value = 3.
                              
    spectrum5.K = 1.58e-13 * fluxUnit
    spectrum5.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrum5.K.fix = True
                              
    spectrum5.piv = 7.0 * u.TeV
    spectrum5.piv.fix = True
                              
    spectrum5.index = -2.446
    spectrum5.index.bounds = (-4., 0.0)
    spectrum5.index.fix = True
                              
                              # - ----------------------------- LS-5039
                              
    spectrum6 = Powerlaw()
    source6  = PointSource("A_LS5039",
                                276.59,
                                -14.616,
                                spectral_shape=spectrum6)
    source6.position.ra.free = False
    source6.position.ra.bounds = (275.5,277.5)
                              
    source6.position.dec.free = False
    source6.position.dec.bounds = (-15.5, -13.0)
                              
                              
    spectrum6.K = 3.2e-14 * fluxUnit
    spectrum6.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit)
    spectrum6.K.fix = True
                              
    spectrum6.index = -2.65
    spectrum6.index.bounds = (-4., 0.0)
    spectrum6.index.fix = True
                              
    spectrum6.piv = 7.0 * u.TeV
    spectrum6.piv.fix = True

 
    #------------------------------------------------------------
    #set the GDE
    shapeA = SpatialTemplate_2D()
    shapeA.load_file('A_2D_pi0_DRAGON_BASE_v4_Amid.fits') #DRAGON, pi0
        
    spectrumA = Powerlaw()
    spectrumA.K = 7.1e-18 # * fluxUnit
    spectrumA.K.fix = True
    # spectrumGDE4.K.bounds = (1e-11, 1e-4)

    spectrumA.index = -2.613
    spectrumA.index.fix = True
    # spectrumGDE4.index.bounds = (-4, -1)

    spectrumA.piv = 7.0e9 #* u.TeV
    spectrumA.piv.fix = True

    # spectrumGDE4.display()
    # shapeGDE4.K.fix = True
    
    shapeA.hash = 1.
    shapeA.hash.fix = True
    
    shapeA.K = 1.
    shapeA.K.fix = True
    sourceA = ExtendedSource("A_pi0", spatial_shape=shapeA, spectral_shape = spectrumA )
    
    # ---------------------------------------------------------------
    
    
    
    
                                                                                     
                                                                                     
                                                                                     
    # ----------
                                                                                     
    spectrumB1 = Powerlaw()
    shapeB1    = Gaussian_on_sphere()
    sourceB1   = ExtendedSource("B_E1831",
                                spatial_shape=shapeB1 ,
                                spectral_shape=spectrumB1)
    shapeB1.lon0 = 277.91  * u.degree
    shapeB1.lon0.fix = True
    shapeB1.lat0 = -9.71 * u.degree
    shapeB1.lat0.fix = True
                                                                                                                                 
    shapeB1.sigma = 0.73 * u.degree
    shapeB1.sigma.fix = True
    shapeB1.sigma.max_value = 5.
                                                                                                                                             
    spectrumB1.K = 9.4e-14 * fluxUnit
    spectrumB1.K.fix = True
    
    spectrumB1.piv = 7.0 * u.TeV
    spectrumB1.piv.fix = True
    
    spectrumB1.index = -2.57
    spectrumB1.index.fix = True
    

   # ----
    spectrumB2 = Powerlaw()
    shapeB2    = Gaussian_on_sphere()
    sourceB2   = ExtendedSource("B_E1837",
                                 spatial_shape=shapeB2,
                                 spectral_shape=spectrumB2)
    shapeB2.lon0 = 279.38  * u.degree
    shapeB2.lon0.fix = True
    shapeB2.lon0.bounds = (278.4, 280.4)
    shapeB2.lat0 =  -6.624   * u.degree
    shapeB2.lat0.fix = True
    shapeB2.lat0.bounds = (-7.54, -5.54)

    shapeB2.sigma = 0.364 * u.degree
    shapeB2.sigma.fix = True
    shapeB2.sigma.max_value = 5.
    
    spectrumB2.K = 1.30e-13 * fluxUnit
    spectrumB2.K.fix = True
    
    spectrumB2.piv = 7.0 * u.TeV
    spectrumB2.piv.fix = True
    
    spectrumB2.index = -2.667
    spectrumB2.index.fix = True
    # --------H1841 extended is prefered
  
  
    spectrum9 = Powerlaw()
    shape9    = Gaussian_on_sphere()
    source9   = ExtendedSource("B_H1841E",
                                 spatial_shape=shape9,
                                 spectral_shape=spectrum9)
    shape9.lon0 = 280.1  * u.degree
    shape9.lon0.fix = True
    shape9.lon0.bounds = (279., 281.)
    shape9.lat0 = -5.477  * u.degree
    shape9.lat0.fix = True
    shape9.lat0.bounds = (-6.45, -4.45)
    
    shape9.sigma = 0.404 * u.degree
    shape9.sigma.fix = True
    shape9.sigma.max_value = 5.
    
    spectrum9.K = 1.03e-13 * fluxUnit
    spectrum9.K.fix = True
    
    spectrum9.piv = 7.0 * u.TeV
    spectrum9.piv.fix = True
    
    spectrum9.index = -2.546
    spectrum9.index.fix = True
   #-------------------------------- extended is prefered

    spectrum10 = Powerlaw()
    shape10 = Gaussian_on_sphere()
    source10 = ExtendedSource("B_Ej184333",
                              spatial_shape=shape10,
                              spectral_shape=spectrum10)
    shape10.lon0 = 280.99 * u.degree
    shape10.lon0.fix = True
    shape10.lat0 = -3.32 * u.degree
    shape10.lat0.fix = True
    
    shape10.sigma = 0.385 * u.degree
    shape10.sigma.fix = True
    shape10.sigma.max_value = 5.
    
    spectrum10.K = 7.7e-14 * fluxUnit
    spectrum10.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrum10.K.fix = True
    
    spectrum10.piv = 7.0 * u.TeV
    spectrum10.piv.fix = True
    
    spectrum10.index = -2.399
    spectrum10.index.bounds = (-4., 0.0)
    spectrum10.index.fix = True


    spectrum11 = Powerlaw()
    shape11 = Gaussian_on_sphere()
    source11 = ExtendedSource("B_Ej1847",
                             spatial_shape=shape11,
                             spectral_shape=spectrum11)
    shape11.lon0 = 281.91 * u.degree
    shape11.lon0.fix = True
    
    shape11.lat0 = -1.79 * u.degree
    shape11.lat0.fix = True
    
    shape11.sigma = 0.36 * u.degree
    shape11.sigma.fix = True
    shape11.sigma.max_value = 5.
    
    spectrum11.K = 3.8e-14 * fluxUnit
    spectrum11.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrum11.K.fix = True

    spectrum11.piv = 7.0 * u.TeV
    spectrum11.piv.fix = True
    
    spectrum11.index = -2.57
    spectrum11.index.bounds = (-4., 0.0)
    spectrum11.index.fix = True
        

    spectrum12 = Powerlaw()
    source12   = PointSource("B_J1843001",
                               280.99,
                               -0.19,
                               spectral_shape=spectrum12)

    spectrum12.K = 4.8e-15  * fluxUnit
    spectrum12.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit)
    spectrum12.K.fix = True

    spectrum12.index = -2.44
    spectrum12.index.bounds = (-4., 0.0)
    spectrum12.index.fix = True

    spectrum12.piv = 7.0 * u.TeV
    spectrum12.piv.fix = True


    spectrum13 = Powerlaw()
    shape13 = Gaussian_on_sphere()
    source13 = ExtendedSource("B_Ej1849",
                                spatial_shape=shape13,
                                spectral_shape=spectrum13)

    shape13.lon0 = 282.35 * u.degree
    shape13.lon0.fix = True
    shape13.lon0.bounds = (281.3, 283.)
    shape13.lat0 = 0.15 * u.degree
    shape13.lat0.fix = True
    shape13.lat0.bounds= (-1.0, 1.0)

    shape13.sigma = 0.11 * u.degree
    shape13.sigma.fix = True
    shape13.sigma.max_value = 1.

    spectrum13.K = 4.5e-14 * fluxUnit
    spectrum13.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrum13.K.fix = True

    spectrum13.piv = 7.0 * u.TeV
    spectrum13.piv.fix = True

    spectrum13.index = -2.28
    spectrum13.index.bounds = (-4., 0.0)
    spectrum13.index.fix = True



    spectrum13_1 = Powerlaw()
    shape13_1 = Gaussian_on_sphere()
    source13_1 = ExtendedSource("B_HessEj1849",
                                spatial_shape=shape13_1,
                                spectral_shape=spectrum13_1)

    shape13_1.lon0 = 282.257625 * u.degree
    shape13_1.lon0.fix = True
    shape13_1.lon0.bounds = (281.3, 283.)
    shape13_1.lat0 = -0.021972 * u.degree
    shape13_1.lat0.fix = True
    shape13_1.lat0.bounds= (-1.0, 1.0)
    
    shape13_1.sigma = 0.69 * u.degree
    shape13_1.sigma.fix = False
    shape13_1.sigma.max_value = 1.

    spectrum13_1.K = 4.9e-13 * fluxUnit
    spectrum13_1.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrum13_1.K.fix = False

    spectrum13_1.piv = 7.0 * u.TeV
    spectrum13_1.piv.fix = True

    spectrum13_1.index = -2.62
    spectrum13_1.index.bounds = (-4., 0.0)
    spectrum13_1.index.fix = False

# ____________________
    spectrum14 = Powerlaw()
    shape14 = Gaussian_on_sphere()
    source14 = ExtendedSource("B_HE1852",
                                 spatial_shape=shape14,
                                 spectral_shape=spectrum14)
    shape14.lon0 = 283.05416667 * u.degree
    shape14.lon0.fix = True
    shape14.lat0 = -0.00638889 * u.degree
    shape14.lat0.fix = True
    
    shape14.sigma = 0.48 * u.degree
    shape14.sigma.fix = True
    shape14.sigma.max_value = 5.

    spectrum14.K = 3.e-14 * fluxUnit
    spectrum14.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrum14.K.fix = True
    
    spectrum14.piv = 7.0 * u.TeV
    spectrum14.piv.fix = True
    
    spectrum14.index = -2.39
    spectrum14.index.bounds = (-4., 0.0)
    spectrum14.index.fix = True


    
    
    
    #=================================
    
    shapeB = SpatialTemplate_2D()
    shapeB.load_file('B_2D_pi0_DRAGON_BASE_v4_Amid.fits') #DRAGON, pi0

    spectrumB = Powerlaw()
    spectrumB.K = 4.85e-18 # * fluxUnit
    spectrumB.K.fix = True
    # spectrumGDE4.K.bounds = (1e-11, 1e-4)
    
    spectrumB.index = -2.608
    spectrumB.index.fix = True
    # spectrumGDE4.index.bounds = (-4, -1)
    
    spectrumB.piv = 7.0e9 #* u.TeV
    spectrumB.piv.fix = True
    
    # spectrumGDE4.display()
    # shapeGDE4.K.fix = True
    
    shapeB.hash = 1.
    shapeB.hash.fix = True
    
    shapeB.K = 1.
    shapeB.K.fix = True
    
    sourceB = ExtendedSource("B_pi0", spatial_shape=shapeB, spectral_shape = spectrumB )
    #---------------------------
        # region C
    #--------------------------------- Extended one is preferred
    

    # this one is selected
    
    spectrumC5 = Powerlaw()
    shapeC5 = Gaussian_on_sphere()
    sourceC5 = ExtendedSource("C_Ej1857",
                               spatial_shape=shapeC5,
                               spectral_shape=spectrumC5)
    shapeC5.lon0 = 284.33 * u.degree
    shapeC5.lon0.fix = True
    shapeC5.lat0 = 2.87 * u.degree
    shapeC5.lat0.fix = True
                               
    shapeC5.sigma = 0.4 * u.degree
    shapeC5.sigma.fix = True
    shapeC5.sigma.max_value = 5.
                               
    spectrumC5.K = 5.8e-14 * fluxUnit
    spectrumC5.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrumC5.K.fix = True
                               
    spectrumC5.piv = 7.0 * u.TeV
    spectrumC5.piv.fix = True
                               
    spectrumC5.index = -2.66
    spectrumC5.index.bounds = (-4., 0.0)
    spectrumC5.index.fix = True
    
    
    
    spectrumC6 = Powerlaw()
    shapeC6 = Gaussian_on_sphere()
    sourceC6 = ExtendedSource("C_Ej1855",
                                spatial_shape=shapeC6,
                                spectral_shape=spectrumC6)
    shapeC6.lon0 = 283.89 * u.degree
    shapeC6.lon0.fix = True
    shapeC6.lat0 = 4.89 * u.degree
    shapeC6.lat0.fix = True
                                
    shapeC6.sigma = 0.1 * u.degree
    shapeC6.sigma.fix = False
    shapeC6.sigma.max_value = 5.
                                
    spectrumC6.K = 2.0e-14 * fluxUnit
    spectrumC6.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrumC6.K.fix = False
                                
    spectrumC6.piv = 7.0 * u.TeV
    spectrumC6.piv.fix = True
                                
    spectrumC6.index = -2.34
    spectrumC6.index.bounds = (-4., 0.0)
    spectrumC6.index.fix = False
    

    #---------------------------
    
    spectrumC1 = Powerlaw()
    sourceC1  = PointSource("C_2HWC_J1902_048",
                             285.51 ,
                             4.86,
                             spectral_shape=spectrumC1)
    sourceC1.position.ra.free = False
    # sourceF13.position.ra.bounds = (287,289)
                             
    sourceC1.position.dec.free = False
     # sourceF13.position.dec.bounds = (10., 11.)
                             
                             
    spectrumC1.K = 3.74e-15* fluxUnit
    spectrumC1.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit)
    spectrumC1.K.fix = False
                             
    spectrumC1.index = -2.3
    spectrumC1.index.bounds = (-4., 0.0)
    spectrumC1.index.fix = False
                             
    spectrumC1.piv = 7.0 * u.TeV
    spectrumC1.piv.fix = True
                             
    
    #---------------


    spectrumC2 = Powerlaw()
    shapeC2= Gaussian_on_sphere()
    sourceC2= ExtendedSource("C_j1908",
                               spatial_shape=shapeC2,
                               spectral_shape=spectrumC2)
    shapeC2.lon0 = 287.05 * u.degree
    shapeC2.lon0.fix = True
    shapeC2.lat0 = 6.39 * u.degree
    shapeC2.lat0.fix = True
                               
    shapeC2.sigma = 0.474 * u.degree
    shapeC2.sigma.fix = True
    shapeC2.sigma.max_value = 5.
                               
    spectrumC2.K = 8.7e-14 * fluxUnit
    spectrumC2.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrumC2.K.fix = True
                               
    spectrumC2.piv = 7.0 * u.TeV
    spectrumC2.piv.fix = True
                               
    spectrumC2.index = -2.316
    spectrumC2.index.bounds = (-4., 0.0)
    spectrumC2.index.fix = True
    
    #--------------------
    
    spectrumC3 = Powerlaw()
    sourceC3  = PointSource("C_2HWC_J1907_084",
                               286.79 ,
                              8.50,
                             spectral_shape=spectrumC3)
    sourceC3.position.ra.free = False
   # sourceF13.position.ra.bounds = (287,289)

    sourceC3.position.dec.free = False
    # sourceF13.position.dec.bounds = (10., 11.)
    
    
    spectrumC3.K = 4.3e-15 * fluxUnit
    spectrumC3.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit)
    spectrumC3.K.fix = True
    
    spectrumC3.index = -2.78
    spectrumC3.index.bounds = (-4., 0.0)
    spectrumC3.index.fix = True
    
    spectrumC3.piv = 7.0 * u.TeV
    spectrumC3.piv.fix = True
    #-------------------------
    
    spectrumC4 = Powerlaw()
    sourceC4  = PointSource("C_j1913_050",
                             288.41,
                             5.08,
                             spectral_shape=spectrumC4)
    sourceC4.position.ra.free = False
    # sourceF13.position.ra.bounds = (287,289)
                             
    sourceC4.position.dec.free = False
    # sourceF13.position.dec.bounds = (10., 11.)
                             
                             
    spectrumC4.K = 5.7e-15 * fluxUnit
    spectrumC4.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit)
    spectrumC4.K.fix = True
                             
    spectrumC4.index = -2.37
    spectrumC4.index.bounds = (-4., 0.0)
    spectrumC4.index.fix = True
                             
    spectrumC4.piv = 7.0 * u.TeV
    spectrumC4.piv.fix = True
    

    shapeC = SpatialTemplate_2D()
    shapeC.load_file('C_2D_pi0_DRAGON_BASE_v4_Amid.fits') #DRAGON, pi0
    
    spectrumC = Powerlaw()
    spectrumC.K = 3.9e-18 # * fluxUnit
    spectrumC.K.fix = True
    # spectrumGDE4.K.bounds = (1e-11, 1e-4)
    
    spectrumC.index = -2.537
    spectrumC.index.fix = True
    # spectrumGDE4.index.bounds = (-4, -1)
  
    spectrumC.piv = 7.0e9 #* u.TeV
    spectrumC.piv.fix = True
    
    shapeC.hash = 1.
    shapeC.hash.fix = True
    
    shapeC.K = 1.
    shapeC.K.fix = True
    
    sourceC = ExtendedSource("C_pi0", spatial_shape=shapeC, spectral_shape = spectrumC )
    



   
    spectrumD1 = Powerlaw()
    shapeD1 = Gaussian_on_sphere()
    sourceD1 = ExtendedSource("D_Ej1912_103",
                            spatial_shape=shapeD1,
                            spectral_shape=spectrumD1)
    shapeD1.lon0 = 288.02 * u.degree
    shapeD1.lon0.fix = True
    shapeD1.lat0 = 10.39 * u.degree
    shapeD1.lat0.fix = True
    
    shapeD1.sigma = 0.36 * u.degree
    shapeD1.sigma.fix = True
    shapeD1.sigma.max_value = 5.
    
    spectrumD1.K = 2.59e-14 * fluxUnit
    spectrumD1.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrumD1.K.fix = True
    
    spectrumD1.piv = 7.0 * u.TeV
    spectrumD1.piv.fix = True
    
    spectrumD1.index = -2.65
    spectrumD1.index.bounds = (-4., 0.0)
    spectrumD1.index.fix = True
    

    
    spectrumD2 = Powerlaw()
    shapeD2 = Gaussian_on_sphere()
    sourceD2 = ExtendedSource("D_Ej1914_118",
                                spatial_shape=shapeD2,
                                spectral_shape=spectrumD2)
    shapeD2.lon0 = 288.63 * u.degree
    shapeD2.lon0.fix = True
    shapeD2.lat0 = 11.83 * u.degree
    shapeD2.lat0.fix = True
                                
    shapeD2.sigma = 0.2 * u.degree
    shapeD2.sigma.fix = True
    shapeD2.sigma.max_value = 5.
                                
    spectrumD2.K = 1.04e-14 * fluxUnit
    spectrumD2.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrumD2.K.fix = True
                                
    spectrumD2.piv = 7.0 * u.TeV
    spectrumD2.piv.fix = True
                                
    spectrumD2.index = -2.48
    spectrumD2.index.bounds = (-4., 0.0)
    spectrumD2.index.fix = True
                                
                                
    # ---------------------------------------------- the fit was not succesful
                                
    spectrumD3 = Powerlaw()
    sourceD3  = PointSource("D_j1914_163",
                                288.72,
                                16.37,
                                spectral_shape=spectrumD3)
    sourceD3.position.ra.free = False
    sourceD3.position.ra.bounds = (288.5,289)
                                
    sourceD3.position.dec.free = False
    sourceD3.position.dec.bounds = (16, 17)
                                
                                
    spectrumD3.K = 3.2e-15 * fluxUnit
    spectrumD3.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit)
    spectrumD3.K.fix = True
                                
    spectrumD3.index = -2.38
    spectrumD3.index.bounds = (-4., 0.0)
    spectrumD3.index.fix = True
                                
    spectrumD3.piv = 7.0 * u.TeV
    spectrumD3.piv.fix = True
                                
    # -------------------------------------------------
                                
                                

                                
    #-----
                                
    spectrumD4 = Powerlaw()
    shapeD4 = Gaussian_on_sphere()
    sourceD4 = ExtendedSource("D_Ej1922_141",
                                spatial_shape=shapeD4,
                                spectral_shape=spectrumD4)
    shapeD4.lon0 = 290.70 * u.degree
    shapeD4.lon0.fix = True
    shapeD4.lat0 = 14.17 * u.degree
    shapeD4.lat0.fix = True
                                
    shapeD4.sigma = 0.17 * u.degree
    shapeD4.sigma.fix = True
    shapeD4.sigma.max_value = 5.
                                
    spectrumD4.K = 9.2e-15 * fluxUnit
    spectrumD4.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrumD4.K.fix = True
                                
    spectrumD4.piv = 7.0 * u.TeV
    spectrumD4.piv.fix = True
                                
    spectrumD4.index = -2.42
    spectrumD4.index.bounds = (-4., 0.0)
    spectrumD4.index.fix = True
                                
                                
                                
    # -----------------------------------------------
                                
    spectrumD5 = Powerlaw()
    shapeD5 = Gaussian_on_sphere()
    sourceD5 = ExtendedSource("D_j1928",
                                spatial_shape=shapeD5,
                                spectral_shape=spectrumD5)
    shapeD5.lon0 = 292.10 * u.degree
    shapeD5.lon0.fix = True
    shapeD5.lat0 = 17.82 * u.degree
    shapeD5.lat0.fix = True
                                
    shapeD5.sigma = 0.21 * u.degree
    shapeD5.sigma.fix = True
    shapeD5.sigma.max_value = 5.
                                
    spectrumD5.K = 1.07e-14 * fluxUnit
    spectrumD5.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrumD5.K.fix = True
                                
    spectrumD5.piv = 7.0 * u.TeV
    spectrumD5.piv.fix = True
                                
    spectrumD5.index = -2.33
    spectrumD5.index.bounds = (-4., 0.0)
    spectrumD5.index.fix = True
                                

    spectrumD6 = Powerlaw()
    shapeD6= Gaussian_on_sphere()
    sourceD6= ExtendedSource("D_Ej1930",
                                spatial_shape=shapeD6,
                                spectral_shape=spectrumD6)
    shapeD6.lon0 = 292.50 * u.degree
    shapeD6.lon0.fix = True
    shapeD6.lat0 = 18.88 * u.degree
    shapeD6.lat0.fix = True
                                
    shapeD6.sigma = 0.48 * u.degree
    shapeD6.sigma.fix = True
    shapeD6.sigma.max_value = 5.
                                
    spectrumD6.K = 1.8e-14 * fluxUnit
    spectrumD6.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrumD6.K.fix = True
                                
    spectrumD6.piv = 7.0 * u.TeV
    spectrumD6.piv.fix = True
                                
    spectrumD6.index = -2.37
    spectrumD6.index.bounds = (-4., 0.0)
    spectrumD6.index.fix = True

    shapeD = SpatialTemplate_2D()
    shapeD.load_file('D_2D_pi0_DRAGON_BASE_v4_Amid.fits') #DRAGON, pi0
    
    spectrumD = Powerlaw()
    spectrumD.K = 2.65e-18 # * fluxUnit
    spectrumD.K.fix = True
    # spectrumGDE4.K.bounds = (1e-11, 1e-4)
    
    spectrumD.index = -2.658
    spectrumD.index.fix = True
    # spectrumGDE4.index.bounds = (-4, -1)
  
    spectrumD.piv = 7.0e9 #* u.TeV
    spectrumD.piv.fix = True
    
    shapeD.hash = 1.
    shapeD.hash.fix = True
    
    shapeD.K = 1.
    shapeD.K.fix = True
    
    sourceD = ExtendedSource("D_pi0", spatial_shape=shapeD, spectral_shape = spectrumD )



 #===================


    spectrumE1 = Powerlaw()
    sourceE1 = PointSource("E_j1939",
                               294.83   ,
                               24.05,
                               spectral_shape=spectrumE1)
    sourceE1.position.ra.free = False
  #  sourceF19.position.ra.bounds = (293.5,295.5)
    
    sourceE1.position.dec.free = False
 #   sourceF19.position.dec.bounds = (23., 25.0)


    spectrumE1.K = 4.6e-15 * fluxUnit
    spectrumE1.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) 
    spectrumE1.K.fix = True

    spectrumE1.index = -2.68
    spectrumE1.index.bounds = (-4., 0.0)
    spectrumE1.index.fix = True

    spectrumE1.piv = 7.0 * u.TeV
    spectrumE1.piv.fix = True

# ---------------------------------------------------

    spectrumE2 = Powerlaw()
    sourceE2  = PointSource("j1947",
                               296.81,
                               24.99,
                               spectral_shape=spectrumE2)
    sourceE2.position.ra.free = False
  #  sourceF20.position.ra.bounds = (295.5,297.5)
    
    sourceE2.position.dec.free = False
  #  sourceF20.position.dec.bounds = (18., 19.5)


    spectrumE2.K =2.6e-15 * fluxUnit
    spectrumE2.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) 
    spectrumE2.K.fix = True

    spectrumE2.index = -2.37
    spectrumE2.index.bounds = (-4., 0.0)
    spectrumE2.index.fix = True

    spectrumE2.piv = 7.0 * u.TeV
    spectrumE2.piv.fix = True
    
    
    
    #----------------------------------------------
    
    
    
    spectrumE3= Powerlaw()
    shapeE3 = Gaussian_on_sphere()
    sourceE3 = ExtendedSource("E2j1949",
                                spatial_shape=shapeE3,
                                spectral_shape=spectrumE3)
    shapeE3.lon0 = 297.42  * u.degree
    shapeE3.lon0.fix = True
    shapeE3.lat0 = 24.46 * u.degree
    shapeE3.lat0.fix = True
    
    shapeE3.sigma = 0.31 * u.degree
    shapeE3.sigma.fix = False
    shapeE3.sigma.max_value = 5.
    
    spectrumE3.K = 2.4e-13 * fluxUnit
    spectrumE3.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrumE3.K.fix = False
    
    spectrumE3.piv = 7.0 * u.TeV
    spectrumE3.piv.fix = True
    
    spectrumE3.index = -2.54
    spectrumE3.index.bounds = (-4., 0.0)
    spectrumE3.index.fix = False


    spectrumE3_1 = Powerlaw()
    sourceE3_1  = PointSource("E2j1949",
                               297.42,
                               24.46,
                               spectral_shape=spectrumE3_1)
    sourceE3_1.position.ra.free = False
  #  sourceF20.position.ra.bounds = (295.5,297.5)
    
    sourceE3_1.position.dec.free = False
  #  sourceF20.position.dec.bounds = (18., 19.5)


    spectrumE3_1.K =2.2e-15 * fluxUnit
    spectrumE3_1.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) 
    spectrumE3_1.K.fix = True

    spectrumE3_1.index = -2.5
    spectrumE3_1.index.bounds = (-4., 0.0)
    spectrumE3_1.index.fix = True

    spectrumE3_1.piv = 7.0 * u.TeV
    spectrumE3_1.piv.fix = True
    
    
    #=================================
    
    shapeE = SpatialTemplate_2D()
    shapeE.load_file('E_2D_pi0_DRAGON_BASE_v4_Amid.fits') #DRAGON, pi0
    
    spectrumE = Powerlaw()
    spectrumE.K = 1.36e-18 # * fluxUnit
    spectrumE.K.fix = True
    # spectrumGDE4.K.bounds = (1e-11, 1e-4)
    
    spectrumE.index = -2.61
    spectrumE.index.fix = True
    # spectrumGDE4.index.bounds = (-4, -1)
  
    spectrumE.piv = 7.0e9 #* u.TeV
    spectrumE.piv.fix = True
    
    shapeE.hash = 1.
    shapeE.hash.fix = True
    
    shapeE.K = 1.
    shapeE.K.fix = True
    
    sourceE = ExtendedSource("E_pi0", spatial_shape=shapeE, spectral_shape = spectrumE )
    



    spectrumF1 = Powerlaw()
    sourceF1  = PointSource("F_j1953",
                               298.30  ,
                               28.84,
                               spectral_shape=spectrumF1)
    sourceF1.position.ra.fix = True
  #  sourceF1.position.ra.bounds = (291.5,293.5)
    
    sourceF1.position.dec.fix = True
   # sourceF1.position.dec.bounds = (18., 19.5)


    spectrumF1.K = 4.4e-15 * fluxUnit
    spectrumF1.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) 
    spectrumF1.K.fix = True

    spectrumF1.index = -2.38
    spectrumF1.index.bounds = (-4., 0.0)
    spectrumF1.index.fix = True

    spectrumF1.piv = 7.0 * u.TeV
    spectrumF1.piv.fix = True

# ----------------------------------------------- this source is to faint
    spectrumF2 = Powerlaw()
    sourceF2  = PointSource("F_j1951",
                               297.99,
                               29.40,
                               spectral_shape=spectrumF2)
    sourceF2.position.ra.free = False
  #  sourceF2.position.ra.bounds = (291.5,293.5)
    
    sourceF2.position.dec.free = False
 #   sourceF2.position.dec.bounds = (18., 29.5)


    spectrumF2.K = 6.4e-15 * fluxUnit
    spectrumF2.K.bounds = (1e-18 * fluxUnit, 1e-10 * fluxUnit) 
    spectrumF2.K.fix = True

    spectrumF2.index = -2.43
    spectrumF2.index.bounds = (-4., 0.0)
    spectrumF2.index.fix = True

    spectrumF2.piv = 7.0 * u.TeV
    spectrumF2.piv.fix = True



# -------------------------------------------------

    spectrumF3 = Powerlaw()
    sourceF3  = PointSource("F_j1957",
                               299.31,
                               29.14,
                               spectral_shape=spectrumF3)
    sourceF3.position.ra.free = False
   # sourceF3.position.ra.bounds = (291.5,293.5)
    
    sourceF3.position.dec.free = False
   # sourceF3.position.dec.bounds = (18., 29.5)


    spectrumF3.K = 2.7e-15 * fluxUnit
    spectrumF3.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) 
    spectrumF3.K.fix = True

    spectrumF3.index = -2.37
    spectrumF3.index.bounds = (-4., 0.0)
    spectrumF3.index.fix = True

    spectrumF3.piv = 7.0 * u.TeV
    spectrumF3.piv.fix = True

# --------------------------------------------------

    spectrumF4 = Powerlaw()
    sourceF4  = PointSource("f_j2005",
                               301.42 ,
                               31.13,
                               spectral_shape=spectrumF4)
    sourceF4.position.ra.free = False
   # sourceF4.position.ra.bounds = (291.5,293.5)
    
    sourceF4.position.dec.free = False
   # sourceF4.position.dec.bounds = (24., 25.5)


    spectrumF4.K = 3.0e-15 * fluxUnit
    spectrumF4.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) 
    spectrumF4.K.fix = True

    spectrumF4.index = -2.42
    spectrumF4.index.bounds = (-4., 0.0)
    spectrumF4.index.fix = True

    spectrumF4.piv = 7.0 * u.TeV
    spectrumF4.piv.fix = True
# ----------------------------------------

    spectrumF5 = Powerlaw()
    sourceF5  = PointSource("F_j2003",
                               300.89 ,
                               32.80,
                               spectral_shape=spectrumF5)
    sourceF5.position.ra.free = False
    #sourceF5.position.ra.bounds = (291.5,293.5)
    
    sourceF5.position.dec.free = False
    #sourceF5.position.dec.bounds = (18., 19.5)


    spectrumF5.K = 4.6e-15 * fluxUnit
    spectrumF5.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) 
    spectrumF5.K.fix = True

    spectrumF5.index = -2.64
    spectrumF5.index.bounds = (-4., 0.0)
    spectrumF5.index.fix =True

    spectrumF5.piv = 7.0 * u.TeV
    spectrumF5.piv.fix = True


    spectrumF6_1 = Powerlaw()
    sourceF6_1 = PointSource("F_j2006",
                               301.69 ,
                               33.96,
                               spectral_shape=spectrumF6_1)
    sourceF6_1.position.ra.free = False
    #sourceF5.position.ra.bounds = (291.5,293.5)
    
    sourceF6_1.position.dec.free = False
    #sourceF5.position.dec.bounds = (18., 19.5)


    spectrumF6_1.K = 6.6e-15 * fluxUnit
    spectrumF6_1.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) 
    spectrumF6_1.K.fix = True

    spectrumF6_1.index = -2.53
    spectrumF6_1.index.bounds = (-4., 0.0)
    spectrumF6_1.index.fix =True

    spectrumF6_1.piv = 7.0 * u.TeV
    spectrumF6_1.piv.fix = True


    shapeF = SpatialTemplate_2D()
    shapeF.load_file('F_2D_pi0_DRAGON_BASE_v4_Amid.fits') #DRAGON, pi0
    
    spectrumF = Powerlaw()
    spectrumF.K = 1.53e-18 # * fluxUnit
    spectrumF.K.fix = True
    # spectrumFDE4.K.bounds = (1e-11, 1e-4)
    
    spectrumF.index = -2.45
    spectrumF.index.fix = True
    # spectrumFDE4.index.bounds = (-4, -1)
  
    spectrumF.piv = 7.0e9 #* u.TeV
    spectrumF.piv.fix = True
    
    shapeF.hash = 1.
    shapeF.hash.fix = True
    
    shapeF.K = 1.
    shapeF.K.fix = True
    
    sourceF = ExtendedSource("F_pi0", spatial_shape=shapeF, spectral_shape = spectrumF )
    

    


    spectrumG1 = Powerlaw()
    shapeG1 = Gaussian_on_sphere()
    sourceG1= ExtendedSource("G_j2019",
                               spatial_shape=shapeG1,
                               spectral_shape=spectrumG1)
    shapeG1.lon0 = 304.85  * u.degree
    shapeG1.lon0.fix = True
    shapeG1.lat0 = 36.80 * u.degree
    shapeG1.lat0.fix = True
                               
    shapeG1.sigma = 0.302 * u.degree
    shapeG1.sigma.fix = True
    shapeG1.sigma.max_value = 2.
                               
    spectrumG1.K = 5.64e-14 * fluxUnit
    spectrumG1.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrumG1.K.fix = True
                               
    spectrumG1.piv = 7.0 * u.TeV
    spectrumG1.piv.fix = True
                               
    spectrumG1.index = -2.256
    spectrumG1.index.bounds = (-4., 0.0)
    spectrumG1.index.fix = True


    spectrumG2 = Powerlaw()
    shapeG2 = Gaussian_on_sphere()
    sourceG2= ExtendedSource("G_j2020",
                               spatial_shape=shapeG2,
                               spectral_shape=spectrumG2)
    shapeG2.lon0 = 305.16  * u.degree
    shapeG2.lon0.fix = True
    shapeG2.lat0 = 40.37 * u.degree
    shapeG2.lat0.fix = True
                               
    shapeG2.sigma = 0.37 * u.degree
    shapeG2.sigma.fix = True
    shapeG2.sigma.max_value = 2.
                               
    spectrumG2.K = 2.6e-14 * fluxUnit
    spectrumG2.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrumG2.K.fix = True
                               
    spectrumG2.piv = 7.0 * u.TeV
    spectrumG2.piv.fix = True
                               
    spectrumG2.index = -2.81
    spectrumG2.index.bounds = (-4., 0.0)
    spectrumG2.index.fix = True



# ------------------------------------------------

    spectrumG3 = Powerlaw() 
    sourceG3   = PointSource("G_J2031_435",
                                307.99,
                                43.56,
                               spectral_shape=spectrumG3)
                               
    sourceG3.position.ra.free = False
  #  sourceG3.position.ra.bounds = (275.5,277.5)
    
    sourceG3.position.dec.free = False
  #  sourceG3.position.dec.bounds = (-14.0, -12.0)
    
    spectrumG3.K = 5.1e-15 * fluxUnit
    spectrumG3.K.fix = True

    spectrumG3.index =  -2.79
    spectrumG3.index.fix = True

    spectrumG3.piv = 7.0 * u.TeV
    spectrumG3.piv.fix = True


# ----------------------------------------------


    spectrumG4 = Powerlaw()
    shapeG4= Gaussian_on_sphere()
    sourceG4= ExtendedSource("G_j2031_415",
                               spatial_shape=shapeG4,
                               spectral_shape=spectrumG4)
    shapeG4.lon0 = 307.93  * u.degree
    shapeG4.lon0.fix = True
    shapeG4.lat0 = 41.51 * u.degree
    shapeG4.lat0.fix = True
                               
    shapeG4.sigma = 0.298 * u.degree
    shapeG4.sigma.fix = True
    shapeG4.sigma.max_value = 2.
                               
    spectrumG4.K = 4.7e-14 * fluxUnit
    spectrumG4.K.bounds = (1e-18 * fluxUnit, 1e-10*fluxUnit) # Factor 1e9 to what's output
    spectrumG4.K.fix = True
                               
    spectrumG4.piv = 7.0 * u.TeV
    spectrumG4.piv.fix = True
                               
    spectrumG4.index = -2.42
    spectrumG4.index.bounds = (-4., 0.0)
    spectrumG4.index.fix = True




    shapeG = SpatialTemplate_2D()
    shapeG.load_file('G_2D_pi0_DRAGON_BASE_v4_Amid.fits') #DRAGON, pi0
    
    spectrumG = Powerlaw()
    spectrumG.K = 3.08e-18 # * fluxUnit
    spectrumG.K.fix = True
    # spectrumFDE4.K.bounds = (1e-11, 1e-4)
    
    spectrumG.index = -2.53
    spectrumG.index.fix = True 
    # spectrumFDE4.index.bounds = (-4, -1)
  
    spectrumG.piv = 7.0e9 #* u.TeV
    spectrumG.piv.fix = True
    
    shapeG.hash = 1.
    shapeG.hash.fix = True
    
    shapeG.K = 1.
    shapeG.K.fix = True
    
    sourceG = ExtendedSource("G_pi0", spatial_shape=shapeG, spectral_shape = spectrumG )
    




    #==============================================================
    # Set up a likelihood model using the source
   
 
 #  lm = Model(source1, source2,source3,source4,source5,source6, sourceA, sourceB1, sourceB2, source9, source10, source11,  source12, source13, source14, sourceB, sourceC,  sourceC5 ,sourceC2,sourceC3,sourceC4,   sourceD1, sourceD2, sourceD3, sourceD4, sourceD5, sourceD6 ,sourceD, sourceE1, sourceE2,sourceE3_1, sourceE, sourceF1, sourceF2,sourceF3,sourceF4,sourceF5,sourceF6_1 ,sourceF,sourceG1,sourceG2,sourceG3, sourceG4 ,sourceG)
 
 # lm = Model(source1, source2,source3,source4,source5,source6, sourceB1, sourceB2, source9, source10, source11,  source12, source13, source14,  sourceC5 ,sourceC2,sourceC3,sourceC4,   sourceD1, sourceD2, sourceD3, sourceD4, sourceD5, sourceD6 , sourceE1, sourceE2,sourceE3_1, sourceE, sourceF1, sourceF2,sourceF3,sourceF4,sourceF5,sourceF6_1 ,sourceG1,sourceG2,sourceG3, sourceG4 )
 
    lm = Model( sourceA, sourceB, sourceC,  sourceD, sourceE ,sourceF, sourceG)
  
    llh = HAWCLike("Extended", mtfile, rsfile,fullsky=True) #if fullsky you need to define a ROI
    llh.set_active_measurements(args.startBin, args.stopBin)
    #SET ROI
    llh.set_strip_ROI(8 , 86 , -6, 6, fixed_ROI=False, galactic=True)
            

    # Double check the free parameters
    print("Likelihood model:\n")
    lm.display(complete=True)

    # Set up the likelihood and run the fit
    print("Performing likelihood fit...\n")
    datalist = DataList(llh)
    jl = JointLikelihood(lm, datalist, verbose=True)

  #  try:
   #     jl.set_minimizer("ROOT")
    #    result = jl.fit()
  #  except AttributeError:
   #     jl.set_minimizer("minuit")
    #    result = jl.fit()


    if NF==1:
        try:
            jl.set_minimizer("ROOT")
        except AttributeError:
            jl.set_minimizer("minuit")

    if NF==0:
        try:
            jl.set_minimizer("ROOT")
            result = jl.fit()
        except AttributeError:
            jl.set_minimizer("minuit")
            result = jl.fit()


    TS = llh.calc_TS()
    sigma = np.sqrt(TS)


    print("Test statistic: %g" % TS)
    print("Significance:   %g\n" % sigma)
    llh.write_residual_map(output+'_residual.root')
    llh.write_model_map(output+'_model.root')

if __name__=="__main__":

    import argparse

    p = argparse.ArgumentParser(description="Example spectral fit with LiFF")
    
    p.add_argument("-n", "--ntransit", dest="ntransit", default=1, type=float,
                   help="Number of transits to use")
    p.add_argument("--startBin", dest="startBin", default=1, type=int,
                   help="Starting analysis bin [0..9]")
    p.add_argument("--stopBin", dest="stopBin", default=9, type=int,
                   help="Stopping analysis bin [0..9]")


    args = p.parse_args()

    go(args)
