# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 17:18:05 2016

@author: stephaniekwan

Quick zero-point magnitude calculations for JHK bands, using 2MASS sources near
IC 10 X-2. See associated .fits files and 160802 H band 2MASS sources.reg. +
Not very quick photometry on WIRC images in the JHK bands, calculating appmag
from counts and then flux from appmag.

Circle = 3 arcsec radius, Annulus = 5 arcsec, 10 arcsec radii.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import *  # Orthogonal Distance Regression
from scipy.stats import norm
randNorm = np.random.normal
#==============================================================================
# H band
#==============================================================================
hMag, hMagSig, hCircMean, hCircNpx, hAnnMean, hAnnStd = np.split(
            np.array([[13.870, 0.042, 59.184453, 451, 2.3665121, 3.9559723], 
                     [13.071, 0.030, 129.92435, 453, 14.036781, 6.3190906],
                     [11.902, 0.031, 342.62851, 453, 30.107537, 9.0172269],
                     [12.995, 0.041, 143.09626, 453, 25.565892, 7.0594376],
                     [13.044, 0.035, 124.8055,  454, 10.677898, 4.8077078],
                     [13.149, 0.038, 127.78717, 453, 11.363273, 8.2873377],
                     [12.934, 0.028, 131.01363, 453, 1.5582469, 3.6310875]]),
                     6, # number of groups to split into along y-axis
                     axis = 1)
                     
# Equivalent data for IC 10 X-2 circle and annulus.
hX2CircMean, hX2CircNpx, hX2AnnMean, hX2AnnStd = \
     np.array([41.788069, 312, 28.881225, 10.546441])
    
# Calculate number of non-background counts in the circle (MJy/sr)
hCircCount     = (hCircMean - hAnnMean) * hCircNpx
hX2CircCount   = (hX2CircMean - hX2AnnMean) * hX2CircNpx
# Scaling factors: (1) WIRC pixel to arcsecond squared, (2) sr to arcsec
scale, srOverArcSec = 0.2487 ** 2, 1/(4.25 * 10**10)
# Convert counts per pixel to counts per steradian
hCircDat   = hCircCount * scale * srOverArcSec  * 10**9
hX2CircDat = hX2CircCount * scale * srOverArcSec * 10**9
# Estimate error of the counts based on annulus std. dev
hCircDatErr = hAnnStd * np.sqrt(hCircNpx) * scale * srOverArcSec  * 10**9
hX2CircDatErr = hX2AnnStd * np.sqrt(hX2CircNpx) * scale * srOverArcSec *10**9

## Use scipy orthogonal fit
def magCounts(param, logCounts):
    '''A linear function: apparent magnitude = zeropoint - 2.5 * log10(counts).
    param is a one-element array containing zeropoint.'''
    return param[0] - 2.5 * logCounts
hLinear = Model(magCounts)             # Create a model object
hData = Data(np.log10(hCircDat),      # x value
             hMag,                    # y value
             wd = np.log10(hCircDatErr), # standard deviation of x
             we = hMagSig )           # standard deviation of y
             
hOdr = ODR(hData, hLinear, beta0 = [15])   
hOutput = hOdr.run()
hOutput.pprint()

hZeroPt, hZeroPtsig = 17.76510759, 0.05394425   # read off of console output
# Now use our relation between magnitude and data points
hX2Mag = hZeroPt - 2.5 * np.log10(hX2CircDat)

#==============================================================================
# H band bootstrap to estimate uncertainty in H mag
#==============================================================================
iters = 1000000
hX2Mag_rand, hX2Flux_rand = np.zeros((2, iters))
for i in range(iters):
    hX2Mag_rand[i] = randNorm(hZeroPt, hZeroPtsig) - 2.5 * \
                     np.log10(randNorm(hX2CircDat, hX2CircDatErr))
hX2Mag_mu, hX2Magsig = norm.fit(hX2Mag_rand)

# Now convert H-band mag to an actual flux, using F_nu 0 mag = 1024 +/- 20.0 Jy
hFnu0, hFnu0sig = 1024., 20. # in Jy
hX2Flux = hFnu0 * 10 ** (- hX2Mag / 2.5) * 10**3
#==============================================================================
# H band bootstrap to estimate uncertainty in H band flux
#==============================================================================
for i in range(iters):              
    hX2Flux_rand[i] = randNorm(hFnu0, hFnu0sig) * 10 ** 3 * \
                            10 ** (-randNorm(hX2Mag, hX2Magsig) / 2.5) 
hX2Flux_mu, hX2Flux_sig = norm.fit(hX2Flux_rand)   # Extract values

print 'H BAND: zero-flux mag of ' + str(hZeroPt) + ' +/-' + str(hZeroPtsig) + \
 ' mags. For IC 10 X-2, flux: ' + str(hX2Flux) + ' +/- ' + str(hX2Flux_sig) + \
 ' mJy, and mag: ' + str(hX2Mag) + ' +/- ' + str(hX2Magsig) + ' mags.' 

#==============================================================================
# K band  (since I have no concept of doing things efficiently)
#==============================================================================
kMag, kMagSig, kCircMean, kCircNpx, kAnnMean, kAnnStd = np.split(
            np.array([[13.778, 0.0044, 38.924658, 455, 0.14424044, 4.7129194], 
                     [12.989, 0.028, 83.745914, 451, 11.262791, 6.8961484],
                     [11.835, 0.029, 239.9615, 453, 24.275672, 10.870496],
                     [12.688, 0.029, 125.76851, 451, 9.2223908, 8.3131534],
                     [12.949, 0.029, 104.21588, 452, 19.995094, 7.6435669],
                     [12.823, 0.034, 96.403545, 453, 7.1562142, 5.5116115],
                     [12.862, 0.032, 87.336419, 454, 0.96671014, 4.2705642]]),
                     6, # number of groups to split into along y-axis
                     axis = 1)
# Equivalent data for IC 10 X-2 circle and annulus.
kX2CircMean, kX2CircNpx, kX2AnnMean, kX2AnnStd = \
    np.array([29.753703, 454, 22.490616, 26.557954])
    
# Calculate number of non-background counts in the circle (MJy/sr)
kCircCount     = (kCircMean - kAnnMean) * kCircNpx
kX2CircCount   = (kX2CircMean - kX2AnnMean) * kX2CircNpx
# Convert data points from counts/sr to just counts
kCircDat   = kCircCount * scale * srOverArcSec * 10**9
kX2CircDat = kX2CircCount * scale * srOverArcSec * 10**9
# Estimate error in the counts based on annulus std. dev
kCircDatErr   = kAnnStd * np.sqrt(kCircNpx) * scale * srOverArcSec * 10**9
kX2CircDatErr = kX2AnnStd * np.sqrt(kX2CircNpx) * scale * srOverArcSec * 10**9

# Use scipy orthogonal fit
kLinear = Model(magCounts)             # Create a model object
kData = Data(np.log10(kCircDat),     # x value
             kMag,                    # y value
             wd = np.log10(kCircDatErr), # standard deviation of x
             we = kMagSig)           # standard deviation of y
             
kOdr = ODR(kData, kLinear, beta0 = [15.] )   
kOutput = kOdr.run()
kOutput.pprint()

kZeroPt, kZeroPtsig    = 17.268, 0.06733548      # read off of console output
kX2Mag    = kZeroPt - 2.5 * np.log10(kX2CircDat)  # Calculate K mag

#==============================================================================
# K band bootstrap to estimate uncertainty in K mag
#==============================================================================
kX2Mag_rand, kX2Flux_rand = np.zeros((2, iters))
for i in range(iters):
    kX2Mag_rand[i] = randNorm(kZeroPt, kZeroPtsig) - 2.5 * \
                     np.log10(randNorm(kX2CircDat, kX2CircDatErr))
kX2Mag_mu, kX2Magsig = norm.fit(kX2Mag_rand)

# Now convert K-band mag to an actual flux, using 2MASS 0-mag Jy fluxes
kFnu0, kFnu0sig = 666.7, 12.6 # in Jy
kX2Flux = kFnu0 * 10 ** (- kX2Mag / 2.5) * 10**3
#==============================================================================
# K band bootstrap to estimate uncertainty in K flux
#==============================================================================
for i in range(iters):              
    kX2Flux_rand[i] = randNorm(kFnu0, kFnu0sig) * 10 ** 3 * \
                            10 ** (-randNorm(kX2Mag, kX2Magsig) / 2.5) 
kX2Flux_mu, kX2Flux_sig = norm.fit(kX2Flux_rand)   # Extract values

print 'K BAND: zero-flux mag of ' + str(kZeroPt) + ' +/-' + str(kZeroPtsig) + \
 ' mags. For IC 10 X-2, flux: ' + str(kX2Flux) + ' +/- ' + str(kX2Flux_sig) + \
 ' mJy, and mag: ' + str(kX2Mag) + ' +/- ' + str(kX2Magsig) + ' mags.' 
 
#==============================================================================
# LAST ONE!: J band  
#==============================================================================  
jMag, jMagSig, jCircMean, jCircNpx, jAnnMean, jAnnStd = np.split(
            np.array([[14.239, 0.027, 34.567618, 453, 0.28655016, 1.3666131],
                      [13.153, 0.021, 93.13718, 454, 4.9306852, 2.3523632],
                      [12.398, 0.021, 180.97938, 454, 10.741922, 3.0262562],
                      [13.422, 0.027, 80.871289, 450, 10.303076, 2.6931108],
                      [13.858, 0.026, 53.657189, 450, 4.1546406, 1.7615609],
                      [14.201, 0.029, 42.78187,  455, 3.8618451, 3.3292666],
                      [13.134, 0.021, 93.067053, 446, 0.61044178, 1.3005565]]),
                     6, # number of groups to split into along y-axis
                     axis = 1)
# Equivalent data for IC 10 X-2 circle and annulus.
jX2CircMean, jX2CircNpx, jX2AnnMean, jX2AnnStd = \
    np.array([18.11532, 316, 11.934676, 4.1923764])
    
# Calculate number of non-background counts in the circle (MJy/sr)
jCircCount     = (jCircMean - jAnnMean) * jCircNpx
jX2CircCount   = (jX2CircMean - jX2AnnMean) * jX2CircNpx
# Convert "counts" (counts/sr) to counts
jCircDat   = jCircCount * scale * srOverArcSec * 10**9
jX2CircDat = jX2CircCount * scale * srOverArcSec * 10**9
# Estimate error in the circle in MJy based on annulus std. dev
jCircDatErr   = jAnnStd * np.sqrt(jCircNpx) * scale * srOverArcSec * 10**9 
jX2CircDatErr = jX2AnnStd * np.sqrt(jX2CircNpx) * scale * srOverArcSec * 10**9

# Use scipy orthogonal fit
jLinear = Model(magCounts)             # Create a model object
jData = Data(np.log10(jCircDat),     # x value
             jMag,                    # y value
             wd = np.log10(jCircDatErr), # standard deviation of x
             we = jMagSig )           # standard deviation of y
             
jOdr = ODR(jData, jLinear, beta0 = [17])   
jOutput = jOdr.run()
jOutput.pprint()

jZeroPt, jZeroPtsig = 17.612, 0.06163881      # read off of console output
jX2Mag = jZeroPt - 2.5 * np.log10(jX2CircDat)  # Calculate J mag

#==============================================================================
# J band bootstrap to estimate uncertainty in J mag
#==============================================================================
jX2Mag_rand, jX2Flux_rand = np.zeros((2, iters))
for i in range(iters):
    jX2Mag_rand[i] = randNorm(jZeroPt, jZeroPtsig) - 2.5 * \
                     np.log10(randNorm(jX2CircDat, jX2CircDatErr))
jX2Mag_mu, jX2Magsig = norm.fit(jX2Mag_rand)

# Now convert J-band mag to an actual flux, using 2MASS 0-mag Jy fluxes
jFnu0, jFnu0sig = 1594., 27.8 # in Jy
jX2Flux = jFnu0 * 10 ** (- jX2Mag / 2.5) * 10**3
#==============================================================================
# J band bootstrap to estimate uncertainty in J flux
#==============================================================================
for i in range(iters):              
    jX2Flux_rand[i] = randNorm(jFnu0, jFnu0sig) * 10 ** 3 * \
                            10 ** (-randNorm(jX2Mag, jX2Magsig) / 2.5) 
jX2Flux_mu, jX2Flux_sig = norm.fit(jX2Flux_rand)   # Extract values

print 'J BAND: zero-flux mag of ' + str(jZeroPt) + ' +/-' + str(jZeroPtsig) + \
 ' mags. For IC 10 X-2, flux: ' + str(jX2Flux) + ' +/- ' + str(jX2Flux_sig) + \
 ' mJy, and mag: ' + str(jX2Mag) + ' +/- ' + str(jX2Magsig) + ' mags.' 
 