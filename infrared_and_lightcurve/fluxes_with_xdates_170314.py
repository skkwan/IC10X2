# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 17:27:53 2016
Updated 2016 Aug 4 to include WIRC points
Updated 2016 Sep 7 to reformat

@author: stephaniekwan
IC-10 X-2 Spitzer IRAC light curves. Updated on August 4th to include July 18th
Palomar JHK band measurements.
Coordinates: 0:20:20.940 +59:17:59.00
Circle radius: 3 arcsec. Annulus radii: 5 arcsec inner, 10 arcsec outer
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from jdcal import gcal2jd
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
# Use LaTeX font
plt.rc({'weight' : 'normal',
        'size'   : 15})
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('text', usetex = True)

index = np.array(['4424960 2004-07-23', '33204224 2010-01-29', 
                  '33203968 2010-02-19', '33203456 2010-03-10', 
              '33202944 2010-09-09' '33202432 2010-10-04', '33201920 2010-10-14',
              '42321152 2011-09-24', '42321408 2012-04-04', '52576256 2015-03-11', 
              '52577280 2016-03-23', '2016-07-18'])
gdates = np.array([(2004,07,23), (2010,01,29), (2010,02,19), (2010,03,10),
                   (2010,9,9), (2010,10,4), (2010,10,14), (2011,9,24),
                   (2012,04,04), (2015,03,11), (2016,03,23)])

# Convert to MJD dates
jdates = np.zeros(len(gdates))
for i in range(len(gdates)):
    jdates[i] = gcal2jd( gdates[i][0], gdates[i][1], gdates[i][2] )[1]
WIRCjdate = gcal2jd(2016,7,18)[1]

#X-ray observation dates: 2003 Mar 12, 2006 Nov 5, 2010 May 21
xdates = np.array([52710.7, 54044.2, 55337.8])
# Dates of non-observations
xnondates = np.array([55140.7, 55190.2, 55238.5, 55290.6, 55397.5, 55444.6])

# Mean counts in circle (units: counts/sr per pixel), for 4.5 and 3.6 microns
# Mean counts in annulus: in units of MJr/sr per pixel
circMean36 = np.array([2.235913,
1.9806753,
1.8627226,
1.9333704,
1.9426107,
1.9242988,
1.8619019,
1.8695578,
1.9416175,
1.8303715,
1.8961317])
                       
annMean36 = np.array([1.502455, 1.4441012, 1.4349068, 1.4300396, 
                      1.4522621, 1.4369512, 1.4367747,
                      1.4509853, 1.4649935, 1.4423924, 1.4682426])    
annSD36 = np.array([0.323036, 0.33284634, 0.30873036, 0.27726872,
                    0.30360679, 0.29375085, 0.31357359,
                    0.32412101, 0.30720197, 0.28204827, 0.28241972])                      
circMean45 = np.array([1.6294469, 1.3514017, 1.2583814, 1.2950296,
                       1.3489466, 1.2898556, 1.2250279,
                     1.2813393, 1.343888, 1.2231404, 1.2529148])
                     
annMean45 = np.array([1.0128354, 0.93392948, 0.94994089, 0.96776315,
                      0.98786045, 0.93146131,0.91232822, 0.96418034,
                      1.0059549, 0.93307992, 0.94233364])  
annSD45 = np.array([0.18814292, 0.19965652, 0.19302296,  0.18062225, 
                    0.18524006, 0.18025225, 0.18849567, 0.19213017,
                    0.18247341, 0.19707077, 0.20098456])                      
circMean58 = np.array([2.4857705])  #only for '4424960 2004-07-23'
circMean80 = np.array([5.6362584])  # " "

annMean58 = np.array([2.2773678])
annMean80 = np.array([5.8670916])

# Standard deviation in annulus counts (counts/sr per pixel)
annSD58 = np.array([0.34377934])
annSD80 = np.array([0.81536177])

# Number of pixels in circle
circNpx36 = np.array([54,52,54,55,52,54,55,56,53,56,55])
circNpx45 = np.array([54,52,54,55,52,54,55,56,53,56,55])
circNpx58, circNpx80 = np.array([54]), np.array([54])

# Calculate number of non-background counts in the circle (counts/sr)
circCounts36 = (circMean36 - annMean36) * circNpx36
circCounts45 = (circMean45 - annMean45) * circNpx45
circCounts58 = (circMean58 - annMean58) * circNpx58
circCounts80 = (circMean80 - annMean80) * circNpx80

# Conversion between steradians and arcsecond. 1 steradian is 4.25 *^ 10 arcseconds
srOverArcSec = 1/(4.25 * 10**10)

# 1 pixel has 0.3600 arcsec^2. Convert "counts" (counts/sr) to counts
circFlux36 = circCounts36 * 0.3600 * srOverArcSec * 10**9
circFlux45 = circCounts45 * 0.3600 * srOverArcSec * 10**9
circFlux58 = circCounts58 * 0.3600 * srOverArcSec * 10**9
circFlux80 = circCounts80 * 0.3600 * srOverArcSec * 10**9

# Estimation of error: standard dev. in annulus counts times area of circle
fluxError36 = annSD36 * np.sqrt(circNpx36) * srOverArcSec * 10**9 * 0.3600
fluxError45 = annSD45 * np.sqrt(circNpx45) * srOverArcSec * 10**9 * 0.3600
fluxError58 = annSD58 * np.sqrt(circNpx58) * srOverArcSec * 10**9 * 0.3600
fluxError80 = annSD80 * np.sqrt(circNpx80) * srOverArcSec * 10**9 * 0.3600

# Palomar P200 JHK fluxes and errors (in mJy)
jFlux2, jErr2 = 0.3822, 0.05623
hFlux2, hErr2 = 0.34596, 0.02698
kFlux2, kErr2 =  0.396159, 0.0773288

circFlux58, circFlux80 = np.array([0.21036669]), np.array([0.19616618])
fluxError58, fluxError80 = np.array([0.03456009]), np.array([0.03161511])


# 2MASS fluxes upper limits (in mJy)
jFlux, hFlux, kFlux = 0.4192, 0.7084, 0.4207
jFluxErr = 0.0593
upperLimDate = gcal2jd(2000,9,16)[1]

dates2 = np.array([upperLimDate, WIRCjdate])

## Plot light curves
#fig, ax = plt.subplots()
#plt.hold(True)
#plt.scatter(dates2, jFluxes, facecolors = 'none', marker = '<', s = 30,
#            edgecolors = 'navy')
#plt.scatter(dates2, hFluxes, facecolors = 'none', marker = 's', s = 30,
#             edgecolors = 'royalblue')
#plt.scatter(dates2, kFluxes, facecolors = 'none', marker = '>', s = 30,
#            edgecolors = 'lightskyblue')
#
#plt.scatter(jdates, circFlux36, color = 'black', marker = 'o', s = 15, 
#            zorder = 3)
#plt.scatter(jdates, circFlux45, color = 'grey', marker = 'v', s = 15,
#            zorder = 3)
#plt.scatter(jdates[0], circFlux58, facecolors = 'none', edgecolors =
#             'darkgrey', marker = 'D', s = 22, zorder = 3)
#plt.scatter(jdates[0], circFlux80, facecolors ='none', edgecolors = 'black',
#            marker = 'o', s = 25, zorder = 3)
#plt.xlim([51500,59500])
#plt.ylim([0.00,0.80])
#plt.legend(('J', 'H', 'K$_s$','[3.6]', '[4.5]', '[5.8]', '[8.0]'),
#            scatterpoints = 1,
#            loc = 'upper right',
#            title = 'Filter/Channel',
#            fontsize = 13,
#            frameon = False)
## Plot time of burst and label it
#plt.axvline(x = 55337.8, color = 'k', ls = 'dashed')
##plt.text(55500, 0.45, "2010 May outburst", rotation=90, fontsize=13)
## Plot error bars
#plt.errorbar(WIRCjdate, kFlux, kErr, color = 'lightskyblue')
#plt.errorbar(WIRCjdate, hFlux, hErr, color = 'royalblue')
#plt.errorbar(WIRCjdate, jFlux, jErr, color = 'navy')
#plt.errorbar(dates2[0], j2Flux, 0.0593, color = 'navy')
#plt.errorbar(jdates, circFlux36, yerr = fluxError36, linestyle = 'None',
#             color = 'black', zorder = 2)
#plt.errorbar(jdates, circFlux45, yerr = fluxError45, linestyle = 'None',
#             color = 'grey', zorder = 2)
#plt.errorbar(jdates[0], circFlux58, yerr = fluxError58, linestyle = 'None',
#             color = 'darkgrey', zorder = 2)
#plt.errorbar(jdates[0], circFlux80, yerr = fluxError80, linestyle = 'None',
#             color = 'black', zorder = 2)
#plt.xlabel('Time (MJD)', fontsize = 14)
#plt.ylabel('Flux density (mJy)', fontsize = 14)
#ax.arrow(dates2[0], kFluxes[0], 0.0, -0.08, head_width = 150,
#         head_length = 0.02, fc ='lightskyblue', ec ='lightskyblue')
#ax.arrow(dates2[0], hFluxes[0], 0.0, -0.08, head_width = 150,
#         head_length = 0.02, fc = 'royalblue', ec ='royalblue')
#ax.arrow(WIRCjdate, 0.7, 0.0, -0.15, head_width = 300, head_length = 0.03,
#         fc = 'k', ec = 'k', linestyle = '-')
#for i in range(len(xdates)):     
#       ax.arrow(xdates[i], 0.03, 0.0, 0.04, head_width = 100, head_length = 0.02,
#         fc = 'darkslategrey', ec = 'darkslategrey', linestyle = '-')     
#for j in range(len(xnondates)):
#       ax.arrow(xnondates[j], 0.03, 0.0, 0.02, head_width = 30, head_length = 0.015,
#        fc = 'lightgrey', ec = 'lightgrey', linestyle = '-')

#x1, x2, y1, y2 = 55100, 56100, 0.01, 0.27
## Plot quiescent levels
#quies36 = 0.19897806
#quies45 = 0.146110673
#plt.axhline(y = quies36, color = 'black', ls = 'dashed')
#plt.axhline(y = quies45, color = 'grey', ls = 'dashed')
#plt.text(51650, 0.11, 'Quiescent levels', fontsize = 10)
#
## Shaded area to denote uncertainty of median (average of mag1sig)
#ax.add_patch(                     
#    patches.Rectangle(
#        (51500, quies36 - np.average(fluxError36)),  # (x, y)
#        59500 - 51500,   # width 
#        2 * np.average(fluxError36),   # height
#        0.0,                       # angle
#        facecolor = 'gainsboro',
#        edgecolor = 'none',
#        zorder = 1
#        ))   
#ax.add_patch(                     
#    patches.Rectangle(
#        (51500, quies45 - np.average(fluxError45)),  # (x, y)
#        59500 - 51500,   # width 
#        2 * np.average(fluxError45),   # height
#        0.0,                       # angle
#        facecolor = 'gainsboro',
#        edgecolor = 'none',
#        zorder = 1
#        ))   
#        
######
### Zoomed inset
#####
#axins = zoomed_inset_axes(ax,1.8,loc=9)
#axins.set_xlim(x1,x2)
#axins.set_ylim(y1,y2)
#plt.scatter(jdates[1:9], circFlux36[1:9], color = 'black', marker = 'o',
#            s = 15, zorder = 2)
#plt.errorbar(jdates[1:9], circFlux36[1:9], yerr = fluxError36[1:9],
#             linestyle = 'None', color = 'black', zorder = 2)
#plt.scatter(jdates[1:9], circFlux45[1:9], color = 'dimgrey', marker = 'v', s = 15,
#            zorder = 3)
#plt.errorbar(jdates[1:9], circFlux45[1:9], yerr = fluxError45[1:9],
#             linestyle = 'None', color = 'dimgrey', zorder = 3)
#plt.axvline(x = 55337.8, color = 'k', ls = 'dashed')
## Plot quiescent levels
#plt.axhline(y = quies36, color = 'k', ls = 'dashed')
#plt.axhline(y = quies45, color = 'grey', ls = 'dashed')
#
#plt.xticks(np.arange(x1, x2, 400))
#mark_inset(ax, axins, loc1 = 3, loc2 = 4, fc = "none", ec = "0.6")
#for i in range(len(xdates)):     
#       axins.arrow(xdates[i], 0.03, 0.0, 0.04, head_width = 100, head_length = 0.02,
#         fc = 'darkslategrey', ec = 'darkslategrey', linestyle = '-')     
#for j in range(len(xnondates)):
#       axins.arrow(xnondates[j], 0.03, 0.0, 0.02, head_width = 30, head_length = 0.015,
#        fc = 'lightgrey', ec = 'lightgrey', linestyle = '-')
#
## Shaded area to denote uncertainty of median (average of mag1sig)
#axins.add_patch(                     
#    patches.Rectangle(
#        (x1, quies36 - np.average(fluxError36)),  # (x, y)
#        x2 - x1,   # width 
#        2 * np.average(fluxError36),   # height
#        0.0,                       # angle
#        facecolor = 'gainsboro',
#        edgecolor = 'none',
#        zorder = 1
#        ))   
#axins.add_patch(                     
#    patches.Rectangle(
#        (x1, quies45 - np.average(fluxError45)),  # (x, y)
#        x2 - x1,   # width 
#        2 * np.average(fluxError45),   # height
#        0.0,                       # angle
#        facecolor = 'gainsboro',
#        edgecolor = 'none',
#        zorder = 1
#        ))   

#fig.savefig("170501_IC10X2_light_curve.pdf")

##################
##### Convert fluxes to magnitudes
##################
m36, m36sig, m45, m45sig = np.array([]), np.array([]), np.array([]), np.array([])
assert(len(circFlux36) == len(circFlux45))
for i in range(0, len(circFlux36)):
    m36 = np.append(m36, -2.5 * np.log10(circFlux36[i] / (280.9 * 10**3)) )
    m45 = np.append(m45, -2.5 * np.log10(circFlux45[i] / (179.7 * 10**3)) )
    m36sig = np.append(m36sig, m36[i] + float(2.5 *
            np.log10((circFlux36[i] + fluxError36[i]) / (280.9 * 10**3)) ))
    m45sig = np.append(m45sig, m45[i] + float(2.5 *
            np.log10((circFlux45[i] + fluxError45[i]) / (179.7 * 10**3))    ))
m58 = -2.5 * np.log10(circFlux58 / (115.0 * 10**3))  
m58sig = m58 + float(2.5 * np.log10((circFlux58 + fluxError58) / 
                                (115.0 * 10**3))  )
m80 = -2.5 * np.log10(circFlux80 / (64.9 * 10**3))  
m80sig = m80 + float( 2.5 * np.log10((circFlux80 + fluxError80) / 
                                (64.9 * 10**3)) )


# Zero-magnitude fluxes in Jy for 2MASS (as posted on website)
jzeroMagFlux, jzeroMagFluxsig = 1594, 27.8 
hzeroMagFlux, hzeroMagFluxsig = 1024, 20.0
kzeroMagFlux, kzeroMagFluxsig = 666.7, 12.6
# Calculate absolute magnitude of 2MASS values
jmag1 = -2.5 * np.log10(jFlux * 10**-3/jzeroMagFlux)
jmag1sig = np.sqrt((jzeroMagFluxsig * (2.5/np.log(10)) / jzeroMagFlux )**2 
                    + (jFluxErr * (2.5/np.log(10)) / jFlux )**2)
# Upper limits on H and K bands (no sig)
hmag1 = -2.5 * np.log10(hFlux * 10**-3 / hzeroMagFlux) 
kmag1 = -2.5 * np.log10(kFlux * 10**-3 / kzeroMagFlux) 

# P200 magnitudes: see 160802 WIRC zero points.py file
jmag2, jmag2sig = 16.47779136, 0.0742855400832
hmag2, hmag2sig = 15.8452643514, 0.0738319761354
kmag2, kmag2sig = 15.5651506508, 0.205591601092

fig, ax = plt.subplots()

jMags = np.array([jmag1, jmag2])     # Concatenate data
hMags = np.array([hmag1, hmag2])
kMags = np.array([kmag1, kmag2])
# Plot MAGNITUDE light curves
#plt.scatter(dates2, jMags, facecolors = 'navy', marker = '<', s = 25,
#            edgecolors = 'navy', zorder = 3)
#plt.scatter(dates2, hMags, facecolors = 'none', marker = 's', s = 25,
#             edgecolors = 'royalblue', zorder = 3)
#plt.scatter(dates2, kMags, facecolors = 'none', marker = '>', s = 25,
#            edgecolors = 'lightskyblue', zorder = 3)
## Add error arrows to the H and K 2MASS images
#ax.arrow(dates2[0], kMags[0], 0.0, 0.3, head_width = 100,
#         head_length = 0.15, fc ='lightskyblue', ec ='lightskyblue', zorder = 2)
#ax.arrow(dates2[0], hMags[0], 0.0, 0.3, head_width = 100,
#         head_length = 0.15, fc = 'royalblue', ec ='royalblue', zorder = 1)
#plt.errorbar(dates2[0], jMags[0], yerr = jmag1sig, linestyle = 'None',
#            color = 'navy', zorder = 3) 
## Plot errorbars for WIRC datapoints
#plt.errorbar(dates2[1], jMags[1], yerr = jmag2sig, linestyle = 'None',
#             color = 'navy', zorder = 3) 
#plt.errorbar(dates2[1], hMags[1], yerr = hmag2sig, linestyle = 'None',
#             color = 'royalblue', zorder = 3)   
#plt.errorbar(dates2[1], kMags[1], yerr = kmag2sig, linestyle = 'None',
#             color = 'lightskyblue', zorder = 3)  
## Plot errorbars for Spitzer datapoints
#plt.errorbar(jdates, m36, yerr = m36sig, linestyle = 'None',
#             color = 'k', zorder = 2)
#plt.errorbar(jdates, m45, yerr = m45sig, linestyle = 'None',
#             color = 'grey', zorder = 2)
#plt.errorbar(jdates[0], m58, yerr = m58sig, linestyle = 'None',
#             color = 'grey', zorder = 1)
#plt.errorbar(jdates[0], m80, yerr = m80sig, linestyle = 'None',
#             color = 'grey', zorder = 1)             
## Plot the Spitzer datapoints
#plt.scatter(jdates, m36, color = 'black', marker = 'o', s = 15, 
#            zorder = 3)
#plt.scatter(jdates, m45, color = 'grey', marker = 'v', s = 15,
#            zorder = 3)
#plt.scatter(jdates[0], m58, facecolors = 'none', edgecolors =
#             'grey', marker = 'D', s = 22, zorder = 3)
#plt.scatter(jdates[0], m80, facecolors ='none', edgecolors = 'grey',
#            marker = 'o', s = 25, zorder = 3)
## Set magnitude plot limits           
#plt.xlim([51500, 58100])
#plt.ylim([13.5, 17.0])
### Plot legend
#lgd = plt.legend(('J', 'H', 'K$_s$','[3.6]', '[4.5]', '[5.8]', '[8.0]'),
#            scatterpoints = 1,
#            loc = 'upper right',
#            title = 'Filter/Channel',
#            fontsize = 13,
#            frameon = True,
#            bbox_to_anchor = (1.26, 1.03))
## Plot non-detections
#for i in range(len(xdates)):     
#       ax.arrow(xdates[i], 17.1, 0.0, -0.2, head_width = 100, head_length = 0.1,
#         fc = 'darkslategrey', ec = 'darkslategrey', linestyle = '-')     
#for j in range(len(xnondates)):
#       ax.arrow(xnondates[j], 17.1, 0.0, -0.2, head_width = 30, head_length = 0.1,
#        fc = 'lightgrey', ec = 'lightgrey', linestyle = '-')            
### Plot quiescent magnitudes
#quiesmag36 = np.mean(m36)
#quiesmag45 = np.mean(m45)
#plt.axhline(y = quiesmag36, color = 'black', ls = 'dashed')
#plt.axhline(y = quiesmag45, color = 'grey', ls = 'dashed')
#plt.text(52650, 15.5, 'Quiescent levels', fontsize = 10)
### Shade rectangles around it
#ax.add_patch(                     
#    patches.Rectangle(
#        (51500, quiesmag36 - np.average(m36sig)),  # (x, y)
#        60000 - 51500,   # width 
#        2 * np.average(m36sig),   # height
#        0.0,                       # angle
#        facecolor = 'gainsboro',
#        edgecolor = 'none',
#        zorder = 1,
#        alpha = 0.95
#        ))   
#ax.add_patch(                     
#    patches.Rectangle(
#        (51500, quiesmag45 - np.average(m45sig)),  # (x, y)
#        60000 - 51500,   # width 
#        2 * np.average(m45sig),   # height
#        0.0,                       # angle
#        facecolor = 'lightgrey',
#        edgecolor = 'none',
#        zorder = 2,
#        alpha = 0.4
#        ))   
#
## Plot time of burst and label it
#plt.axvline(x = 55337.8, color = 'k', ls = 'dashed')
#plt.text(55500, 15.6, "2010 May outburst", rotation=90, fontsize=11)
#
## Reverse y axis
#plt.gca().invert_yaxis()
#
#plt.xlabel('Date (MJD)', fontsize = 14)
#plt.ylabel('Magnitude', fontsize = 14)
#fig.savefig("170522_IRmags_lc.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
#
#np.savetxt('170502_Spitzer_mags_TENTATIVE.txt',
#           np.column_stack((jdates, m36, m36sig, m45, m45sig)),
#           header = 'Spitzer IC 10 X-2. Col 1: Date (MJD),' +
#           ' Cols 2-3: [3.6] mags and uncertainties, Cols 3-4: [4.5] mags and uncertainties. For the first epoch, '+
#           '[5.8] mag: %f +/- %f and [8.0] mag: %f +/- %f' % (m58, m58sig, m80, m80sig),
#           fmt = '%10.5f')   
#
#print '\n 2MASS magnitudes (date:', dates2[0], ') \n J band =', jmag1, '+/-', jmag1sig,\
#                         '. Upper limits on H and K bands =', hmag1, \
#                         'and', kmag1, '.', \
#            '\n Palomar magnitudes (Date:', dates2[1], ') \n J =', jmag2, '+/-', jmag2sig,\
#            '\n H =', hmag2, '+/-', hmag2sig, \
#            '\n K =', kmag2, '+/-', kmag2sig
#plt.tight_layout()     
#plt.show()           
