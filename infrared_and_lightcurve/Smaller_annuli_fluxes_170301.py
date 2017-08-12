# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 20:21:21 2017

@author: stephaniekwan
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 21:53:55 2017

@author: stephaniekwan
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 17:27:53 2016
Updated 2016 Aug 4 to include WIRC points
Updated 2016 Sep 7 to reformat

@author: stephaniekwan
IC-10 X-2 Spitzer IRAC light curves.Updated on August 4th to include July 18th
Palomar JHK band measurements.
Coordinates: 0:20:20.940 +59:17:59.00
Circle radius: 3 arcsec. Annulus radii: 5 arcsec inner, 10 arcsec outer
"""

import numpy as np
import matplotlib.pyplot as plt
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
    jdates[i] = gcal2jd(gdates[i][0], gdates[i][1], gdates[i][2])[1]                       
WIRCjdate = gcal2jd(2016,7,18)[1]
# Mean counts in circle (units: counts/sr per pixel), for 4.5 and 3.6 microns
# Mean counts in annulus: in units of MJr/sr per pixel
circMean36 = np.array([2.235913,
1.9806753,
1.8627226,
1.9333704,
1.9806753,
1.9242988,
1.8619019,
1.8695578,
1.9416175,
1.8303715,
1.8961317])
                       
annMean36 = np.array([1.502455, 1.4441012, 1.4349068, 1.4300396, 
                      1.4441012, 1.4369512, 1.4367747,
                      1.4509853, 1.4649935, 1.4423924, 1.4682426])    
annSD36 = np.array([0.323036, 0.33284634, 0.30873036, 0.27726872,
                    0.33284634, 0.29375085, 0.31357359,
                    0.32412101, 0.30720197, 0.28204827, 0.28241972])                      
circMean45 = np.array([1.6294469, 1.3514017, 1.2583814, 1.2950296,
                       1.3514017, 1.2898556, 1.2250279,
                     1.2813393, 1.343888, 1.2231404, 1.2529148])
                     
annMean45 = np.array([1.0128354, 0.93392948, 0.94994089, 0.96776315,
                      0.93392948, 0.93146131,0.91232822, 0.96418034,
                      1.0059549, 0.93307992, 0.94233364])  
annSD45 = np.array([0.18814292, 0.19965652, 0.19302296,  0.18062225, 
                    0.19965652, 0.18025225, 0.18849567, 0.19213017,
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

# JHK fluxes and errors (in mJy)
jFlux, jErr = 0.3822, 0.05623
hFlux, hErr = 0.34596, 0.02698
kFlux, kErr =  0.396159, 0.0773288

circFlux58, circFlux80 = np.array([0.21036669]), np.array([0.19616618])
fluxError58, fluxError80 = np.array([0.03456009]), np.array([0.03161511])

# JHK fluxes and errors (in mJy)
jFlux, jErr = 0.3822, 0.05623
hFlux, hErr = 0.34596, 0.02698
kFlux, kErr =  0.396159, 0.0773288

# 2MASS fluxes upper limits (in mJy)
j2Flux, h2Flux, k2Flux = 0.4192, 0.7084, 0.4207
j2FluxErr = 0.0593
upperLimDate = gcal2jd(2000,9,16)[1]

dates2 = np.array([upperLimDate, WIRCjdate])
jFluxes = np.array([j2Flux, jFlux])
hFluxes = np.array([h2Flux, hFlux])
kFluxes = np.array([k2Flux, kFlux])
# Plot light curves
fig, ax = plt.subplots()
plt.hold(True)
plt.scatter(dates2, jFluxes, facecolors = 'none', marker = '<', s = 30,
            edgecolors = 'navy')
plt.scatter(dates2, hFluxes, facecolors = 'none', marker = 's', s = 30,
             edgecolors = 'royalblue')
plt.scatter(dates2, kFluxes, facecolors = 'none', marker = '>', s = 30,
            edgecolors = 'lightskyblue')

plt.scatter(jdates, circFlux36, color = 'black', marker = 'o', s = 15)
plt.scatter(jdates, circFlux45, color = 'grey', marker = 'v', s = 15)
plt.scatter(jdates[0], circFlux58, facecolors = 'none', edgecolors =
             'darkgrey', marker = 'D', s = 22)
plt.scatter(jdates[0], circFlux80, facecolors ='none', edgecolors = 'black',
            marker = 'o', s = 25)
plt.xlim([51500,59500])
plt.ylim([0.00,0.80])
plt.legend(('J', 'H', 'K$_s$','[3.6]', '[4.5]', '[5.8]', '[8.0]'),
            scatterpoints = 1,
            loc = 'upper right',
            title = 'Filter/Channel',
            fontsize = 13,
            frameon = False)
# Plot time of burst and label it
plt.axvline(x = 55337.8, color = 'k', ls = 'dashed')
#plt.text(55500, 0.45, "2010 May outburst", rotation=90, fontsize=13)
# Plot error bars
plt.errorbar(WIRCjdate, kFlux, kErr, color = 'lightskyblue')
plt.errorbar(WIRCjdate, hFlux, hErr, color = 'royalblue')
plt.errorbar(WIRCjdate, jFlux, jErr, color = 'navy')
plt.errorbar(dates2[0], j2Flux, 0.0593, color = 'navy')
plt.errorbar(jdates, circFlux36, yerr = fluxError36, linestyle = 'None',
             color = 'black')
plt.errorbar(jdates, circFlux45, yerr = fluxError45, linestyle = 'None',
             color = 'grey')
plt.errorbar(jdates[0], circFlux58, yerr = fluxError58, linestyle = 'None',
             color = 'darkgrey')
plt.errorbar(jdates[0], circFlux80, yerr = fluxError80, linestyle = 'None',
             color = 'black')
plt.xlabel('Time (MJD)', fontsize = 14)
plt.ylabel('Flux density (mJy)', fontsize = 14)
ax.arrow(dates2[0], kFluxes[0], 0.0, -0.08, head_width = 150,
         head_length = 0.02, fc ='lightskyblue', ec ='lightskyblue')
ax.arrow(dates2[0], hFluxes[0], 0.0, -0.08, head_width = 150,
         head_length = 0.02, fc = 'royalblue', ec ='royalblue')
ax.arrow(WIRCjdate, 0.7, 0.0, -0.15, head_width = 300, head_length = 0.03,
         fc = 'k', ec = 'k', linestyle = '-')
#plt.text(57250, 0.75, 'TripleSpec spectroscopy', rotation = 'vertical', 
#         fontsize = 12, color = 'red')
x1, x2, y1, y2 = 55100, 56100, 0.01, 0.27

axins = zoomed_inset_axes(ax,1.8,loc=9)
axins.set_xlim(x1,x2)
axins.set_ylim(y1,y2)
plt.scatter(jdates[1:9], circFlux36[1:9], color = 'black', marker = 'o',
            s = 15)
plt.errorbar(jdates[1:9], circFlux36[1:9], yerr = fluxError36[1:9],
             linestyle = 'None', color = 'black')
plt.scatter(jdates[1:9], circFlux45[1:9], color = 'grey', marker = 'v', s = 15)
plt.errorbar(jdates[1:9], circFlux45[1:9], yerr = fluxError45[1:9],
             linestyle = 'None', color = 'grey')
plt.axvline(x = 55337.8, color = 'k', ls = 'dashed')
plt.xticks(np.arange(x1, x2, 400))
mark_inset(ax, axins, loc1 = 3, loc2 = 4, fc = "none", ec = "0.6")

fig.savefig("170215_IC10_X2_smoothed_lc.pdf")