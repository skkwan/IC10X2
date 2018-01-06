# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:45:32 2017


subscript 1 = R-band (has more entries)
subscript 2 = g-band
@author: stephaniekwan
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator
from astropy.table import Table

# Use LaTeX font
plt.rc({'weight' : 'normal',
        'size'   : 15})
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('text', usetex = True)

table1 = np.genfromtxt('../optical/PTF_fulltable_phot1.txt', comments = '|',
                        skip_header = 43, skip_footer = 2)

table2 = np.genfromtxt('../optical/PTF_fulltable_phot2.txt', comments = '|',
                       skip_header = 43, skip_footer = 2)                  
                  
flux1, flux2 = table1[:, 7], table2[:, 7] 

sigflux1, sigflux2 = table1[:, 8], table2[:, 8]
mjd1, mjd2 = table1[:, 15], table2[:, 15]
snr1, snr2 = table1[:, 9], table2[:, 9]
zp1, zp2 = table1[:, 5], table2[:, 5]

snt = 3
snu = 5
flux_ref1, sigflux_ref1 = 2771.08, 205.304 #from bottom of table
flux_ref2, sigflux_ref2 = 587.622, 46.0016

mag1, mag1sig, mag1date = np.array([]), np.array([]), np.array([])
upperlim1, upperlim1date = np.array([]), np.array([])
for i in range(len(flux1)):
    # Section 9: define new sigflux DC
    if (sigflux1[i] > sigflux_ref1):
        sigflux_DC = np.sqrt(sigflux1[i] ** 2 - sigflux_ref1 ** 2)
    else:
        sigflux_DC = np.sqrt(sigflux1[i] ** 2 + sigflux_ref1 ** 2)
    # Section 9: redefine SNR
    newSnr = (flux1[i] + flux_ref1) / sigflux_DC
    if (newSnr > snt):   # we have a confident detection
        mag1 = np.append(mag1, 
                         zp1[i] - 2.5 * np.log10(flux1[i] + flux_ref1))
        mag1sig = np.append(mag1sig, 
                            1.0875 * sigflux1[i] / (flux1[i] + flux_ref1))
        mag1date = np.append(mag1date, mjd1[i])
    else:
        # compute upper flux limit and plot as arrow or triangle
        upperlim1 = np.append(upperlim1, 
                              zp1[i] - 2.5 * np.log10(snu * sigflux1[i]))
        upperlim1date = np.append(upperlim1date, mjd1[i])

toosmall = []
for i in range(0, len(mag1)):
    if mag1[i] < 10.0:
        toosmall.append(i)

for i in toosmall[::-1]:
    mag1 = np.delete(mag1, i)
    mag1date = np.delete(mag1date, i)
    mag1sig = np.delete(mag1sig, i)

gs = gridspec.GridSpec(2, 2, width_ratios = [5, 1]) 
ax1 = plt.subplot(gs[0])
ax1 = plt.gca() 
ax1.text(-0.07, 1.1, '\\textbf{(a)}', transform=ax1.transAxes,
      fontsize = 15, fontweight = 'bold', va = 'top', ha = 'right')
      
plt.scatter(mag1date, mag1, marker = 'o', s = 2, color = 'black', zorder = 3)
#plt.scatter(upperlim1date, upperlim1, color = 'grey', marker = 'v',
#            facecolors = 'grey', s = 15, zorder = 4)
for i in range(0, len(upperlim1)):
    ax1.arrow(upperlim1date[i], upperlim1[i],
              0.0, 0.3, head_width = 20, head_length = 0.15,
              fc = 'grey', ec = 'grey', linestyle = '-')        
            
plt.errorbar(mag1date, mag1, yerr = mag1sig, linestyle = 'None',
             color = 'grey', linewidth = 1, zorder = 2)
             
plt.axhline(y = np.median(mag1), color = 'k', ls = ':')

xlowerlim1, xupperlim1, ylowerlim1, yupperlim1  = 56400, 57800, 17.0, 20.0

ax1.set_xlim([xlowerlim1, xupperlim1])  
ax1.set_ylim([ylowerlim1, yupperlim1])    
ax1.invert_yaxis()    
plt.xlabel('Date (MJD)', fontsize = 14)
plt.ylabel('R Magnitude', fontsize = 14)
minorLocator = AutoMinorLocator()
minorLocator2= AutoMinorLocator()
ax1.xaxis.set_minor_locator(minorLocator)
ax1.yaxis.set_minor_locator(minorLocator2)

# Shaded area to denote uncertainty of median (average of mag1sig)
ax1.add_patch(                     
    patches.Rectangle(
        (xlowerlim1, np.median(mag1) - 5 * np.average(mag1sig)),  # (x, y)
        xupperlim1 - xlowerlim1,   # width 
        10 * np.average(mag1sig),   # height
        0.0,                       # angle
        facecolor = 'lightgrey',
        edgecolor = 'none',
        zorder = 1
        ))


# Add inset
i1, i2 = 4, 12
x1, x2 = mag1date[i1], mag1date[i2]

axins = plt.subplot(gs[1])
axins.set_xlim(x1 + 7, x2 + 10)
axins.set_xticks(np.arange(56480, 56514, 10))
axins.set_ylim(18.35, 18.85)
plt.scatter(mag1date[i1:i2], mag1[i1:i2], color = 'black', marker = 'o',
            s = 4, zorder = 3)
plt.errorbar(mag1date[i1:i2], mag1[i1:i2], yerr = mag1sig[i1:i2],
             linestyle = 'None', color = 'black', zorder = 2)
plt.axhline(y = np.median(mag1), color = 'k', ls = ':')
axins.invert_yaxis()
minorLocator3 = AutoMinorLocator()
minorLocator4 = AutoMinorLocator()
axins.xaxis.set_minor_locator(minorLocator3)
axins.yaxis.set_minor_locator(minorLocator4)

# Inset: Shaded area to denote uncertainty of median (average of mag1sig)
axins.add_patch(                     
    patches.Rectangle(
        (x1 - 10, np.median(mag1) - 5 * np.average(mag1sig)),  # (x, y)
        xupperlim1 - xlowerlim1,   # width 
        10 * np.average(mag1sig),   # height
        0.0,                       # angle
        facecolor = 'lightgrey',
        edgecolor = 'none',
        zorder = 1
        ))

###########################
## g-band data
###########################
mag2, mag2sig, mag2date = np.array([]), np.array([]), np.array([])
upperlim2, upperlim2date = np.array([]), np.array([])
for i in range(len(flux2)):
    # Section 9: define new sigflux DC
    if (sigflux2[i] > sigflux_ref2):
        sigflux_DC = np.sqrt(sigflux2[i] ** 2 - sigflux_ref2 ** 2)
    else:
        sigflux_DC = np.sqrt(sigflux2[i] ** 2 + sigflux_ref2 ** 2)
    # Section 9: redefine SNR
    newSnr = (flux2[i] + flux_ref2) / sigflux_DC
    if (newSnr > snt):   # we have a confident detection
        mag2 = np.append(mag2, 
                         zp2[i] - 2.5 * np.log10(flux2[i] + flux_ref2))
        mag2sig = np.append(mag2sig, 
                            1.0875 * sigflux2[i] / (flux2[i] + flux_ref2))
        mag2date = np.append(mag2date, mjd2[i])
    else:
        # compute upper flux limit and plot as arrow or triangle
        upperlim2 = np.append(upperlim2, 
                              zp2[i] - 2.5 * np.log10(snu * sigflux2[i]))
        upperlim2date = np.append(upperlim2date, mjd1[i])
toosmall2 = []
for i in range(0, len(mag2)):
    if mag2[i] < 10.0:
        toosmall2.append(i)
        
for i in toosmall2[::-1]:
    mag2 = np.delete(mag2, i)
    mag2date = np.delete(mag2date, i)
    mag2sig = np.delete(mag2sig, i)

ax2 = plt.subplot(gs[2])

ax2.text(-0.07, 1.15, '\\textbf{(b)}', transform = ax2.transAxes,
      fontsize = 15, fontweight = 'bold', va = 'top', ha = 'right')
xlowerlim2, xupperlim2, ylowerlim2, yupperlim2 = 55000, 58000, 18.5, 21.0
ax2.set_xlim([xlowerlim2, xupperlim2])
ax2.set_ylim([ylowerlim2, yupperlim2])
plt.axhline(y = np.median(mag2), color = 'k', ls = ':')
plt.scatter(mag2date, mag2, marker = 'o', s = 5, color = 'black', zorder = 3)
plt.errorbar(mag2date, mag2, yerr = mag2sig, linestyle = 'None', 
             color = 'grey', zorder = 2)
ax2.invert_yaxis()
# g-band: Shaded area to denote uncertainty of median 
ax2.add_patch(                     
    patches.Rectangle(
        (xlowerlim2, np.median(mag2) - 5 * np.average(mag2sig)),  # (x, y)
        xupperlim2 - xlowerlim2,   # width 
        10 * np.average(mag2sig),   # height
        0.0,                       # angle
        facecolor = 'lightgrey',
        edgecolor = 'none',
        zorder = 1
        ))             

#plt.scatter(upperlim2date, upperlim2, color = 'grey', marker = 'v',
#            facecolors = 'none',
#            s = 15)
# Plot the upper limits of non-detections
for i in range(0, len(upperlim2)):
    ax2.arrow(upperlim2date[i], upperlim2[i],
              0.0, 0.3, head_width = 30, head_length = 0.15,
              fc = 'grey', ec = 'grey', linestyle = '-') 
              
# Plot time of burst and label it
plt.axvline(x = 55337.8, color = 'k', ls = 'dashed')
plt.text(55380, 19.0, 'X-ray outburst', rotation = 'vertical', 
         fontsize = 12, color = 'black')          
plt.xlabel('Date (MJD)', fontsize = 14)
plt.ylabel('g Magnitude', fontsize = 14)

minorLocator5 = AutoMinorLocator()
minorLocator6 = AutoMinorLocator()
ax2.xaxis.set_minor_locator(minorLocator5)
ax2.yaxis.set_minor_locator(minorLocator6)

#plt.show()
#plt.tight_layout()
#plt.savefig("170705_Optical_lc.pdf")


######################
### Write out data to file
######################
Rdata_sorted = np.dstack((mag1date, mag1, mag1sig))
gdata_sorted = np.dstack((mag2date, mag2, mag2sig))

np.savetxt('170501_R_mags.txt', np.column_stack((mag1date, mag1, mag1sig)),
           header = 'Date(MJD) Magnitude uncertainty', fmt = '%10.5f',
            delimiter = '\t') 
np.savetxt('170501_g_mags.txt', np.column_stack((mag2date, mag2, mag2sig)),
          header = 'Date(MJD) Magnitude uncertainty', fmt = '%10.5f',
            delimiter = '\t') 
# Optional arguments:
# header = 'PTF IC 10 X-2 R band data: median = %f +/- %f magnitudes. Date (MJD), magnitude, and uncertainty ' 
#           % (np.median(mag1), np.average(mag1sig)  
# header = 'PTF IC 10 X-2 g band data: median = %f +/- %f magnitudes. Date (MJD), magnitude, and uncertainty ' 
#           % (np.median(mag2), np.average(mag2sig))           
#######################
## Write out data to LaTeX file
#######################
t1 = Table.read('170501_R_mags.txt', format = 'ascii')
t2 = Table.read('170501_g_mags.txt', format = 'ascii')
t1.write('R band LaTeX.txt', format = 'latex', overwrite = True)           
t2.write('g band LaTeX.txt', format = 'latex', overwrite = True)   
               
