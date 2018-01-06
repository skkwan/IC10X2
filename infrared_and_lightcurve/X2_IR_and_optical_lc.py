# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 14:43:57 2017
Updated Fri Aug 11 2017 to use github repository directories.

@author: stephaniekwan

Code for producing IR and optical combined light curves for IC 10 X-2.
"""
import sys
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
import matplotlib.patches as patches

# Get optical fluxes
sys.path.insert(0, '../optical')
import Optical_lc_041217 as opt
# Get infrared fluxes
sys.path.insert(0, '../infrared_and_lightcurve')
import fluxes_with_xdates_170314 as ir

fig, ax = plt.subplots()

gs = gridspec.GridSpec(1, 2, width_ratios = [5, 3]) 
ax = plt.subplot(gs[0])

# Plot IR MAGNITUDE light curves
plt.scatter(ir.dates2, ir.jMags, facecolors = 'navy', marker = '<', s = 25,
            edgecolors = 'navy', zorder = 3, label = 'J')
plt.scatter(ir.dates2, ir.hMags, facecolors = 'none', marker = 's', s = 25,
             edgecolors = 'royalblue', zorder = 3, label = 'H')
plt.scatter(ir.dates2, ir.kMags, facecolors = 'none', marker = '>', s = 25,
            edgecolors = 'lightskyblue', zorder = 3, label = 'K')
# Add error arrows to the H and K 2MASS images
ax.arrow(ir.dates2[0], ir.kMags[0], 0.0, 0.3, head_width = 100,
         head_length = 0.15, fc ='lightskyblue', ec ='lightskyblue', zorder = 2)
ax.arrow(ir.dates2[0], ir.hMags[0], 0.0, 0.3, head_width = 100,
         head_length = 0.15, fc = 'royalblue', ec ='royalblue', zorder = 1)
plt.errorbar(ir.dates2[0], ir.jMags[0], yerr = ir.jmag1sig, linestyle = 'None',
            color = 'navy', zorder = 3) 
# Plot errorbars for WIRC datapoints
plt.errorbar(ir.dates2[1], ir.jMags[1], yerr = ir.jmag2sig, linestyle = 'None',
             color = 'navy', zorder = 3) 
plt.errorbar(ir.dates2[1], ir.hMags[1], yerr = ir.hmag2sig, linestyle = 'None',
             color = 'royalblue', zorder = 3)   
plt.errorbar(ir.dates2[1], ir.kMags[1], yerr = ir.kmag2sig, linestyle = 'None',
             color = 'lightskyblue', zorder = 3)  
# Plot errorbars for Spitzer datapoints
plt.errorbar(ir.jdates, ir.m36, yerr = ir.m36sig, linestyle = 'None',
             color = 'k', zorder = 2)
plt.errorbar(ir.jdates, ir.m45, yerr = ir.m45sig, linestyle = 'None',
             color = 'grey', zorder = 2)
plt.errorbar(ir.jdates[0], ir.m58, yerr = ir.m58sig, linestyle = 'None',
             color = 'grey', zorder = 1)
plt.errorbar(ir.jdates[0], ir.m80, yerr = ir.m80sig, linestyle = 'None',
             color = 'grey', zorder = 1)             
# Plot the Spitzer datapoints
plt.scatter(ir.jdates, ir.m36, color = 'black', marker = 'o', s = 15, 
            zorder = 3, label = '3.6')
plt.scatter(ir.jdates, ir.m45, color = 'grey', marker = 'v', s = 15,
            zorder = 3, label = '4.5')
plt.scatter(ir.jdates[0], ir.m58, facecolors = 'none', edgecolors =
             'grey', marker = 'D', s = 22, zorder = 3, label = '5.8')
plt.scatter(ir.jdates[0], ir.m80, facecolors ='none', edgecolors = 'grey',
            marker = 'o', s = 25, zorder = 3, label = '8.0')
            
# Plot non-detections
for i in range(len(ir.xdates)):     
       ax.arrow(ir.xdates[i], 17.1, 0.0, -0.2, head_width = 100, head_length = 0.1,
         fc = 'darkslategrey', ec = 'darkslategrey', linestyle = '-')     
for j in range(len(ir.xnondates)):
       ax.arrow(ir.xnondates[j], 17.1, 0.0, -0.2, head_width = 30, head_length = 0.1,
        fc = 'lightgrey', ec = 'lightgrey', linestyle = '-')            
# Plot quiescent magnitudes
quiesmag36, quiesmag45 = np.mean(ir.m36), np.mean(ir.m45)
plt.axhline(y = quiesmag36, color = 'black', ls = 'dashed')
plt.axhline(y = quiesmag45, color = 'grey', ls = 'dashed')
#plt.text(52650, 15.5, 'Quiescent levels', fontsize = 10)

entirelowerlim = 51500
entireupperlim = 58100   
# Shade rectangles around it
ax.add_patch(                     
    patches.Rectangle(
        (entirelowerlim, quiesmag36 - np.average(ir.m36sig)),  # (x, y)
        entireupperlim - entirelowerlim,   # width 
        2 * np.average(ir.m36sig),   # height
        0.0,                       # angle
        facecolor = 'gainsboro',
        edgecolor = 'none',
        zorder = 1,
        alpha = 0.95
        ))   
ax.add_patch(                     
    patches.Rectangle(
        (entirelowerlim, quiesmag45 - np.average(ir.m45sig)),  # (x, y)
        entireupperlim - entirelowerlim,   # width 
        2 * np.average(ir.m45sig),   # height
        0.0,                       # angle
        facecolor = 'lightgrey',
        edgecolor = 'none',
        zorder = 2,
        alpha = 0.4
        ))   

# Plot time of burst and label it
plt.axvline(x = 55337.8, color = 'k', ls = 'dashed')
#plt.text(55500, 15.6, "2010 May outburst", rotation=90, fontsize=11)
 

#################
## Optical points: R band
#################
ax = plt.gca() 
plt.scatter(opt.mag1date, opt.mag1, marker = 'o', s = 2, color = 'black',
            zorder = 4, label = 'R')
for i in range(0, len(opt.upperlim1)):          # Upper limits
    ax.arrow(opt.upperlim1date[i], opt.upperlim1[i],
              0.0, 0.3, head_width = 20, head_length = 0.15,
              fc = 'grey', ec = 'grey', linestyle = '-')        
            
plt.errorbar(opt.mag1date, opt.mag1, yerr = opt.mag1sig, linestyle = 'None',
             color = 'grey', linewidth = 1, zorder = 3)
             
plt.axhline(y = np.median(opt.mag1), color = 'k', ls = ':')

# Shaded area to denote uncertainty of median (average of mag1sig)
ax.add_patch(                     
    patches.Rectangle(
        (entirelowerlim, np.median(opt.mag1) - 5 * np.average(opt.mag1sig)),  # (x, y)
        entireupperlim - entirelowerlim,   # width 
        10 * np.average(opt.mag1sig),   # height
        facecolor = 'lightgrey',
        edgecolor = 'none',
        zorder = 1
        ))
ax.add_patch(                     
    patches.Rectangle(
        (entirelowerlim, np.median(opt.mag1) - np.average(opt.mag1sig)),  # (x, y)
        entireupperlim - entirelowerlim,   # width 
        2 * np.average(opt.mag1sig),   # height
        0.0,                       # angle
        facecolor = 'lightblue', edgecolor = 'none', zorder = 2))      

# Labels
plt.text(55380, 18.4, '5$\sigma$', fontsize = 10)
plt.text(55380, 18.6, '1$\sigma$', fontsize = 10)
               
#################
## Optical points: g band
#################
plt.scatter(opt.mag2date, opt.mag2, marker = 'o', s = 5, color = 'red',
            zorder = 3, label = 'g')
plt.errorbar(opt.mag2date, opt.mag2, yerr = opt.mag2sig, linestyle = 'None', 
             color = 'grey', zorder = 2)
ax.add_patch(                     
    patches.Rectangle(
        (entirelowerlim, np.median(opt.mag2) - 5 * np.average(opt.mag2sig)),  # (x, y)
        entireupperlim - entirelowerlim,   # width 
        10 * np.average(opt.mag2sig),   # height
        facecolor = 'lightgrey',
        edgecolor = 'none',
        zorder = 1
        ))             
             
ax.add_patch(                     
    patches.Rectangle(
        (entirelowerlim, np.median(opt.mag2) - np.average(opt.mag2sig)),  # (x, y)
        entireupperlim - entirelowerlim,   # width 
        2 * np.average(opt.mag2sig),   # height
        0.0,                       # angle
        facecolor = 'lightblue',  edgecolor = 'none', zorder = 1))             

plt.axhline(y = np.median(opt.mag2), color = 'k', ls = ':')

##############
## Finishing touches
##############
# Set magnitude plot limits       
plt.xlim([entirelowerlim, entireupperlim])
plt.ylim([13.5, 21.0])


plt.xlabel('Date (MJD)', fontsize = 14)
plt.ylabel('Magnitude', fontsize = 14)

# Reverse y axis
plt.gca().invert_yaxis()

#################
## Add inset
#################
axins = plt.subplot(gs[1])
# sub region of the original image
x1, x2, y1, y2 = 56450, 57500, 17.25, 19.7
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
plt.scatter(opt.mag1date, opt.mag1, marker = 'o', s = 2, color = 'black',
            zorder = 4)
plt.gca().invert_yaxis()

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
# Broken in Python 3 now :(
#mark_inset(ax, axins, loc1=2, loc2=3, fc="none", ec="0.5")

# Create a rectangle patch in the original plot
rect = patches.Rectangle((x1, y1), x2 - x1, y2 - y1,
                         linewidth = 1, edgecolor = 'grey', facecolor = 'none')
ax.add_patch(rect)

plt.errorbar(opt.mag1date, opt.mag1, yerr = opt.mag1sig, linestyle = 'None',
             color = 'grey', linewidth = 1, zorder = 3)
             
plt.axhline(y = np.median(opt.mag1), color = 'k', ls = ':')

# Shaded area to denote uncertainty of median (average of mag1sig)
axins.add_patch(                     
    patches.Rectangle(
        (entirelowerlim, np.median(opt.mag1) - 5 * np.average(opt.mag1sig)),  # (x, y)
        entireupperlim - entirelowerlim,   # width 
        10 * np.average(opt.mag1sig),   # height
        facecolor = 'lightgrey',
        edgecolor = 'none',
        zorder = 1
        ))
axins.add_patch(                     
    patches.Rectangle(
        (entirelowerlim, np.median(opt.mag1) - np.average(opt.mag1sig)),  # (x, y)
        entireupperlim - entirelowerlim,   # width 
        2 * np.average(opt.mag1sig),   # height
        0.0,                       # angle
        facecolor = 'lightblue', edgecolor = 'none', zorder = 2))   

## Plot legend
handles, labels = ax.get_legend_handles_labels()
labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
lgd = ax.legend(handles, labels, loc = 'lower left', title = 'Filter/Channel',
                scatterpoints = 1,
                fontsize = 10, frameon = True)

plt.xlabel('Date (MJD)', fontsize = 14)
plt.ylabel('Magnitude', fontsize = 14)


plt.tight_layout() 
plt.show()           

fig.savefig("X2_IR_and_optical_lc.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')


