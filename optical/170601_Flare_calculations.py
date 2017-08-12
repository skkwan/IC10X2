# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 13:21:57 2017

Generates spaced-out gridded plot for counting flares. Plot parameters need to be adjusted by hand
and output file needs to be manually renamed each time (sorry)

@author: stephaniekwan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
#from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib import gridspec
#----------------------------------
#Get fluxes
#----------------------------------
import Optical_lc_041217 as f

rMag = f.mag1
rDate = f.mag1date
rMagSig = f.mag1sig
upperlimRDate , upperlimR = f.upperlim1date, f.upperlim1

gMag = f.mag2
gDate = f.mag2date
gMagSig = f.mag2sig
upperlimGDate, upperlimG = f.upperlim2date, f.upperlim2

gs = gridspec.GridSpec(1, 1, width_ratios = [5, 1]) 
ax1 = plt.subplot(gs[0])
ax1 = plt.gca() 
      
plt.scatter(rDate, rMag, marker = 'o', s = 2, color = 'black', zorder = 3)
#plt.scatter(upperlim1date, upperlim1, color = 'grey', marker = 'v',
#            facecolors = 'grey', s = 15, zorder = 4)
for i in range(0, len(upperlimRDate)):
    ax1.arrow(upperlimRDate[i], upperlimR[i],
              0.0, 0.3, head_width = 20, head_length = 0.15,
              fc = 'grey', ec = 'grey', linestyle = '-')        
            
plt.errorbar(rDate, rMag, yerr = rMagSig, linestyle = 'None',
             color = 'grey', linewidth = 1, zorder = 2)
             
plt.axhline(y = np.median(rMag), color = 'k', ls = ':')

xlowerlim1, xupperlim1, ylowerlim1, yupperlim1  = 57000, 57300, 17.0, 19.5

ax1.set_xlim([xlowerlim1, xupperlim1])  
ax1.set_ylim([ylowerlim1, yupperlim1]) 
ax1.set_xticks(np.arange(xlowerlim1, xupperlim1, 10))  
plt.xticks(rotation = 45, fontsize = 6) 
ax1.tick_params(axis = 'x', direction = 'out')
ax1.grid()
ax1.invert_yaxis()    
plt.xlabel('Date (MJD)', fontsize = 14)
plt.ylabel('R Magnitude', fontsize = 14)


# Shaded area to denote uncertainty of median (average of mag1sig)
ax1.add_patch(                     
    patches.Rectangle(
        (xlowerlim1, np.median(rMag) - 5 * np.average(rMagSig)),  # (x, y)
        xupperlim1 - xlowerlim1,   # width 
        10 * np.average(rMagSig),   # height
        facecolor = 'lightgrey',
        edgecolor = 'none',
        zorder = 1
        ))
 
# Shaded area to denote uncertainty of median (average of mag1sig)
ax1.add_patch(                     
    patches.Rectangle(
        (xlowerlim1, np.median(rMag) - 1 * np.average(rMagSig)),  # (x, y)
        xupperlim1 - xlowerlim1,   # width 
        2 * np.average(rMagSig),   # height
        facecolor = 'lightblue',
        edgecolor = 'none',
        zorder = 1
        ))       
        
plt.show()
plt.savefig('Counting flares/Rfig3.png')
