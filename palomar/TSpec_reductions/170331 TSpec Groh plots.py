# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 15:48:38 2017

@author: stephaniekwan
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 17:14:20 2016

@author: stephaniekwan

Plot prominent emission lines of IC 10 X-2 from TripleSpec chop-subtracted 
data. He I and Pa-Gamma lines fit on one plot, Pa-Beta line goes in a separate
plot (comment/uncomment blocks to plot each set).

DEPENDENCIES: Subfolder titled "groh_data", with .txt files of additional 
reference spectra to plot.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


plt.rc('text', usetex=True)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})

plt.rc({'weight' : 'normal',
        'size'   : 15})


plt.clf()
plt.close()
table = np.genfromtxt('IC10X2_JHK_modified.rtf', delimiter = '    ', 
                  comments = '\p', skip_header = 2, skip_footer = 4)
wl    = table[:, 0] 
counts = table[:, 1]
fig = plt.figure()
normFlux = counts / 0.024               # Normalize to height of He I being 1

lbv1_table = np.genfromtxt('groh_data/w243.txt', delimiter = ' ', 
                  comments = '#', skip_header = 1)
lbv1_wl = lbv1_table[:, 0] * 10**(-4)           # Convert angstrom to microns
lbv1_normFlux = lbv1_table[:, 2] / 2.87104      # Normalize    
lbv2_table = np.genfromtxt('groh_data/hd316285.txt', delimiter = ' ',
                          comments = '#', skip_header = 1)
lbv2_wl = lbv2_table[:, 0] * 10**(-4)
lbv2_normFlux = lbv2_table[:, 2] / 22.7456
lbv3_table = np.genfromtxt('groh_data/agcar.txt', delimiter = ' ', 
                  comments = '#', skip_header = 1)
lbv3_wl = lbv3_table[:, 0] * 10**(-4)
lbv3_normFlux = lbv3_table[:, 2] / 16.28              

wr1_table = np.genfromtxt('groh_data/wr6_2.txt', delimiter = ' ', 
                  comments = '#', skip_header = 1)
wr1_wl = wr1_table[:, 0] * 10**(-4)           # Convert angstrom to microns
wr1_normFlux = wr1_table[:, 2]   
wr2_table = np.genfromtxt('groh_data/wr136.txt', delimiter = ' ',
                          comments = '#', skip_header = 1)
wr2_wl = wr2_table[:, 0] * 10**(-4)
wr2_normFlux = wr2_table[:, 2] 
wr3_table = np.genfromtxt('groh_data/wr90.txt', delimiter = ' ',
                          comments = '#', skip_header = 1)
wr3_wl = wr3_table[:, 0] * 10**(-4)
wr3_normFlux = wr3_table[:, 2]                          

ob1_table =  np.genfromtxt('groh_data/hd80077.txt', delimiter = ' ', 
                  comments = '#', skip_header = 1)     
ob1_wl = ob1_table[:, 0] * 10**(-4)           # Convert angstrom to microns
ob1_normFlux = ob1_table[:, 2]    
ob2_table = np.genfromtxt('groh_data/hd169454.txt', delimiter = ' ',
                          comments = '#', skip_header = 1)
ob2_wl = ob2_table[:, 0] * 10**(-4)                          
ob2_normFlux = ob2_table[:, 2]   
ob3_table = np.genfromtxt('groh_data/hd152408.txt', delimiter = ' ',
                          comments = '#', skip_header = 1)  
ob3_wl = ob3_table[:, 0] * 10**(-4)
ob3_normFlux = ob3_table[:, 2]                         

be1_table = np.genfromtxt('groh_data/xper.txt', delimiter = ' ',
                          comments = '#', skip_header = 1)
be1_wl = be1_table[:, 0]* 10**(-4)
be1_normFlux = be1_table[:, 2]    
be2_table = np.genfromtxt('groh_data/delta_sco.txt', delimiter = ' ',
                          comments = '#', skip_header = 1)
be2_wl = be2_table[:, 0] * 10**(-4)
be2_normFlux = be2_table[:, 2]
be3_table = np.genfromtxt('groh_data/chi_oph.txt', delimiter = ' ',
                          comments = '#', skip_header = 1)
be3_wl = be3_table[:, 0] * 10**(-4)
be3_normFlux = be3_table[:, 2]                          

 # Two subplots, unpack the axes array immediately
f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2) 
axlargest = plt.gca()

###########################################
# IC10X2 plots
###########################################
ax1.set_title('IC 10 X-2', fontsize = 14)
ax1.text(0.03, 1.3, '\\textbf{(a)}', transform=ax1.transAxes,
      fontsize = 14, fontweight = 'bold', va = 'top', ha = 'right')
# Plot IC10X2 spectra
ax1.plot(wl[7100:7500], normFlux[7100:7500], color = 'black')
ax1.invert_xaxis()
ax1.set_xlim([1.075, 1.100])
ax1.set_ylim([0, 1.9])

# Plot and label the original He I line in red
ax1.axvline(x = 1.08303398, color = 'red', ls = 'dashed')
# Plot and label the peak He I emission line in green 1.08239
ax1.axvline(x = 1.08236, color = 'green', ls = 'dashed')
ax1.text(1.08, 1.2, 'He I', color = 'green', rotation = 90, fontsize = 12)
# Paschen-gamma lines
ax1.axvline(x = 1.093817, color = 'red', ls = 'dashed')
ax1.axvline(x = 1.0931, color = 'green', ls = 'dashed')
ax1.text(1.091, 1.2, 'Pa $\gamma$', rotation = 90, fontsize = 12, color = 'green') 


# Remove walls and ticks between the IC 10 X-2 plots
ax1.spines['right'].set_visible(False) 
ax2.spines['left'].set_visible(False)
ax1.yaxis.tick_left()
ax1.tick_params(labeltop = 'off') #don't put tick labels at the top
ax2.yaxis.tick_right()
d =  .04        # Parameter for the slashes on the top and bottom frames
kwargs = dict(transform= ax1.transAxes, color='k', clip_on=False)
ax1.plot((1-d,1+d),(-d,+d), **kwargs) # top-left diagonal
ax1.plot((1-d,1+d),(1-d,1+d), **kwargs) # bottom-left diagonal
# Switch to the bottom axes
kwargs.update(transform=ax2.transAxes) # switch to the bottom axes
ax2.plot((-d,d),(-d,+d), **kwargs) # top-right diagonal
ax2.plot((-d,d),(1-d,1+d), **kwargs) # bottom-right diagonal
   
# Paschen-beta lines
# Plot the original emission line in red
ax2.plot(wl[5200:5389], normFlux[5200:5389], color = 'black')
ax2.axvline(x = 1.282, color = 'red', ls = 'dashed')
#ax2.text(1.283- 0.0005, 1.5, 'Pa $\\beta$ (1.2815)', rotation = 90, 
#         fontsize = 12, color = 'red')
# Plot the peak emission line in green
ax2.axvline(x = 1.28103, color = 'green', ls = 'dashed')
ax2.text(1.2786, 1.2, 'Pa $\\beta$', rotation = 90,
         fontsize = 12, color = 'green')
ax2.invert_xaxis()
ax2.set_xlim([1.270, 1.2939]) 
ax2.set_ylim([0, 1.9])

offset = 0.000717

###########################################
# LBV/LBV candidate plot
###########################################
min2, max2   = 570, 570 + 394
ax3.plot(lbv1_wl[min2:max2], lbv1_normFlux[min2:max2]  + 0.4,
         color = 'grey', linewidth = 0.75)
ax3.plot(lbv2_wl[434:828], lbv2_normFlux[434:828] + 0.4, '-',
         color = 'grey', dashes= (2, 0.4), linewidth = 0.5)
ax3.plot(lbv3_wl[441:835], lbv3_normFlux[441:835] + 1, color = 'black', linewidth = 0.5)
ax3.plot(wl[7100:7500] + 0.00064397,
         normFlux[7100:7500], color = 'black') # IC 10 X-2
ax3.set_yticks(np.arange(0.0, 2.0, 0.5))
ax3.set_xlim([1.075, 1.100])
ax3.set_ylim([0.0, 2.0])
ax3.text(0.03, 1.3, '\\textbf{(b)}', transform = ax3.transAxes,
      fontsize = 14, fontweight = 'bold', va = 'top', ha = 'right')
ax3.set_title('LBV/LBV candidate', fontsize = 14)


# Wolf-Rayet star plot
min3, max3 = 552, 946
ax4.plot(wr1_wl[529:924], wr1_normFlux[529:924]+1.2, color = 'black', linewidth = 0.5)
ax4.plot(wr2_wl[513:908], wr2_normFlux[513:908]+0.0, '-', linewidth = 0.5,
         color = 'grey', dashes=(2, 0.4))
ax4.plot(wr3_wl[503:898], wr3_normFlux[503:898],color = 'grey', linewidth = 0.75)
ax4.plot(wl[7100:7500] + offset,
         normFlux[7100:7500] + 1, color = 'black') # IC 10 X-2
ax4.set_ylim([0.0, 7.0])
ax4.set_xlim([1.075, 1.100])
ax4.set_yticks(np.arange(0, 7.0, 2.0))
ax4.set_title("WR star", fontsize = 14)
ax4.text(0.03, 1.3, '\\textbf{(c)}', transform = ax4.transAxes,
      fontsize = 14, fontweight = 'bold', va = 'top', ha = 'right')
      
# OB star plot
ax5.plot(ob1_wl[502:897], ob1_normFlux[502:897] + 0.9, 
         '-', linewidth = 0.5,
         color = 'grey', dashes=(2, 0.4))
#'-', color = 'grey', dashes = (5, 1))
ax5.plot(ob2_wl[512:964], ob2_normFlux[512:964] + 1.8, color = 'black', linewidth = 0.5)
ax5.plot(ob3_wl[570:964], ob3_normFlux[570:964], color = 'grey', linewidth = 0.75 )
ax5.plot(wl[7100:7500] + offset,
         normFlux[7100:7500] + 0.5, color = 'black') # IC 10 X-2
ax5.set_ylim([0.5, 4.0])
ax5.set_title('OB star', fontsize = 14) 
ax5.set_xlim([1.075, 1.100])
ax5.text(0.03, 1.3, '\\textbf{(d)}', transform = ax5.transAxes,
      fontsize = 14, fontweight = 'bold', va = 'top', ha = 'right')
ax5.set_yticks(np.arange(0.5, 4.0, 1.0))      

# Be star plot
ax6.plot(be1_wl[529:924], be1_normFlux[529:924] + 1.5,  
         color = 'black', linewidth = 0.5)
ax6.plot(be2_wl[502:897], be2_normFlux[502:897] + 0.6, '-', linewidth = 0.5,
         color = 'grey', dashes=(2, 0.4))
ax6.plot(be3_wl[502:897], be3_normFlux[502:897], color = 'grey', linewidth = 0.75)
ax6.plot(wl[7100:7500] + offset,
         normFlux[7100:7500] + 0.5, color = 'black', linewidth = 1.0 ) # IC 10 X-2
ax6.set_ylim([0.5, 3.4])
ax6.set_xlim([1.075, 1.100])
ax6.set_title('Be star', fontsize = 14)
ax6.text(0.03, 1.3, '\\textbf{(e)}', transform = ax6.transAxes,
      fontsize = 14, fontweight = 'bold', va = 'top', ha = 'right')
ax6.set_yticks(np.arange(0.5, 3.5, 1.0))

# Set common labels
f.text(0.5, 0.03, 'Wavelength ($\mu$m)', ha = 'center', va = 'center',
       fontsize = 14)
f.text(0.05, 0.5, 'Flux (arbitrary units)', ha = 'center',
       va = 'center', rotation = 'vertical', fontsize = 14)

plt.subplots_adjust(wspace = 0.2, hspace= 0.7)
#plt.tight_layout()
 
plt.savefig('171111 TSpec and Groh plot.pdf')    