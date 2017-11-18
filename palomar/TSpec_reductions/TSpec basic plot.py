# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 17:14:20 2016

@author: stephaniekwan

Plot prominent emission lines of IC 10 X-2 from TripleSpec chop-subtracted 
data. He I and Pa-Gamma lines fit on one plot, Pa-Beta line goes in a separate
plot (comment/uncomment blocks to plot each set).
"""

import numpy as np
import matplotlib.pyplot as plt

plt.clf()
plt.close()
table = np.genfromtxt('IC10X2_JHK_modified.rtf', delimiter = '    ', 
                  comments = '\p', skip_header = 2, skip_footer = 4)
wl    = table[:, 0] - 0.0005
counts = table[:, 1]
fig = plt.figure()
normFlux = counts / 0.024

 # Two subplots, unpack the axes array immediately
f, (ax1, ax2) = plt.subplots(1,2, sharey = True)
ax1.plot(wl[7100:7500], normFlux[7100:7500], color = 'black')
ax1.invert_xaxis()
ax1.set_xlim([1.075, 1.100])
ax1.set_ylim([0, 2.5])
# Plot and label the original He I line in red
ax1.axvline(x = 1.08303398 - 0.0005, color = 'red', ls = 'dashed')
ax1.text(1.084- 0.0005, 1.5, 'He I (1.0830)', color = 'red', rotation = 90, 
         fontsize = 12)
# Plot and label the peak emission line in green
ax1.axvline(x = 1.08239- 0.0005, color = 'green', ls = 'dashed')
ax1.text(1.08- 0.0005, 1.88, 'He I blueshifted (1.0819)', color = 'green', 
rotation = 90, fontsize = 12)
# Paschen-gamma lines
ax1.axvline(x = 1.093817- 0.0005, color = 'red', ls = 'dashed')
ax1.text(1.095- 0.0005, 1.5, 'Pa$\gamma$ (1.0933)', rotation = 90, 
         fontsize = 12, color = 'red')
ax1.axvline(x = 1.0931- 0.0005, color = 'green', ls = 'dashed')
ax1.text(1.091- 0.0005, 1.5, 'Pa$\gamma$ (1.0926)', 
         rotation = 90, fontsize = 12, color = 'green')
   
# Paschen-beta lines
# Plot the original emission line in red
ax2.plot(wl[5200:5389], normFlux[5200:5389], color = 'black')
ax2.axvline(x = 1.282- 0.0005, color = 'red', ls = 'dashed')
ax2.text(1.283- 0.0005, 1.5, 'Pa $\\beta$ (1.2815)', rotation = 90, 
         fontsize = 12, color = 'red')
# Plot the peak emission line in green
ax2.axvline(x = 1.28103- 0.0005, color = 'green', ls = 'dashed')
ax2.text(1.278- 0.0005, 1.5, 'Pa $\\beta$ (1.2805)', rotation = 90,
         fontsize = 12, color = 'green')
ax2.invert_xaxis()
ax2.set_xlim([1.270, 1.2939]) 
ax2.set_ylim([0, 2.0])

# Set common labels
f.text(0.5, 0.04, 'Wavelength ($\mu$m)', ha = 'center', va = 'center',
       fontsize = 13)
f.text(0.06, 0.5, 'Relative strength to He I line', ha = 'center',
       va = 'center', rotation = 'vertical', fontsize = 13)
 
plt.savefig('170331 TSpec plot.pdf')    