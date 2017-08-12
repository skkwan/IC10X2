# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 13:34:22 2016
Updated 2016 Sep 7 to reformat

@author: stephaniekwan

Calculate the color and magnitude (3.6 - 4.5 microns) of IC 10 X-2,
and compare it to the population of massive stars in the LMC with
IRAC and SAGE measurements.

Source for colormap: http://stackoverflow.com/questions/12965075/matplotlib-sca
tter-plot-colour-as-function-of-third-variable

Source for zoom-in inset plot: http://akuederle.com/matplotlib-zoomed-up-inset
"""
import numpy as np
import matplotlib.pyplot as plt
import Flux_Calculations_160621 as flc
from math import log10
from array import *
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib import cm

# Rename arrays from photometry calculations
f_36     = flc.circFlux36
f_36_sig = flc.fluxError36
f_45     = flc.circFlux45
f_45_sig = flc.fluxError45
f_58     = flc.circFlux58
f_58_sig = flc.fluxError58
f_80     = flc.circFlux80
f_80_sig = flc.fluxError80
jdates   = flc.jdates

# Open the file
f = 'Bonanos 2009 Table 3/table3_only_Ch1_Ch2.dat'
ref_mu = 18.41
dat = np.genfromtxt(f, dtype=(float,float, '|S10'),
                           delimiter='\t', 
                           skip_header=43,
                           missing_values='      ',
                           filling_values='0',
                           usecols=(2, 3, 4))
              
# Make copy of data from which we will remove invalid stars (hence list type)
data     = dat.tolist()
badindex = []          # list of bad indices as ints
# Check if entries have a zero measurement for one of their fluxes
for i in range(len(dat)):
    if 0 in data[i]:
        badindex.append(i)
# Now remove the indices from the data array, in reverse order
for b in badindex[::-1]:
    del data[b]
# Calculate absolute magnitudes based on given distance modulus. 
ref_abs_Mag36, ref_abs_Mag45 = np.zeros((2, len(data)))
for i in range(len(data)):
    ref_abs_Mag36[i] = np.asarray(data[i][0]) - ref_mu
    ref_abs_Mag45[i] = np.asarray(data[i][1]) - ref_mu
ref_color = ref_abs_Mag36 - ref_abs_Mag45

# Now create lists of the different spectral classes.
nclasses = 10   # There are 9 in the paper but make a 10th one for rejects
Ostar, earlyB, lateB, AFG, RSG, WR, sgB_e, LBV, BeXray, reject =  \
[[] for _ in range(nclasses)]
  
# Using the spectral type, sort the objects  
for i in range(len(data)):
    if 'O' in data[i][2] and 'W' not in data[i][2]:
        Ostar.append(data[i])
    elif any(s in data[i][2] for s in ('B0','B1','B2','B3','B4','earl')) and \
        not any(s in data[i][2] for s in ('Be', 'Fe', '[e]', 'K', 'XRB')):
            earlyB.append(data[i])
    elif any(s in data[i][2] for s in ('B5','B6','B7','B8','B9','B extr')) and\
    not any(s in data[i][2] for s in ('[e]')):
        lateB.append(data[i])
    elif any(s in data[i][2] for s in ('A','F','G')) and 'Be' not in \
    data[i][2] and 'LBV' not in data[i][2]:
        AFG.append(data[i])
    elif 'W' in data[i][2]:
        WR.append(data[i])
    elif 'LBV' in data[i][2]:
        LBV.append(data[i])
    elif 'K' in data[i][2] or 'M' in data[i][2]:
        RSG.append(data[i])
    elif '[e]' in data[i][2]:  # or 'Be'?
        sgB_e.append(data[i])
    elif 'XRB' in data[i][2]:  # B extr?
        BeXray.append(data[i])
    elif any(s in data[i][2] for s in ('0','1','2','3','4')):
        earlyB.append(data[i])
    else:
        lateB.append(data[i])

O_col, O_Mag           = np.zeros((2, 1, len(Ostar)))
earlyB_col, earlyB_Mag = np.zeros((2, 1, len(earlyB)))
lateB_col, lateB_Mag   = np.zeros((2, 1, len(lateB)))
AFG_col, AFG_Mag       = np.zeros((2, 1, len(AFG)))
RSG_col, RSG_Mag       = np.zeros((2, 1, len(RSG)))
WR_col, WR_Mag         = np.zeros((2, 1, len(WR)))
sgB_e_col, sgB_e_Mag   = np.zeros((2, 1, len(sgB_e)))
LBV_col, LBV_Mag       = np.zeros((2, 1, len(LBV)))
BeXray_col, BeXray_Mag = np.zeros((2, 1, len(BeXray)))

col_array = [O_col, earlyB_col, lateB_col, AFG_col, RSG_col, WR_col, sgB_e_col,
             LBV_col, BeXray_col]
abs_mag_35array = [O_Mag, earlyB_Mag, lateB_Mag, AFG_Mag, RSG_Mag, WR_Mag,
                   sgB_e_Mag, LBV_Mag, BeXray_Mag]
cat_list = (Ostar, earlyB, lateB, AFG, RSG, WR, sgB_e, LBV, BeXray)
# Calculate the x- and y-values for each star category
for cat in (Ostar, earlyB, lateB, AFG, RSG, WR, LBV, sgB_e, BeXray):
    ind = cat_list.index(cat)
    for i in range(len(cat)):
        col_array[ind][0][i] = np.asarray(cat[i][0]) - np.asarray(cat[i][1])
        abs_mag_35array[ind][0][i] = np.asarray(cat[i][0]) - ref_mu
# Scatter plot the categories using unique shapes and colors   
markershapes = '^', 's', '*', '*', 'v', 'o', 'D', 'o', 'o'
colors = 'black', '#00e5ee','#FF4FA7', '#0000FF', '#FF0000', '#ffffff', \
         '#120A8F', 'black', 'yellow'
edgecolorlist = 'White', 'None', 'None', 'None', 'None', 'black', 'None',\
                'None', 'Black'
catnames = 'O stars', 'Early B', 'Late B', 'AFG', 'RSG', 'WR','sgB[e]', 'LBV',\
           'Be-Xray'
linewidthlist = 0.3, 0, 0, 0, 0, 0.5, 0, 0, 1
sizelist = 20, 10, 20, 30, 12, 12, 12, 30, 50
# Plot results
fig, ax = plt.subplots()
plt.hold(True)
for cat in (earlyB, Ostar, WR, AFG, RSG, sgB_e, lateB, LBV, BeXray):
    ind = cat_list.index(cat)
    ax.scatter(col_array[ind], -abs_mag_35array[ind], c = colors[ind],
               marker = markershapes[ind], linewidths = linewidthlist[ind],
               s = sizelist[ind], edgecolors = edgecolorlist[ind],
               label = catnames[ind])
plt.legend(loc='lower right',fontsize = 10, scatterpoints = 1)
plt.xlabel('$[3.6] - [4.5]$', fontsize = 12)
plt.ylabel('- 3.6$\mu m$ Absolute magnitude', fontsize = 12)
ax.set_xlim(-0.5, 2)
ax.set_ylim(0, 16)

# Handle the IC 10 X-2 fluxes
abs_Mag36, abs_Mag45, app_Mag36, app_Mag45, color = np.zeros((5,len(jdates)))
mu    = 24.1    # Distance modulus of IC 10-X2
f0_36 = 280.9   # zero-point flux for the 3.6 channel
f0_45 = 179.7   # zero-point flux for the 4.5 channel
for j in range(len(jdates)):
    app_Mag36[j] = -2.5 * log10(f_36[j] * 10**-3 / f0_36)
    app_Mag45[j] = -2.5 * log10(f_45[j] * 10**-3 / f0_45)
    abs_Mag36[j] = app_Mag36[j] - mu
    abs_Mag45[j] = app_Mag45[j] - mu
    color[j]     = abs_Mag36[j] - abs_Mag45[j]

# Plot IC 10 X-2's journey across the CMD with a colormap
ic10x2 = plt.scatter(color, -abs_Mag36, marker = '+', c = jdates, 
                     cmap = plt.cm.cool, linewidth = 2, s = 40)
# Add a colorbar
cb = plt.colorbar(ic10x2, shrink = 0.6,
                  ticks = [53500, 54375, 55250, 56125, 57000],
                  orientation = 'vertical')
cb.set_label('Epoch (MJD) of IC 10 X-2', fontsize = 12)
cb.ax.invert_yaxis()

# Make a zoomed inset: zoom-factor: 2.5, location: upper-left
axins = zoomed_inset_axes(ax, 2.5, loc = 1) 
for category in (earlyB, Ostar,WR, AFG, RSG, sgB_e, lateB, LBV, BeXray):
    ind = cat_list.index(category)
    axins.scatter(col_array[ind], -abs_mag_35array[ind], c=colors[ind],
                  marker = markershapes[ind], linewidth = linewidthlist[ind],
                  s = sizelist[ind], edgecolors = edgecolorlist[ind])
axins.scatter(color, -abs_Mag36, marker = '+', c = jdates, cmap = plt.cm.cool,
              s = 40, linewidth = 2)
x1, x2, y1, y2 = 0.0, 0.4, 8, 10
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
plt.yticks(visible = False)
plt.xticks(visible = False)
# position the inset
mark_inset(ax, axins, loc1 = 2, loc2 = 4, fc = "none", ec = "0.5")

plt.savefig('160907 CMD from Bonanos with IC 10 X2 variability.pdf')