# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 13:55:16 2016 as a copy+paste of the previous version of
this file. For making a journal figure.

Calculate the color and magnitude (3.6 - 4.5 microns) of IC 10 X-2,
and compare it to the population of massive stars in the LMC with
IRAC and SAGE measurements.

Source for colormap: http://stackoverflow.com/questions/12965075/matplotlib-sca
tter-plot-colour-as-function-of-third-variable

Source for zoom-in inset plot: http://akuederle.com/matplotlib-zoomed-up-inset
"""
import numpy as np
import matplotlib.pyplot as plt
import fluxes_with_xdates_170314_copy as flc
from math import log10
#from array import *
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex = True)
from matplotlib.ticker import MultipleLocator

# Rename arrays from photometry calculations
f_36, f_37_sig = flc.circFlux36, flc.fluxError36
f_45, f_45_sig = flc.circFlux45, flc.fluxError45
f_58, f_58_sig = flc.circFlux58, flc.fluxError58
f_80, f_80_sig = flc.circFlux80, flc.fluxError80
jdates   = flc.jdates

# Open the file
f = 'table3_only_Ch1_Ch2.dat'
ref_mu = 18.91
dat = np.genfromtxt(f, dtype=(float,float, '|S10'),
                           delimiter   = '\t', 
                           skip_header = 43,
                           missing_values = '      ',
                           filling_values ='0',
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
refAbsMag36, refAbsMag45 = np.zeros((2, len(data)))
for i in range(len(data)):
    refAbsMag36[i] = np.asarray(data[i][0]) - ref_mu
    refAbsMag45[i] = np.asarray(data[i][1]) - ref_mu
refColor = refAbsMag36 - refAbsMag45

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
    elif any(s in data[i][2] for s in ('B5','B6','B7','B8','B9')) and\
    '[e]' not in data[i][2] and 'Be' not in data[i][2]:
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
    elif '[e]' in data[i][2]: 
        sgB_e.append(data[i])
    elif any(s in data[i][2] for s in ('XRB', 'Fe', 'extr', 'Be')):  # B extr?
        BeXray.append(data[i])
    elif any(s in data[i][2] for s in ('B0.5', 'B1')):
        earlyB.append(data[i])
    elif any(s in data[i][2] for s in ('B5', 'B6', 'B9')):
        lateB.append(data[i])
    elif any(s in data[i][2] for s in ('1', '2')):   # BC 1 and BC2 and BN?2
        earlyB.append(data[i])
    elif any(s in data[i][2] for s in ('BN6')):
        lateB.append(data[i])
    else:
        reject.append(data[i])

O_col, O_Mag           = np.zeros((2, 1, len(Ostar)))
earlyB_col, earlyB_Mag = np.zeros((2, 1, len(earlyB)))
lateB_col, lateB_Mag   = np.zeros((2, 1, len(lateB)))
AFG_col, AFG_Mag       = np.zeros((2, 1, len(AFG)))
RSG_col, RSG_Mag       = np.zeros((2, 1, len(RSG)))
WR_col, WR_Mag         = np.zeros((2, 1, len(WR)))
sgB_e_col, sgB_e_Mag   = np.zeros((2, 1, len(sgB_e)))
LBV_col, LBV_Mag       = np.zeros((2, 1, len(LBV)))
BeXray_col, BeXray_Mag = np.zeros((2, 1, len(BeXray)))

sortColors = [O_col, earlyB_col, lateB_col, AFG_col, RSG_col, WR_col, sgB_e_col,
             LBV_col, BeXray_col]
sortAbsMag36 = [O_Mag, earlyB_Mag, lateB_Mag, AFG_Mag, RSG_Mag, WR_Mag,
                   sgB_e_Mag, LBV_Mag, BeXray_Mag]
cat_list = (Ostar, earlyB, lateB, AFG, RSG, WR, sgB_e, LBV, BeXray)
# Calculate the x- and y-values for each star category
for cat in (Ostar, earlyB, lateB, AFG, RSG, WR, LBV, sgB_e, BeXray):
    ind = cat_list.index(cat)
    for i in range(len(cat)):
        sortColors[ind][0][i] = np.asarray(cat[i][0]) - np.asarray(cat[i][1])
        sortAbsMag36[ind][0][i] = np.asarray(cat[i][0]) - ref_mu
# Scatter plot the categories using unique shapes and colors   
markershapes = '^', 's', 'v', 'p', 'v', 'o', 'D', '<', 'o'
colors = 'k', 'silver','white', 'dimgrey', 'lightgrey', '#ffffff', \
         'g', 'deepskyblue', 'gainsboro'
edgecolorlist = 'White', 'None', 'k', 'darkgrey', 'k', 'k', 'k',\
                'k', 'k'
catnames = 'O stars', 'Early B', 'Late B', 'AFG/comp', 'RSG', 'WR','sgB[e]', 'LBV',\
           'Be-Xray'
linewidthlist = 0.3, 0.3, 0.3, 0.3, 0.3, 0.5, 0.3, 0.3, 1
sizelist = 20, 10, 20, 30, 12, 12, 25, 30, 45
# Plot results
fig, ax = plt.subplots()
plt.hold(True)
for cat in (earlyB, Ostar, WR, AFG, RSG, sgB_e, lateB, LBV, BeXray):
    ind = cat_list.index(cat)
    ax.scatter(sortColors[ind], sortAbsMag36[ind], c = colors[ind],
               marker = markershapes[ind], linewidths = linewidthlist[ind],
               s = sizelist[ind], edgecolors = edgecolorlist[ind],
               label = catnames[ind])
axGrab = plt.gca()
axGrab.set_ylim(ax.get_ylim()[::-1])  

# Add minor tick labels. 0.1 and 1 intervals
minorLocatorx = MultipleLocator(0.1)
minorLocatory = MultipleLocator(1)
ax.xaxis.set_minor_locator(minorLocatorx)
ax.yaxis.set_minor_locator(minorLocatory)

plt.legend(loc='lower right',fontsize = 10, scatterpoints = 1, frameon = False)
plt.xlabel('m$_{[3.6]}$ - m$_{[4.5]}$', fontsize = 14)
plt.ylabel('M$_{[3.6]}$', fontsize = 14)


# Handle the IC 10 X-2 fluxes
abs_Mag36, abs_Mag45, app_Mag36, app_Mag45, color = np.zeros((5,len(jdates)))
mu_X2 = 24.1    # Distance modulus of IC 10-X2
f0_36 = 280.9   # zero-point flux for the 3.6 channel
f0_45 = 179.7   # zero-point flux for the 4.5 channel
for j in range(len(jdates)):
    app_Mag36[j] = -2.5 * log10(f_36[j] * 10**-3 / f0_36)
    app_Mag45[j] = -2.5 * log10(f_45[j] * 10**-3 / f0_45)
    abs_Mag36[j] = app_Mag36[j] - mu_X2
    abs_Mag45[j] = app_Mag45[j] - mu_X2
    color[j]     = app_Mag36[j] - app_Mag45[j]

# Plot IC 10 X-2's journey across the CMD with a colormap
ic10x2 = plt.scatter(color, abs_Mag36, marker = 'o', c = jdates, 
                     cmap = plt.cm.plasma, linewidth = 0.3, s = 30,
                     edgecolor = 'k')
#sn2010da = plt.scatter(0.76, -10.36, marker = 's', c = 'red', linewidth = 0.3,
#                       s = 35, edgecolor = 'k')  
#ax.text(0.85,-10,'(SN)2010da')                       
# Add a colorbar
cb = plt.colorbar(ic10x2, shrink = 0.6,
                  ticks = [53500, 54375, 55250, 56125, 57000],
                  orientation = 'vertical')
cb.set_label('IC 10 X-2 \n Date (MJD)', y = 1.25, rotation = 0, fontsize = 14)
cb.ax.xaxis.set_label_position('top')
cb.ax.invert_yaxis()

# Make a zoomed inset: zoom-factor, location: upper center
axins = zoomed_inset_axes(ax, 3, loc = 1) 
ax.set_xlim(-0.5, 3)
ax.set_ylim(1, -16)
plt.gca().invert_yaxis()

for category in (earlyB, Ostar, WR, AFG, RSG, sgB_e, lateB, LBV, BeXray):
    ind = cat_list.index(category)
    axins.scatter(sortColors[ind], sortAbsMag36[ind], c = colors[ind],
                  marker = markershapes[ind], linewidth = linewidthlist[ind],
                  s = sizelist[ind], edgecolors = edgecolorlist[ind])
                  
axins.scatter(color, abs_Mag36, marker = 'o', c = jdates, cmap = plt.cm.plasma,
              s = 40, linewidth = 0.3, edgecolor = 'k')

# Create labels for points in the inset
labels = ['{0}'.format(i + 1) for i in range(len(jdates))]    
#for label, x, y in zip(labels, color, abs_Mag36):          
#    axins.annotate(label, xy = (x, y), xytext = (-20, 20),
#                   textcoords = 'offset points', ha='right', va='bottom',
#                   bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
#        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))   
# Don't know how to automate this so I'm just going to do it by hand
axins.annotate('1', xy = (color[0], abs_Mag36[0]), xytext = (17, -20),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.15', fc = 'yellow', alpha = 0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
axins.annotate('2', xy = (color[1], abs_Mag36[1]), xytext = (30, -20),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.15', fc = 'yellow', alpha = 0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))   
axins.annotate('3', xy = (color[2], abs_Mag36[2]), xytext = (-15, -30),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.15', fc = 'yellow', alpha = 0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))                
axins.annotate('4', xy = (color[3], abs_Mag36[3]), xytext = (-5, -30),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.15', fc = 'yellow', alpha = 0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
axins.annotate('5', xy = (color[4], abs_Mag36[4]), xytext = (20,12),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.15', fc = 'yellow', alpha = 0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
axins.annotate('6', xy = (color[5], abs_Mag36[5]), xytext = (7,12),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.15', fc = 'yellow', alpha = 0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
axins.annotate('7', xy = (color[6], abs_Mag36[6]), xytext = (-10,16),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.15', fc = 'yellow', alpha = 0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
axins.annotate('8', xy = (color[7], abs_Mag36[7]), xytext = (32,-30),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.15', fc = 'yellow', alpha = 0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
axins.annotate('9', xy = (color[8], abs_Mag36[8]), xytext = (-20,12),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.15', fc = 'yellow', alpha = 0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
axins.annotate('10', xy = (color[9], abs_Mag36[9]), xytext = (15,-30),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.15', fc = 'yellow', alpha = 0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
axins.annotate('11', xy = (color[10], abs_Mag36[10]), xytext = (2,-30),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.15', fc = 'yellow', alpha = 0.5),
                arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
                
#axins.text(0.20, -9.70, 'IC 10 X-2')
x1, x2, y1, y2 = -0.05, 0.38, -7.75, -9.55
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
plt.yticks(visible = False)
plt.xticks(visible = False)
# position the inset
mark_inset(ax, axins, loc1 = 2, loc2 = 4, fc = "none", ec = "0.5")
plt.savefig('170917 Bonanos CMD.pdf')
