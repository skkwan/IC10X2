# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 14:30:51 2017

Perform Lomb-Scargle algorithm on PTF optical data.

Dependencies: PTF_fulltable_phot{1/2}.txt

@author: stephaniekwan
"""

# Author: Jake VanderPlas
# License: BSD
#   The figure produced by this code is published in the textbook
#   "Statistics, Data Mining, and Machine Learning in Astronomy" (2013)
#   For more information, see http://astroML.github.com
#   To report a bug or issue, use the following forum:
#    https://groups.google.com/forum/#!forum/astroml-general
import numpy as np
from matplotlib import pyplot as plt
from astroML.time_series import\
    lomb_scargle, lomb_scargle_BIC, lomb_scargle_bootstrap

#----------------------------------------------------------------------
# This function adjusts matplotlib settings for a uniform feel in the textbook.
# Note that with usetex=True, fonts are rendered with LaTeX.  This may
# result in an error if LaTeX is not installed on your system.  In that case,
# you can set usetex to False.
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize = 8, usetex=True)

#------------------------------------------------------------
# Generate Data
table1 = np.genfromtxt('PTF_fulltable_phot1.txt', 
                  comments = '|', skip_header = 43, skip_footer = 2)
table2 = np.genfromtxt('PTF_fulltable_phot2.txt',
                       comments = '|', skip_header = 43, skip_footer = 2)                  
                  
flux1, flux2 = table1[:, 7], table2[:, 7] 

sigflux1, sigflux2 = table1[:, 8], table2[:, 8]
mjd1, mjd2 = table1[:, 15], table2[:, 15]
snr1, snr2 = table1[:, 9], table2[:, 9]
zp1, zp2 = table1[:, 5], table2[:, 5]

snt, snu = 3, 5
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
        
#------------------------------------------------------------
# Compute periodogram
period = np.linspace(1, 50, 100000)
omega = 2 * np.pi / period
#PS = lomb_scargle(t, y_obs, dy, omega, generalized=True)

PS = lomb_scargle(mag1date, mag1, mag1sig, omega, generalized = True)

#------------------------------------------------------------
# Get significance via bootstrap
D = lomb_scargle_bootstrap(mag1date, mag1, mag1sig, omega, generalized = True,
                           N_bootstraps=500, random_state=0)
sig1, sig5 = np.percentile(D, [99, 95])

#------------------------------------------------------------
# Plot the results
fig = plt.figure(figsize=(5, 3.75))
fig.subplots_adjust(left=0.1, right=0.9, hspace=0.25)

# First panel: the data
ax = fig.add_subplot(211)
ax.errorbar(mag1date, mag1, mag1sig, fmt='.k', lw=1, ecolor='gray')
ax.set_xlabel('time (days)')
ax.set_ylabel('flux')
ax.set_xlim(mag1date[0] - 50, mag1date[-1] + 50)
ax.set_ylim(16, max(mag1) + 0.5)

# Second panel: the periodogram & significance levels
ax1 = fig.add_subplot(212, xscale='log')
ax1.plot(period, PS, '-', c='black', lw=1, zorder=1)
ax1.plot([period[0], period[-1]], [sig1, sig1], ':', c='black')
ax1.plot([period[0], period[-1]], [sig5, sig5], ':', c='black')

#ax1.annotate("", (0.3, 0.65), (0.3, 0.85), ha='center',
#             arrowprops=dict(arrowstyle='->'))

ax1.set_xlim(period[0], period[-1])
ax1.set_ylim(-0.05, 0.85)

ax1.set_xlabel(r'period (days)')
ax1.set_ylabel('power')

# Twin axis: label BIC on the right side
ax2 = ax1.twinx()
ax2.set_ylim(tuple(lomb_scargle_BIC(ax1.get_ylim(), mag1, mag1sig)))
ax2.set_ylabel(r'$\Delta BIC$')

ax1.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
ax1.xaxis.set_minor_formatter(plt.FormatStrFormatter('%.1f'))
ax1.xaxis.set_major_locator(plt.LogLocator(10))
ax1.xaxis.set_major_formatter(plt.FormatStrFormatter('%.3g'))

plt.savefig("042817_Lomb_scargle_R_band.pdf")

