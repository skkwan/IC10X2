# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 10:28:45 2016

@author: stephaniekwan

Interpolate astronomical silicate emissivity to create it as a function of
wavelength.

Updated August 9th to create it as a function of frequency.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.constants import h,k,c

# Open the file
f = open('AySi_emissivity.dat', 'r')
# Read all the lines into an array
lines = f.readlines()
# Save headers for quick reference
date, dustname, authors, header, header2, dustspec, colnames = \
    lines[0], lines[1], lines[2], lines[3], lines[4], lines[6], lines[7]
# Close the file
f.close

# Initialize arrays for emissivity coefficients and corresponding wavelengths 
wl   = np.zeros(len(lines) - 8)       # microns
q_abs = np.zeros(len(lines) - 8)     # no units
# Loop through "lines" and save values, skipping the first 8 lines
for i in range(8, len(lines)):
    wl[i - 8]   = lines[i].split()[0]
    q_abs[i - 8] = lines[i].split()[1]
fr = c / (wl * 1e-6)
wl = wl[::-1]
q_abs_for_wl = q_abs[::-1]


# Interpolate the function qWav(wl) 
qWav = interp1d(wl, q_abs_for_wl)
qFrq = interp1d(fr, q_abs)
# Plot results
wlnew   = np.arange(1.00000e-3, 10.00000000e+01, 1e-5)
frnew   = c / (wlnew * 1e-6)
fig = plt.figure()
plt.semilogy(fr, qFrq(fr), 'o', frnew, qFrq(frnew), '-')
# Label the axes
plt.xlabel('Frequency (Hz)', fontsize=13)
plt.ylabel('Q(freq)', fontsize=13)
fig.suptitle('Interpolation of emissivity of 0.1 um astronomical silicate dust (Draine)')
# Scale the axes to show important features
#plt.xlim([0.0, 1e17])
#plt.xlim(0.0, 2e14)
#plt.ylim(0.0,0.002)
#plt.xlim(0.1,20.00)
#plt.clf()
#plt.close()
plt.savefig("160809 Astronomical silicate Q interpolation function of freq.pdf")
