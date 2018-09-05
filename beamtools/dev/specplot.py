#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 8 17:47:28 2018

@author: cpkmanchee
"""
import numpy as np
import matplotlib.pyplot as plt
import beamtools as bt
from matplotlib.gridspec import GridSpec

show_plt = True
save_plt = False

dpi=600

wlim = [1005,1065]

f_sp ='/Users/cpkmanchee/Google Drive/PhD/Data/2016-05-13 REGEN spectra/2016-05-18 REGEN Cavity Dumped/20160518-002.csv'

sp = bt.import_data_file(f_sp, 'oo_spec')
sp.intensity = bt.normalize(sp.intensity)

# flimit, _ = bt.pulse.spectrumFT([spec.wavelength,spec.intensity])
# AiFT = bt.pulse.autocorr(np.abs(flimit.et)**2)

plot_grid = [1,1]
plot_w = 3
asp_ratio = 0.8*plot_grid[0]/plot_grid[1]
plot_h = plot_w*asp_ratio

#create grid
fig = plt.figure(figsize=(plot_w,plot_h), facecolor='w')
gs = GridSpec(plot_grid[0], plot_grid[1])

ax1 = fig.add_subplot(gs[0,0])
#plot spectrum
ax1.plot(sp.wavelength, sp.intensity, '-', c='xkcd:brick red')
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Intensity (arb.)')
#ax1.set_xlim(wlim)
ax1.tick_params('both', labelsize='x-small')


fig.tight_layout()
if show_plt:
    plt.show()
    
if save_plt:
    fig.savefig('spec.png', dpi=dpi, bbox_inches='tight')
