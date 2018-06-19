#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 19:49:28 2018

@author: cpkmanchee
"""
import numpy as np
import matplotlib.pyplot as plt
import beamtools as bt

from matplotlib.gridspec import GridSpec

show_plt = True
save_plt = False

dpi=600

wlim = [1000,1100]

spec_file = '/Users/cpkmanchee/Documents/Code/beamtools/beamtools/dev/test data/pulse_interp/20180103-004-000.csv'
ac_file = '/Users/cpkmanchee/Documents/Code/beamtools/beamtools/dev/test data/pulse_interp/20180103-004-000ac001.txt'

sp = bt.import_data_file(spec_file, 'oo_spec')
ac = bt.import_data_file(ac_file, 'bt_ac')

plot_grid = [1,1]
plot_w = 4.50
asp_ratio = 1.25*plot_grid[0]/plot_grid[1]
plot_h = plot_w*asp_ratio

fig = plt.figure(figsize=(plot_w,plot_h), facecolor='w')
gs = GridSpec(plot_grid[0], plot_grid[1])

ax1 = fig.add_subplot(gs[0,0])


ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('Intensity (arb.)')
ax1.set_xlim(wlim)
ax1.tick_params('both')

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
# fake up the array of the scalar mappable. Urgh...
sm.set_array([])
cbx = plt.colorbar(sm,ticks=[0,1])

cbx.ax.set_yticklabels(['40 W','210 W'])
cbx.ax.set_ylabel('Pump power',labelpad=-20)

fig.tight_layout()

if show_plt:
    plt.show()
    
if save_plt:
    fig.savefig('rod-spectra.png', dpi=dpi, bbox_inches='tight')
