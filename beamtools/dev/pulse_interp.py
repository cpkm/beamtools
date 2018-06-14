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

wlim = [1010,1060]

spectrum_file = '/Users/cpkmanchee/Google Drive/PhD/Data/2017-12-13 Rod power output/'
file_pre =  source_dir+'spectra/'+'20171213-rod_'
ext = '.csv'

files = ['40W','60W','100W','120W','160W','180W','220W']


plot_grid = [1,1]
plot_w = 4.50
asp_ratio = 1.25*plot_grid[0]/plot_grid[1]
plot_h = plot_w*asp_ratio

fig = plt.figure(figsize=(plot_w,plot_h), facecolor='w')
gs = GridSpec(plot_grid[0], plot_grid[1])

ax1 = fig.add_subplot(gs[0,0])

cmap = LinearSegmentedColormap.from_list('viridis_half',plt.cm.viridis(np.linspace(0.25,.75,256)))
N=len(files)

for i,f in enumerate(files):
    file = file_pre+f+ext
    spec = bt.import_data_file(file, 'oo_spec')
    spcr,intr = roi(spec.wavelength, spec.intensity,wlim)
    intr = bt.normalize(intr, method='area')
    a,v = bt.moments(spcr,intr)
    s = v**(1/2)
    ax1.plot(spcr, intr+.03*i, '-', c=cmap(i/(N-1)))
    ax1.annotate('{:.0f} nm\n{:.1f} nm'.format(a,4*s), xy=(1050,i*0.03+.002),
           va = "bottom", ha="left", fontsize = 'small')

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
