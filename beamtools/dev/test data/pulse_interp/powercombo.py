#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 12:40:28 2017

@author: cpkmanchee
"""
import numpy as np
import beamtools as bt
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from matplotlib.gridspec import GridSpec

show_plt = True
save_plt = False

dpi = 600

source_dir = '/Users/cpkmanchee/Google Drive/PhD/Data/'

#import time-power data
file = source_dir + '2018-01-03 NALM Oscillator/2018-01-12 PM_LogData.txt'
pdata = bt.import_data_file(file,'thorlabs_pm')
tformat = bt.file_formats[pdata.filetype]['time_format']
time = [datetime.strptime(d,tformat) for d in pdata.time]
time = np.array([(t-time[0]).total_seconds() for t in time])

s=8000
f=10400

t = time[s:f]-time[s]
p = pdata.power[s:f]
p0 = np.mean(p)
std = np.std(p)

pfit = np.poly1d(np.polyfit(t,p,24))
res = p-pfit(t)
perr = np.sqrt(np.sum(res**2)/(len(res)-1))


#set up plot grid
plot_grid = [1,1]
plot_w = 3.25
asp_ratio = 0.67*plot_grid[0]/plot_grid[1]
plot_h = plot_w*asp_ratio

fig = plt.figure(figsize=(plot_w,plot_h), facecolor='w')
gs = GridSpec(plot_grid[0], plot_grid[1])

# identical to ax1 = plt.subplot(gs.new_subplotspec((0,0), colspan=3))
ax3 = plt.subplot(gs[0,0])

#plot time-power
ax3.plot(t,p/p0,'.', c='xkcd:ocean blue', marker='.',markersize=2)
ax3.plot(t,pfit(t)/p0)
ax3.set_xlabel('Time (sec.)')
# Make the y-axis label, ticks and tick labels match the line color.
ax3.set_ylabel('Power $(P/P_{avg})$')
ax3.tick_params('y', colors='k')

t1 = '$P_{avg} = $' + '{:.1f} mW'.format(1000*p0)
t2 = '\n$\sigma_{total}$ = ' + '{:.2f}%'.format(100*std/p0) 
t3 = '\n$\sigma_{high}$ = ' + '{:.2f}%'.format(100*perr/p0) 
text = t1+t2+t3

ax3.text(1000,.992, text, fontsize='small')
ax3.set_ylim([0.99,1.005])

gs.tight_layout(fig)

if show_plt:
    plt.show()
    
if save_plt:
    fig.savefig('power.png', dpi=dpi, bbox_inches='tight')
