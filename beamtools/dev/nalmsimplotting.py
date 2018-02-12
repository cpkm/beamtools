import numpy as np
import matplotlib.pyplot as plt

from beamtools import upp, h, c

from tqdm import tqdm, trange

import sys
import shutil
import glob
import os
from datetime import datetime

from string import ascii_lowercase as abet

folder = ('/Users/cpkmanchee/Documents/Code/Code Output/'
            'beamtools/beamtools/dev/nalmoscillator_output/20180209nalmoscillator-02/')
save_folder = folder + 'single pass/'

locs = abet[:17].upper()
num='000'
ext='.pkl'

main_file = [save_folder+'main'+A+num+ext for A in locs]
nalm_file = [save_folder+'nalm'+A+num+ext for A in locs[6:13]]
out_file = save_folder+'outZ'+num+ext

smf1=0.72142857
smf2=1.44285714
smf3=5.0
ydf1=0.65

fiber = np.array([0,smf1,ydf1,smf1,smf2,smf2,0,smf1,ydf1,smf1,smf1,smf3,smf1,0,smf2,0,smf2])
length = np.cumsum(fiber)

mp = [upp.load_obj(f) for f in main_file]

t_rms = [upp.rms_width(p.time,p.getPt())[1] for p in mp]
f_rms = [upp.rms_width(p.freq,p.getPf())[1] for p in mp]

plot_grid = [2,len(main_file)]
plot_w = 12.5
asp_ratio = plot_grid[0]/plot_grid[1]
plot_h = plot_w*asp_ratio


f, ax = plt.subplots(2,len(main_file), sharey='row',figsize=(plot_w,plot_h))
plt.show()

f, ax = plt.subplots(2, sharex=True)
ax[0].plot(length,t_rms)
ax[1].plot(length,f_rms)
plt.show()
