import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.transforms as transforms
import matplotlib.colors as mcolors

from beamtools import upp, h, c, normalize

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
nalm_locs = locs[6:13][::-1]
num='000'
ext='.pkl'

main_file = [save_folder+'main'+A+num+ext for A in locs]
nalm_file = [save_folder+'nalm'+A+num+ext for A in nalm_locs]
out_file = save_folder+'outZ'+num+ext

#Plotting
#fiber lengths
smf1=0.72142857
smf2=1.44285714
smf3=5.0
ydf1=0.65

fiber = np.array([0,smf1,ydf1,smf1,smf2,smf2,0,smf1,ydf1,smf1,smf1,smf3,smf1,0,smf2,0,smf2])
length = np.cumsum(fiber)
nalm_fiber = np.array([0,smf1,smf3,smf1,smf1,ydf1,smf1])
nalm_length = length[6]+np.cumsum(nalm_fiber)

#sorting labels for pulses
main_ptlabels = ['A','C','F','I','L','P']
nalm_ptlabels = ['L','I']
plot_pts = [locs.index(i) for i in main_ptlabels]
nalm_pts = [nalm_locs.index(i) for i in nalm_ptlabels]
n2m_pts = [main_ptlabels.index(n) for n in [nalm_locs[i] for i in nalm_pts]]

#Label names for plot locations
main_ptnames = ['I','II','III','IV','V','VI']
nalm_ptnames = ["V'","IV'"]

#Noted fiber locations
fiber_sect = {
            3:{'name':'Yb fiber','ang':45}, 
            9:{'name':'Yb fiber','ang':45}, 
            11:{'name':'SMF','ang':0}
            }

#Fiber components (vertical lines)
fiber_comp = {0:'OC',6:'50:50',13:'50:50',15:'BPF',16:'OC'}

#all pulses
mp = [upp.load_obj(f) for f in main_file]
lp = [upp.load_obj(f) for f in nalm_file]

#pulses to individually plot
plot_mp = [mp[i] for i in plot_pts]
plot_lp = [lp[i] for i in nalm_pts]

#calculate RMS
t_rms = np.array([upp.rms_width(p.time,p.getPt())[1] for p in mp])
f_rms = np.array([upp.rms_width(p.freq,p.getPf())[1] for p in mp])
nt_rms = np.array([upp.rms_width(p.time,p.getPt())[1] for p in lp])
nf_rms = np.array([upp.rms_width(p.freq,p.getPf())[1] for p in lp])

hr=[0.2,4]+[1]*len(plot_pts)
plot_grid = [len(hr),2]
plot_w = 6.25
asp_ratio = 1.25
plot_h = plot_w*asp_ratio

c_yb = 'xkcd:greenish'
c_sm = 'xkcd:sea blue'
c_t = 'xkcd:ochre'
c_f = 'xkcd:eggplant'
c_c = 'xkcd:ocean blue'

fig = plt.figure(figsize=(plot_w,plot_h), facecolor='w')
gs = GridSpec(plot_grid[0], plot_grid[1], height_ratios=hr)

#Generate pulse width plot axes
ax1 = plt.subplot(gs[1,:])
ax2 = ax1.twinx()
a = ax1.get_position().bounds
[ax.set_position([a[0],a[1]+0.04,a[2],a[3]-0.04]) for ax in (ax1,ax2)]

ax1.plot(length,t_rms*1E12, c=c_t)
ax1.plot(nalm_length,nt_rms*1E12,'--',c=c_t)
ax1.set_xlabel('Fiber length (m)')
ax1.set_ylabel('Pulse width RMS (ps)')
ax1.set_ylim([0.6,1.8])

ax2.plot(length,f_rms*1E-12,c=c_f)
ax2.plot(nalm_length,nf_rms*1E-12,'--',c=c_f)
ax2.set_ylabel('Spectral width RMS (ps$^{-1}$)')
ax2.set_ylim([0,5])

#Fiber section highlights
ax1.axvspan(length[1], length[2], facecolor=c_yb, alpha=0.25)
ax1.axvspan(length[7], length[8], 0.5,1, facecolor=c_yb, alpha=0.25)
ax1.axvspan(length[10], length[11], 0.5,1, facecolor=c_sm, alpha=0.25)
ax1.axvspan(nalm_length[1], nalm_length[2], 0,0.5, 
            facecolor=c_sm, alpha=0.25, hatch='\\\\',
            edgecolor='k', linewidth=0)
ax1.axvspan(nalm_length[4], nalm_length[5], 0,0.5, 
            facecolor=c_yb, alpha=0.25, hatch='\\\\',
            edgecolor='k', linewidth=0)

#Plot location markers
bbox_props = dict(boxstyle='circle', fc='xkcd:cool grey', 
                ec='xkcd:dark grey', lw=1, alpha=0.5)

trans = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
[ax1.text(length[pt],0.94, main_ptnames[i], 
        ha='center', va='center', 
        bbox=bbox_props, 
        transform=trans) for i,pt in enumerate(plot_pts)]
[ax1.text(nalm_length[pt],0.03, nalm_ptnames[i], 
        ha='center', va='center', 
        bbox=bbox_props,
        transform=trans) for i,pt in enumerate(nalm_pts)]

#Generate top plot
ax0 = plt.subplot(gs[0,:], sharex=ax1)
a = ax0.get_position().bounds
ax0.set_position([a[0],a[1]-0.01,a[2],a[3]+0.01])
[ax0.spines[k].set_visible(False) for k in ax0.spines]
ax0.set_xticklabels([])
ax0.set_yticklabels([])
ax0.set_xticks([])
ax0.set_yticks([])

#Fiber section highlights
ax0.axvspan(length[1], length[2], facecolor=c_yb, alpha=0.25)
ax0.axvspan(length[7], length[8], facecolor=c_yb, alpha=0.25)
ax0.axvspan(length[10], length[11], facecolor=c_sm, alpha=0.25)

#Fiber component top labels
trans = transforms.blended_transform_factory(ax0.transData, ax0.transAxes)
[(ax1.axvline(length[k], color='xkcd:dark grey'),
    ax0.text(length[k],0, v, 
        transform = trans, 
        va='bottom',rotation=45)) for k,v in fiber_comp.items()]

#Section top labels
[ax0.text((length[k]+length[k-1])/2,0, v['name'], 
        transform = trans, 
        va='bottom',ha='center',
        rotation=v['ang']) for k,v in fiber_sect.items()]

#Generate pulse plots
at=[]
af=[]
ac=[]

#Main pulse plots
for i, p in enumerate(plot_mp):
    at.append(plt.subplot(gs[(i+2),0]))
    ac.append(at[i].twinx())
    at[i].plot(p.time*1E12,normalize(p.getPt()),c=c_t,zorder=2)
    ac[i].plot(p.time*1E12,p.chirp()*1E-12,'-',
        c=c_c,linewidth=0.75,zorder=1)
    at[i].set_xlim([-4,4])

    af.append(plt.subplot(gs[(i+2),1]))
    af[i].plot(np.fft.fftshift(p.freq*1E-12),
        np.fft.fftshift(normalize(p.getPf())),c=c_f)
    af[i].set_xlim([-12,12])
    af[i].set_yticklabels([])

    ac[i].set_ylim([-8,8])
    ac[i].set_yticklabels([])

    af[i].text(-0.1,0.5, main_ptnames[i], 
        ha="center", transform=af[i].transAxes, bbox=bbox_props)

    if i < len(plot_mp)-1:
        at[i].set_xticklabels([])
        af[i].set_xticklabels([])

    else:
        at[i].set_xlabel('Tau (ps)')
        af[i].set_xlabel('Omega (ps$^{-1}$)')

#Nalm pulse plots
for i, p in enumerate(plot_lp):
    at[n2m_pts[i]].plot(p.time*1E12,normalize(p.getPt()),'--',
        c=c_t,zorder=2, linewidth=1.0)
    ac[n2m_pts[i]].plot(p.time*1E12,p.chirp()*1E-12,'--', 
        c=c_c,linewidth=0.75,zorder=1)
    af[n2m_pts[i]].plot(np.fft.fftshift(p.freq*1E-12),
        np.fft.fftshift(normalize(p.getPf())),'--',c=c_f,linewidth=1.0)
   
plt.show()




