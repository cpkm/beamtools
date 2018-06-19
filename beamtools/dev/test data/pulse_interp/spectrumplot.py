#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 8 17:47:28 2018

@author: cpkmanchee
"""
import numpy as np
import matplotlib.pyplot as plt
import beamtools as bt

show_plt = True
save_plt = False

dpi=600

wlim = [1020,1040]
tlim = [-10,20]


file = '20180103-004-000ac001.txt'

ac = bt.import_data_file(file, 'bt_ac')
ac.intensity = ac.power
ac.power = bt.normalize(ac.power)

ac_fit,_ = bt.pulse.fit_ac([ac.delay,ac.power], bgform='linear')


file = '20180103-004-000.csv'
spec = bt.import_data_file(file, 'oo_spec')
spec.intensity = bt.normalize(spec.intensity)

flimit, _ = bt.pulse.spectrumFT([spec.wavelength,spec.intensity])

AiFT = bt.pulse.autocorr(np.abs(flimit.et)**2)

pfit, _ = bt.pulse.fit_ac([flimit.time,AiFT])

fig = plt.figure(figsize=(4,3.25), facecolor='w')
#ax2 = ax1.twinx()
#ax2.scatter(spec.pump, spec.energy, s=0)
#ax2.set_ylabel('Pulse energy (uJ)')
#ax2.tick_params('y', colors='k')

#ax1.text(20,65, '{:.0f}% slope \nefficiency'.format(slope*100))

#scale = np.around(spec.energy.values[-1]/spec.out.values[-1], decimals=2)

#ax1.set_ylim([0,np.around(spec.intensity.max(),decimals=-1)])
#ax2.set_ylim(np.asarray(ax1.get_ylim())*scale)

f2p = 'sech2'

ax2 = fig.add_subplot(111)

for fit in pfit:
    if fit.ftype.lower() in bt.alias_dict[f2p]:
        ax2.plot(flimit.time*1E12,bt.normalize(AiFT), '-', 
        label='TL pulse', c='xkcd:ocean blue')

        label = 'TL - {} fit\n{:.2f} ps'.format(fit.ftype, bt.pulse.sigma_fwhm(fit.popt[0]*1E12,fit.ftype))
        ax2.plot(flimit.time*1E12, bt.normalize(fit.subs(flimit.time)), '--',
            label=label, c='xkcd:deep purple')

for fit in ac_fit:
    if fit.ftype.lower() in bt.alias_dict[f2p]:
        ax2.plot(ac.delay-fit.popt[2],bt.normalize(bt.rmbg([ac.delay,ac.power],fit=fit)), '-', 
        label='Autocorrelation', c='xkcd:chartreuse')

        label = 'AC - {} fit\n{:.2f} ps'.format(fit.ftype, bt.pulse.sigma_fwhm(fit.popt[0],fit.ftype))
        ax2.plot(ac.delay-fit.popt[2], bt.normalize(bt.rmbg([ac.delay,fit.subs(ac.delay)],fit=fit)), '--',
            label=label, c='xkcd:dark green')

ax2.set_xlim(tlim)
ax2.set_xlabel('Delay (ps)')
ax2.set_ylabel('Intensity (arb.)')
ax2.legend(loc=[0.55,0.1], fontsize='small')
fig.tight_layout()

ax1 = plt.axes([.66, .64, .25, .25], facecolor='y')
ax1.plot(spec.wavelength, spec.intensity, '-', c='xkcd:ocean blue')
#ax1.set_xlabel('Wavelength (nm)')
# Make the y-axis label, ticks and tick labels match the line color.
#ax1.set_ylabel('Intensity')
ax1.set_xticks([1020,1030,1040])
ax1.tick_params('both', labelsize='x-small')
ax1.set_xlim(wlim)

if show_plt:
    plt.show()
    
if save_plt:
    fig.savefig('spectrum.png', dpi=dpi, bbox_inches='tight')
