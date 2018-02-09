# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 21:13:47 2017

@author: cpkmanchee

Simulate pulse propagation in oscillator

Schematic:



Notes:
This file requires pulsemodel.py
This file uses the functions and classes defined in pulsemodel.py (used via 
import)

Everything is in SI units: m,s,W,J etc.
"""

import numpy as np
import matplotlib.pyplot as plt

from beamtools import upp, h, c

from tqdm import tqdm, trange

import sys
import shutil
import glob
import os
from datetime import datetime


#create save name and folders for output
start_date = datetime.now().strftime("%Y%m%d")
start_time = datetime.now().strftime("%H%M%S")

#output_folder is outside of git repository in: code_folder/Code Output/...
#if file is in code_folder/X/Y/Z, results are in code_folder/Code Output/X/Y/Z
code_folder = '/Users/cpkmanchee/Documents/Code'
output_folder = (code_folder + '/Code Output' 
                + os.path.dirname(__file__).split(code_folder)[-1] + '/' 
                + os.path.splitext(os.path.basename(__file__))[0] + '_output')

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

result_folder = (output_folder + '/' 
                + start_date + os.path.splitext(os.path.basename(__file__))[0])

dataset_num = 0
while not not glob.glob((result_folder+'-'+str(dataset_num).zfill(2)+'*')):
    dataset_num = dataset_num + 1
result_folder =  result_folder + '-' + str(dataset_num).zfill(2)
os.makedirs(result_folder)

filebase = (result_folder + '/' + start_date + '-' 
            + start_time + '-' + str(dataset_num).zfill(2))
fileext = '.pkl'
output_num = 0
filename =  filebase + 'pulse' + str(output_num).zfill(3) + fileext

shutil.copy(__file__, result_folder + '/' + os.path.basename(__file__))


def savepulse(pulse, name='pulse'):
    '''
    to be used locally only
    all file/folder names must be previously defined
    '''
    global output_num, filename    
    
    while not not glob.glob(filename):
        output_num = output_num + 1
        filename = filebase + name + str(output_num).zfill(3) + fileext
    upp.save_obj(pulse,filename)


def cavity(pulse,auto_z_step=False):
    '''Define cavity round trip
    NOTE: pulse object is modified!!!
    returns:
        pulse.At = current pulse profile
        output_At = cavity output (outcoupled) profile
    '''
    pulse.At = upp.propagate_fiber(pulse,smf1,autodz=auto_z_step)

    Ps = np.sum(np.abs(pulse.At)**2)*pulse.dt/tau_rt
    ydf1.gain = upp.calc_gain(ydf1,p1P,Ps)
    pulse.At = upp.propagate_fiber(pulse,ydf1,autodz=False)

    pulse.At = upp.propagate_fiber(pulse,smf1,autodz=auto_z_step)
    pulse.At,_ = upp.power_tap(pulse, 0, loss=0.06)
    pulse.At = upp.propagate_fiber(pulse,smf2,autodz=auto_z_step)
    pulse.At,_ = upp.power_tap(pulse, 0, loss=0.06)
    pulse.At = upp.propagate_fiber(pulse,smf2,autodz=auto_z_step)
    
    nalmp = pulse.copyPulse()
    pulse.At, nalmp.At = upp.coupler_2x2(pulse, None, tap=50, loss=0.06)

    Psn = np.sum(np.abs(pulse.At)**2)*pulse.dt/tau_rt
    ydf2.gain = upp.calc_gain(ydf2,p2P,Psn)

    #main pulse
    pulse.At = upp.propagate_fiber(pulse,smf1,autodz=auto_z_step)
    pulse.At = upp.propagate_fiber(pulse,ydf2,autodz=False)
    pulse.At = upp.propagate_fiber(pulse,smf1,autodz=auto_z_step)
    pulse.At,_ = upp.power_tap(pulse, 0, loss=0.06)
    pulse.At = upp.propagate_fiber(pulse,smf1,autodz=auto_z_step)
    pulse.At = upp.propagate_fiber(pulse,smf3,autodz=auto_z_step)

    #nalmp
    nalmp.At = upp.propagate_fiber(nalmp,smf1,autodz=auto_z_step)
    nalmp.At = upp.propagate_fiber(nalmp,smf3,autodz=auto_z_step)
    nalmp.At = upp.propagate_fiber(nalmp,smf1,autodz=auto_z_step)
    nalmp.At,_ = upp.power_tap(nalmp, 0, loss=0.06)
    nalmp.At = upp.propagate_fiber(nalmp,smf1,autodz=auto_z_step)
    nalmp.At = upp.propagate_fiber(nalmp,ydf2,autodz=False)
    nalmp.At = upp.propagate_fiber(nalmp,smf1,autodz=auto_z_step)

    pulse.At,_  = upp.coupler_2x2(pulse, nalmp, tap=50, loss=0.06)
    pulse.At = upp.propagate_fiber(pulse,smf2,autodz=auto_z_step)
    pulse.At = upp.optical_filter(pulse, filter_type='bpf', bandwidth=2E-9, loss=0.06, order=2)
    pulse.At = upp.propagate_fiber(pulse,smf2,autodz=auto_z_step)

    pulse.At, output_At = upp.coupler_2x2(pulse, None, tap=75, loss=0.06)

    return pulse.At, output_At


def run_sim(
        pulse, max_iter=100, err_thresh=1E-6, 
        save_pulse = 1, auto_z_step=False):

    if save_pulse is not False:
        N = min(max_iter,save_pulse)
        savepulse(pulse,name='cavity')
    else:
        N = max_iter + 1

    f,ax = initialize_plot(10)
    update_pulse_plot(pulse,f,ax)

    t = trange(max_iter, desc='Total progress')
    t.set_postfix(str='{:.1e}'.format(0))
    for i in t:
        input_At = pulse.At
        cavity_At, output_At = cavity(pulse, auto_z_step)

        update_pulse_plot(pulse,f,ax)

        if (i+1)%N == 0:
            savepulse(pulse, name='cavity')
            savepulse(pulse.copyPulse(output_At), name='output')

        power_in = np.abs(input_At)**2
        power_out = np.abs(pulse.At)**2

        test = check_residuals(power_in,power_out,
            integ_err=err_thresh, p2p_err=err_thresh)

        if test[0]:
            if save_pulse is not False:
                savepulse(pulse,name='cavity')
                savepulse(pulse.copyPulse(output_At), name='output')

            break

        t.set_postfix(str='{:.1e},{:.1e}'.format(test[1],test[2]))

    if save_pulse is not False:
        savepulse(pulse,name='cavity')
        savepulse(pulse.copyPulse(output_At), name='output')


def check_residuals(initial, final, integ_err=1E-4, p2p_err=1E-4):
    '''Check residuals for covergence test.
    Return True if pass. False if fail.
    '''
    res = (initial-final)
    p2p = np.abs(res).max()/initial.max()
    integ = (np.sum(np.abs(res)**2)**(1/2))/np.sum(initial)

    if p2p < p2p_err and integ < integ_err:
        return True,integ,p2p
    else:
        return False,integ,p2p


def initialize_plot(N,colormap=plt.cm.viridis):
    '''Set up plotting figure for sim.'''
    plt.ion()
    cm = colormap

    f, ax = plt.subplots(1,2)
    [[a.plot([],[],c=cm(i/(N-1)),zorder=i)[0] for _,i in enumerate(range(N))] for a in ax]
    ax[0].set_xlim([-2E-11,2E-11])
    ax[0].set_xlabel('Tau (s)')
    ax[1].set_xlim([-5E12,5E12])
    ax[1].set_xlabel('Omega ($s^{-1}$)')

    plt.show()
    plt.pause(0.001)

    return f,ax

def update_pulse_plot(pulse,fig,ax):
    '''Update pulse plot'''
    #Update time plots
    lines = ax[0].lines
    N = len(lines)
    i = np.argmin([li.zorder for li in lines])
    lines[i].set_data(pulse.time,pulse.getPt())
    lines[i].zorder += N
    ax[0].set_ylim([0,np.max(pulse.getPt())])
    plt.pause(0.001)

    lines = ax[1].lines
    N = len(lines)
    i = np.argmin([li.zorder for li in lines])
    lines[i].set_data(np.fft.fftshift(pulse.freq),np.fft.fftshift(pulse.getPf()))
    lines[i].zorder += N
    ax[1].set_ylim([0,np.max(pulse.getPf())])
    plt.pause(0.001)



#Define Pulse Object
pulse = upp.Pulse(1.03E-6)
pulse.initializeGrid(16, 1.5E-9)
T0 = 200E-15
mshape = 1
chirp0 = 0
P_peak = 1E1   #peak power, 10kW-->1ps pulse, 400mW avg @ 40MHz
pulse.At = np.sqrt(P_peak)*(
            np.exp(-((1+1j*chirp0)/(2*T0**2))*pulse.time**(2*mshape)))

#Define fiber components
smf1 = upp.Fiber(0.65)
smf1.alpha = 0.000576
smf1.beta = np.array([
            0.0251222977, 
            4.5522276126132602e-05, 
            -5.0542788517531417e-08])*(1E-12)**(np.array([2,3,4]))
smf1.gamma = 0.00045
smf1.core_d = 5.5E-6

smf2 = smf1.copyFiber(length=1.3)
smf3 = smf1.copyFiber(length=5.0)


#gain fiber, nufern ysf-HI
ydf1 = upp.FiberGain(0.65, grid_type='rel',z_grid=100)
ydf1.alpha = 0.00345
ydf1.beta = np.array([
            0.0251222977, 
            4.5522276126132602e-05, 
            -5.0542788517531417e-08])*(1E-12)**(np.array([2,3,4]))
ydf1.gamma = 0.00045
ydf1.sigma_a = np.array([3.04306,0.04966])*1E-24
ydf1.sigma_e = np.array([3.17025,0.59601])*1E-24
ydf1.lambdas = np.array([0.976,1.030])*1E-6
ydf1.core_d = 6.0E-6
ydf1.N = 1.891669E25

ydf2 = ydf1.copyFiber()

#Pump parameters
pP = 0.25    #pump power, CW
Rp = 0.5
p1P = Rp*pP
p2P = (1-Rp)*pP

#Cavity
tau_rt = 1/(13.1E6)
