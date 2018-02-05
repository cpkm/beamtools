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
import scipy as sp
import matplotlib.pyplot as plt

from beamtools import upp

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
output_folder = code_folder + '/Code Output' 
                + os.path.dirname(__file__).split(code_folder)[-1] + '/' 
                + os.path.splitext(os.path.basename(__file__))[0] + '_output'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

result_folder = output_folder + '/' 
                + start_date + os.path.splitext(os.path.basename(__file__))[0]

dataset_num = 0
while not not glob.glob((result_folder+'-'+str(dataset_num).zfill(2)+'*')):
    dataset_num = dataset_num + 1
result_folder =  result_folder + '-' + str(dataset_num).zfill(2)
os.makedirs(result_folder)

filebase = result_folder + '/' + start_date + '-' 
            + start_time + '-' + str(dataset_num).zfill(2)
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
    upp.saveObj(pulse,filename)


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

    pulse.At = upp.propagate_fiber(pulse,smf2,autodz=auto_z_step)

    pulse.At, output_At = upp.coupler_2x2(pulse,None,tap=25)

    return pulse.At, output_At


def run_sim(
        pulse, max_iter=100, err_thresh=1E-6, 
        save_pulse = 1, auto_z_step=False):

    if save_pulse is not False:
        N = min(max_iter,save_pulse)
        savepulse(pulse,name='cavity')
    else:
        N = max_iter + 1

    t = trange(max_iter, desc='Total progress')
    t.set_postfix(str='{:.1e}'.format(0))
    for i in t:
        input_At = pulse.At
        cavity_At, output_At = cavity(pulse, auto_z_step)

        if (i+1)%N == 0
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


#constants
h = upp.h  #J*s
c = upp.c     #m/s

#Define Pulse Object
pulse = upp.Pulse(1.03E-6)
pulse.initializeGrid(18, 1.5E-9)
T0 = 1000E-15
mshape = 1
chirp0 = 0
P_peak = 10E3   #peak power, 10kW-->1ps pulse, 400mW avg @ 40MHz
pulse.At = np.sqrt(P_peak)*(
            sp.exp(-(1/(2*T0**2))*(1+1j*chirp0)*pulse.time**(2*mshape)))

input_pulse = pulse.copyPulse()

#Define fiber components
smf1 = upp.Fiber(2.0)
smf1.alpha = 0.000576
smf1.beta = np.array([
            0.0251222977, 
            4.5522276126132602e-05, 
            -5.0542788517531417e-08])*(1E-12)**(np.array([2,3,4]))
smf1.gamma = 0.00045
smf1.core_d = 5.5E-6

smf2 = upp.Fiber(1.0)
smf2.alpha = 0.000576
smf2.beta = np.array([
            0.0251222977, 
            4.5522276126132602e-05, 
            -5.0542788517531417e-08])*(1E-12)**(np.array([2,3,4]))
smf2.gamma = 0.00045
smf2.core_d = 5.5E-6

smf3 = upp.Fiber(1.0)
smf3.alpha = 0.000576
smf3.beta = np.array([
            0.0251222977, 
            4.5522276126132602e-05, 
            -5.0542788517531417e-08])*(1E-12)**(np.array([2,3,4]))
smf3.gamma = 0.00045
smf3.core_d = 5.5E-6


#gain fiber, nufern ysf-HI
ydf1 = upp.FiberGain(0.6, grid_type='rel',z_grid=100)
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

#Pump parameters
p1P = 0.1    #pump power, CW


'''
#Plotting
tau = pulse.time
omega = pulse.freq

#create plot figure
fieldPlot, (t_ax, f_ax) = plt.subplots(2)   #set up plot figure
fieldPlot.suptitle('Pulse propagation profile')

#plot input
t_input, = t_ax.plot(tau,np.abs(pulse.At)**2, 'b--')    #plot time profile
t_ax.set_xlabel('Time (s)')
f_ax.plot(np.fft.fftshift(pulse.freq)/(2*np.pi),np.fft.fftshift(np.abs(pulse.getAf())**2), 'b--')  #plot freq profile
f_ax.set_xlabel('Frequency shift (Hz)')


t01, sig1 = upp.rmswidth(pulse.time,np.abs(pulse.At)**2)
output = upp.gratingPair(pulse, 1.0, 1500, 45)
t02, sig2 = upp.rmswidth(pulse.time,np.abs(output)**2)

input = pulse.At
pulse.At = output

#plot output
t_output, = t_ax.plot(tau,np.abs(pulse.At)**2, 'b-')    #plot time profile
t_ax.set_xlabel('Time (s)')
f_ax.plot(np.fft.fftshift(pulse.freq)/(2*np.pi),np.fft.fftshift(np.abs(pulse.getAf())**2), 'b-')  #plot freq profile
f_ax.set_xlabel('Frequency shift (Hz)')

plt.figlegend((t_input,t_output), ('Input', 'Output'), 'center right')
'''