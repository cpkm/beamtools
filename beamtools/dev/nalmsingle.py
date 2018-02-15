import numpy as np
import matplotlib.pyplot as plt

from beamtools import upp, h, c

from tqdm import tqdm, trange

import sys
import shutil
import glob
import os
from datetime import datetime


def single_cavity(pulse, auto_z_step=False, save_folder=None):
    '''Define cavity round trip
    NOTE: pulse object is modified!!!
    returns:
        pulse.At = current pulse profile
        output_At = cavity output (outcoupled) profile
    '''
    if save_folder is not None:
        save = True
    else:
        save = False

    #Main section
    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='A')

    pulse.At = upp.propagate_fiber(pulse,smf1,autodz=auto_z_step)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='B')

    Ps = np.sum(np.abs(pulse.At)**2)*pulse.dt/tau_rt
    ydf1.gain = upp.calc_gain(ydf1,p1P,Ps)
    pulse.At = upp.propagate_fiber(pulse,ydf1,autodz=False)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='C')

    pulse.At = upp.propagate_fiber(pulse,smf1,autodz=auto_z_step)
    pulse.At,_ = upp.power_tap(pulse, 0, loss=0.06)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='D')

    pulse.At = upp.propagate_fiber(pulse,smf2,autodz=auto_z_step)
    pulse.At,_ = upp.power_tap(pulse, 0, loss=0.06)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='E')

    pulse.At = upp.propagate_fiber(pulse,smf2,autodz=auto_z_step)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='F')
    

    #NALM section
    nalmp = pulse.copyPulse()
    pulse.At, nalmp.At = upp.coupler_2x2(pulse, None, tap=50, loss=0.06)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='G')
        savepulse(nalmp,basename=save_folder,name='nalm',loc='M')


    Psn = np.sum(np.abs(pulse.At)**2)*pulse.dt/tau_rt
    ydf2.gain = upp.calc_gain(ydf2,p2P,Psn)

    #main pulse
    pulse.At = upp.propagate_fiber(pulse,smf1,autodz=auto_z_step)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='H')

    pulse.At = upp.propagate_fiber(pulse,ydf2,autodz=False)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='I')

    pulse.At = upp.propagate_fiber(pulse,smf1,autodz=auto_z_step)
    pulse.At,_ = upp.power_tap(pulse, 0, loss=0.06)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='J')

    pulse.At = upp.propagate_fiber(pulse,smf1,autodz=auto_z_step)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='K')

    pulse.At = upp.propagate_fiber(pulse,smf3,autodz=auto_z_step)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='L')

    pulse.At = upp.propagate_fiber(pulse,smf1,autodz=auto_z_step)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='M')

    #nalmp
    nalmp.At = upp.propagate_fiber(nalmp,smf1,autodz=auto_z_step)
    if save:
        savepulse(nalmp,basename=save_folder,name='nalm',loc='L')

    nalmp.At = upp.propagate_fiber(nalmp,smf3,autodz=auto_z_step)

    if save:
        savepulse(nalmp,basename=save_folder,name='nalm',loc='K')

    nalmp.At = upp.propagate_fiber(nalmp,smf1,autodz=auto_z_step)

    if save:
        savepulse(nalmp,basename=save_folder,name='nalm',loc='J')
    
    nalmp.At,_ = upp.power_tap(nalmp, 0, loss=0.06)
    nalmp.At = upp.propagate_fiber(nalmp,smf1,autodz=auto_z_step)

    if save:
        savepulse(nalmp,basename=save_folder,name='nalm',loc='I')

    nalmp.At = upp.propagate_fiber(nalmp,ydf2,autodz=False)

    if save:
        savepulse(nalmp,basename=save_folder,name='nalm',loc='H')

    nalmp.At = upp.propagate_fiber(nalmp,smf1,autodz=auto_z_step)

    if save:
        savepulse(nalmp,basename=save_folder,name='nalm',loc='G')

    #Main section
    pulse.At,_  = upp.coupler_2x2(pulse, nalmp, tap=50, loss=0.06)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='N')

    pulse.At = upp.propagate_fiber(pulse,smf2,autodz=auto_z_step)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='O')

    pulse.At = upp.optical_filter(pulse, filter_type='bpf', bandwidth=2E-9, loss=0.06, order=2)
    
    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='P')

    pulse.At = upp.propagate_fiber(pulse,smf2,autodz=auto_z_step)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='Q')

    pulse.At, output_At = upp.coupler_2x2(pulse, None, tap=75, loss=0.06)

    if save:
        savepulse(pulse,basename=save_folder,name='main',loc='A')
        savepulse(pulse.copyPulse(output_At),basename=save_folder,name='out',loc='Z')

    return pulse.At, output_At


def savepulse(pulse, basename, name='pulse', loc=''):
    '''Save pulse object'''
    fileext = '.pkl'
    output_num = 0
    filename = basename + name + loc + str(output_num).zfill(3) + fileext
    while not not glob.glob(filename):
        output_num = output_num + 1
        filename = basename + name + loc + str(output_num).zfill(3) + fileext
    upp.save_obj(pulse,filename)


folder = ('/Users/cpkmanchee/Documents/Code/Code Output/'
            'beamtools/beamtools/dev/nalmoscillator_output/20180209nalmoscillator-02/')
cavity_file= '20180209-235005-02cavity005.pkl'
output_file = '20180209-235005-02output006.pkl'

save_folder = folder + 'single pass/'
if not os.path.exists(save_folder):
    os.makedirs(save_folder)

pulse = upp.load_obj(folder+cavity_file)


#Define fiber components
smf1 = upp.Fiber(0.72142857)
smf1.alpha = 0.000576
smf1.beta = np.array([
            0.0251222977, 
            4.5522276126132602e-05, 
            -5.0542788517531417e-08])*(1E-12)**(np.array([2,3,4]))
smf1.gamma = 0.00045
smf1.core_d = 5.5E-6

smf2 = smf1.copyFiber(length=1.44285714)
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
pP = 0.15    #pump power, CW
Rp = 0.5
p1P = Rp*pP
p2P = (1-Rp)*pP

#Cavity
tau_rt = 1/(13.1E6)


