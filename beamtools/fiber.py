# -*- coding: utf-8 -*-
"""
Created on Thu Mar 08 20:47:42 2018

@author: cpkmanchee

Fiber optic calculations
"""

import numpy as np
import sympy as sym

from beamtools.constants import h,c,pi

__all__ = ['numapp','vnum','mfr','mfa']

def f2f():
    return 0

def numapp(n_core, n_clad, n_0=1):
    '''Calculate numerical aperature of optical fiber
    n_core = refractive index of core
    n_clad = refractive index of clad
    n_0 = surrounding index, usually air = 1
    '''
    return (1/n_0)*np.sqrt(n_core**2-c_clad**2)


def vnum(a_core,lambda_0, 
        na=None, n_core=None, n_clad=None, n_0=1):
    '''Calculate V-number of fiber
    a_core

    '''
    if all([i is not None for i in [n_core,n_clad]]):
        na = numapp(n_core,n_clad,n_0)

    return 2*pi*a_core*na/lambda_0


def mfr(a_core,v):
    '''Calculate mode radius
    marcuse = Marcuse calculation
    pete_corr = Petermann II correction
    '''
    marcuse = 0.65 + 1.619*v**(-3/2) + 2.879*v**(-6)
    pete_cor = 0.016 + 1.561*v**(-7)

    return a*(marcuse-pett_cor)

def mfa(w=None, a_core=None, v=None):
    '''Calculate mode field area
    '''
    if w is not None:
        return pi*w**2

    elif all([i is not None for i in [a_core,v]]):
        return pi*mfr(a_core,v)**2

    else:
        raise ValueError('Must provide either mode radius (w) '+
                'or V-number and core radius.')
