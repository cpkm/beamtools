# -*- coding: utf-8 -*-
"""
Created on Thu Mar 08 20:47:42 2018

@author: cpkmanchee

Fiber optic calculations
"""

import numpy as np
import sympy as sym

from beamtools.constants import h,c,pi

__all__ = ['butt_eff','numapp','vnum','mfr','mfa']

def butt_eff(w1_0, w2_0, z, dx, theta, lambda0, 
        dy=0, phi=0):
    '''Butt coupling efficiency for two fibers
    w1 = mfr of fiber 1
    w2 = mfr of fiber 2
    z = separation along optical axis
    d = separation perp to optical axis
    theta = angle x in rad
    '''
    k = 2*pi/lambda0
    A = (1/w1**2+1/w2**2)**2 + (k**2/4)*(1/r1-1/r2)**2

    w1 = w1_0*(1+(2*z/(k*w1_0))**2)**(1/2)
    w2 = w2_0*(1+(2*z/(k*w2_0))**2)**(1/2)
    r1 = z*(1+(k*w1_0**2/(2*z))**2)
    r2 = z*(1+(k*w2_0**2/(2*z))**2)

    nz = (2/(w1*w2))*A**(-1/2)

    ndx = np.exp(-2*dx**2*(
            ((1/(w1**2*w2**2))*(1/w1**2+1/w2**2)
            + (k**2/4)*(1/(w1*r1)**2-1/(w2*r2)**2))/A))

    ndy = np.exp(-2*dy**2*(
            ((1/(w1**2*w2**2))*(1/w1**2+1/w2**2)
            + (k**2/4)*(1/(w1*r1)**2-1/(w2*r2)**2))/A))

    nth = np.exp((-k**2*theta**2/2)*(
            (1/w1**2+1/w2**2)/A))

    nph = np.exp((-k**2*phi**2/2)*(
            (1/w1**2+1/w2**2)/A))

    ndxth = np.exp((-k**2*theta*dx)*(
            (1/w1**2+1/w2**2)/A))

    ndyph = np.exp((-k**2*phi*dy)*(
            (1/w1**2+1/w2**2)/A))

    return nz*ndx*nth*ndxth*nz*ndy*nph*ndyph


def numapp(n_core, n_clad, n_0=1):
    '''Calculate numerical aperature of optical fiber
    n_core = refractive index of core
    n_clad = refractive index of clad
    n_0 = surrounding index, usually air = 1
    '''
    return (1/n_0)*np.sqrt(n_core**2-c_clad**2)


def vnum(a_core,lambda0, 
        na=None, n_core=None, n_clad=None, n_0=1):
    '''Calculate V-number of fiber
    a_core

    '''
    if all([i is not None for i in [n_core,n_clad]]):
        na = numapp(n_core,n_clad,n_0)

    return 2*pi*a_core*na/lambda0


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
