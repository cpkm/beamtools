# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 20:28:42 2016

@author: cpkmanchee

Dispersion calculation from Sellmeier eqn
"""

import numpy as np
import sympy as sym

from beamtools.constants import h,c,pi

__all__ = ['sellmeier','ior']


def sellmeier(b_coefs, c_coefs, orders, l0=1.03E-6):
    '''Calculate dispersion parameters from Sellmeier coefficients.

    Inputs:
    b_coefs, B coefficients (unitless), should be 1x3 array (can be arb though)
    c_coefs, C coefficients (um^2), same shape as b_coefs
    orders, # of dispersion orders to calculate
    l0, central wavelength, m

    Output:
    beta, dispersion parameters, 1xm array, m=orders
    '''
    B = b_coefs
    C = c_coefs*1E-12
    w0 = 2*pi*c/l0

    if np.isscalar(l0):
        w0=np.array([w0])

    beta = np.zeros((orders+1,len(w0)))

    n, w = sym.symbols('n, w')  
    n = (1 + (B/(1-C*(w/(2*pi*c))**2)).sum() )**(1/2)

    for i in range(orders+1):
        #beta[i] = (1/c)*(i*sym.diff(n,w,i-1).subs(w,w0) + w0*sym.diff(n,w,i).subs(w,w0))
        bt = (1/c)*(i*sym.diff(n,w,i-1) + w*sym.diff(n,w,i))
        beta[i] = eval(str(bt),{'w':w0})

    return beta


def ior(b_coefs, c_coefs, l0=1.03E-6):
    '''Calculate index of refraction at given wavelength for Sellmeier coefs.
    '''

    B = b_coefs
    C = c_coefs*1E-12
    w0 = 2*pi*c/l0

    if np.isscalar(l0):
        return (1 + (B/(1-C*(w0/(2*pi*c))**2)).sum() )**(1/2)
 
    else:
        return np.array([(1 + (B/(1-C*(W/(2*pi*c))**2)).sum() )**(1/2) for W in w0])
