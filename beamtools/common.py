'''
Common methods for beamtools package

Created Fri May 12

@author: cpkmanchee
'''

import numpy as np

__all__ = ['normalize','gaussian','sech2', 'gaussian2D', 'alias_dict']

def normalize(f, offset=0):
    '''Normalize array of data. Optional offset.
    '''
    return (f-f.min())/(f.max()-f.min()) + offset

def gaussian(x,sigma,amp=1,x0=0,const=0):
    '''Gaussian distribution.
    x = independent variable
    sigma = sd (width parameter)
    x0 = centre position
    amp = amplitude
    const = y-offset
    '''
    return amp*np.exp(-(x-x0)**2/(2*sigma**2)) + const

def sech2(x,sigma,amp=1,x0=0,const=0):
    '''Hyperbolic secant-squared distribution.
    x = independent variable
    sigma = width parameter
    x0 = centre position
    amp = amplitude
    const = y-offset
    '''
    return amp*(1/np.cosh((x-x0)/sigma))**2 + const


def gaussian2D(xy_meshgrid,x0,y0,sigx,sigy,amp,const,theta=0):
    '''Generates a 2D gaussian surface of size (n x m).
    
    Inputs:
        xy_meshgrid = [x,y]
        x = meshgrid of x array
        y = meshgrid of y array
        
    where x and y are of size (n x m)
    n = y.shape[0] (or x.) = number of rows
    m = x.shape[1] (or y.) = number of columns
    
        x0,y0 = peak location
        sig_ = standard deviation in x and y, gaussian 1/e radius
        amp = amplitude
        const = offset (constant)
        theta = rotation parameter, 0 by default
    
    Output:
        g.ravel() = flattened array of gaussian amplitude data
    
    where g is the 2D array of gaussian amplitudes of size (n x m)
    '''
    x = xy_meshgrid[0]
    y = xy_meshgrid[1]

    a = np.cos(theta)**2/(2*sigx**2) + np.sin(theta)**2/(2*sigy**2)
    b = -np.sin(2*theta)/(4*sigx**2) + np.sin(2*theta)/(4*sigy**2)
    c = np.sin(theta)**2/(2*sigx**2) + np.cos(theta)**2/(2*sigy**2)

    g = amp*np.exp(-(a*(x-x0)**2 -b*(x-x0)*(y-y0) + c*(y-y0)**2)) + const
       
    return g.ravel()


alias_dict = {
        'gaus': ('gaus','gaussian','g'),
        'sech2': ('sech2','secant squared','hyperbolic secant squared','s'),
        'GausFit': ('GausFit','gaus','gaussian'),
        'Sech2Fit': ('Sech2Fit', 'sech2','secant squared','hyperbolic secant squared'),
        'constant': ('constant','const','c'),
        'linear': ('linear','lin','l'),
        'quadratic': ('quadratic', 'quad', 'q'),
        'symmetric': ('symmetric','sym','s'),
        'asymmetric': ('asymmetric','asym','a')
        }