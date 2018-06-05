'''
Common methods for beamtools package

Created Fri May 12

@author: cpkmanchee
'''

import numpy as np
from beamtools.file_formats import file_formats

__all__ = ['normalize','rmbg','gaussian','sech2','lorentzian',
            'gaussian2D','rk4','moments','d4sigma','roi','alias_dict']


class Func:
    def __init__(self, value=None, index=None):
        self.val = value
        self.ind = index

    def at(self,x):
        return np.interp(x, self.ind, self.val)

    def diff(self):
        self.gradient = np.gradient(self.val)/np.gradient(self.ind)
            
    def diff_at(self,x):
        return np.interp(x,self.ind,self.gradient)


class FitResult():
    def __init__(self, ffunc, ftype, popt, pcov=0, indep_var='time', bgform='constant'):
        self.ffunc = ffunc
        self.ftype = ftype
        self.popt = popt
        self.pcov = pcov
        self.iv = indep_var
        self.bgform = bgform

    def subs(self,x):
        return self.ffunc(x,*self.popt)

    def get_args(self):
        return inspect.getargspec(self.ffunc)


class DataObj(dict):
    def __init__(self,d):
        self.__dict__ = d
    def fields(self):
        return self.__dict__.keys()
    def properties(self):
        [print(k,v) for k,v in file_formats[self.filetype].items()]
        return


def normalize(f, offset=0, method='normal'):
    '''Normalize array of data. Optional offset.
    '''
    norm = (f-f.min())/(f.max()-f.min()) + offset

    if method.lower() in ['area']:
        norm = norm/np.sum(norm)

    return norm

def rmbg(data, fit=None, form='constant'):
    '''Removes background from data
    data = [x,y]
    if sending poly fit params: p[0]*x**(N-1) + ... + p[N-1]

    return --> y - background
    '''
    if fit is None:
        #estimate background from given form
        if form.lower() in alias_dict['constant']:
            p = min(y)
        
        elif form.lower() in alias_dict['linear']:
            p = np.linalg.solve([[1,x[0]],[1,x[-1]]], [y[0],y[-1]])
            p = np.flipud(p)

        elif form.lower() in alias_dict['quadratic']:
            index = np.argmin(y)
            if index == 0:
                x3 = 2*x[0]-x[-1]
                y3 = y[-1]    
            elif index == len(y)-1:
                x3 = 2*x[-1]-x[0]
                y3 = y[0]   
            else:
                x3 = x[index]
                y3 = y[index]
            
            a = [[1,x[0],x[0]**2],[1,x[-1],x[-1]**2],[1,x3,x3**2]]
            b = [y[0],y[-1],y3]
            p = np.linalg.solve(a,b)
            p = np.flipud(p)
            
        else:
            print('Unknown background form')
            p = np.zeros((3))

    elif isinstance(fit,FitResult):
        #get background from FitResult object
        if fit.bgform.lower() in alias_dict['constant']:
            p=1
        elif fit.bgform.lower() in alias_dict['linear']:
            p=2
        elif fit.bgform.lower() in alias_dict['quadratic']:
            p=3
        else:
            p=1

        bg = np.polyval(fit.popt[-p:], data[0])

    elif any([type(fit) is z for z in [list,np.array]]):
        #background polynomial parameters supplied
        bg = np.polyval(fit,data[0])

    else:
        #Unknown or error
        print('Unknown fit argument.')
        bg = 0

    return data[1]-bg



def gaussian(x,sigma,amp=1,x0=0,const=0,sg=1):
    '''Gaussian distribution.
    x = independent variable
    sigma = sd (width parameter)
    x0 = centre position
    amp = amplitude
    const = y-offset
    '''
    return amp*np.exp(-((x-x0)**2/(2*sigma**2))**sg) + const

def sech2(x,sigma,amp=1,x0=0,const=0):
    '''Hyperbolic secant-squared distribution.
    x = independent variable
    sigma = width parameter
    x0 = centre position
    amp = amplitude
    const = y-offset
    '''
    return amp*(1/np.cosh((x-x0)/sigma))**2 + const

def lorentzian(x,sigma,amp=1,x0=0,const=0):
    '''Lorentzian distribution.
    x = independent variable
    sigma = width parameter
    x0 = centre position
    amp = amplitude
    const = y-offset
    '''
    return amp*(sigma**2/((x-x0)**2+sigma**2)) + const


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

def rk4(f, x, y0, const_args=[], abs_x=False):
    '''
    functional form
    y'(x) = f(x,y,constants)

    f must be function, f(x,y,const_args)
    x = array
    y0 = initial condition,
    cont_args = additional constants required for f

    returns y, integrated array
    '''

    N = np.size(x)
    y = np.zeros(np.shape(x))
    y[0] = y0
    dx = np.gradient(x)

    if abs_x:
        dx = np.abs(dx)

    for i in range(N-1):
        k1 = f(x[i], y[i], *const_args)
        k2 = f(x[i] + dx[i]/2, y[i] + k1*dx[i]/2, *const_args)
        k3 = f(x[i] + dx[i]/2, y[i] + k2*dx[i]/2, *const_args)
        k4 = f(x[i] + dx[i], y[i] + k3*dx[i], *const_args)

        y[i+1] = y[i] + (k1 + 2*k2 + 2*k3 + k4)*dx[i]/6

    return y

def moments(x,y):
    '''calculate 1st and second moments of distribution
        returns first and second moments

    first moments are averages in each direction
    second moments are variences in x, y and diagonal
    '''
    dx = np.gradient(x)

    A = np.sum(y*dx)
    
    #first moments (averages)
    avgx = np.sum(y*x*dx)/A

    #calculate second moments if required

    sig2x = np.sum(y*(x-avgx)**2*dx)/A
        
    return np.array([avgx,sig2x])

def d4sigma(x,y):
    '''calculate d4sigma width
    '''
    [av,sig2] = moments(x,y)

    return 4*(sig2)**(1/2)

def roi(x,y,x_lim):
    '''return subset of x,y arrays
    '''
    ind = [np.argmin(np.abs(x-xl)) for xl in x_lim]
    x_out = x[ind[0]:ind[1]]
    y_out = y[ind[0]:ind[1]]

    return np.array([x_out,y_out])


alias_dict = {
        'gaus': ('gaus','gaussian','g'),
        'sech2': ('sech2','secant squared','hyperbolic secant squared','s'),
        'lorentz': ('lorentz','lorentzian','cauchy','cauchy-lorentz','l-c'),
        'GausFit': ('GausFit','gaus','gaussian'),
        'Sech2Fit': ('Sech2Fit', 'sech2','secant squared','hyperbolic secant squared'),
        'constant': ('constant','const','c'),
        'linear': ('linear','lin','l'),
        'quadratic': ('quadratic', 'quad', 'q'),
        'symmetric': ('symmetric','sym','s'),
        'asymmetric': ('asymmetric','asym','a')
        }