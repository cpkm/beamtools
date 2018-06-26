'''
Pulse characterization

Created Fri May 12 2017

@author: cpkmanchee
'''

import numpy as np
import os.path

from inspect import getargspec
from warnings import warn

from beamtools.constants import h,c,pi
from beamtools.common import normalize, gaussian, sech2, alias_dict, FitResult
from beamtools.import_data_file import import_data_file as _import
from beamtools.import_data_file import DataObj
from scipy.optimize import curve_fit


__all__ = ['autocorr','spectrumFT', 'fit_ac', 'ac_x2t', 'sigma_fwhm']


def autocorr(x):
    '''Calculate autocorrelation of function
    '''
    return np.correlate(x, x, mode='same')


def spectrumFT(data,from_file = False, file_type='oo_spec', units_wl='nm', n_interp=0):
    '''Compute transform limited pulse from spectrum.
    data = wavelength vs. PSD (intensity) if from_file=False
        = filename of spectrum file to be imported if from_file=True
    Units assumed to be nm for wavelength.
    If from_file is set True, data should be filename
    Optional file_format, default is oceanoptics_spectrometer. Currently
    can not change this (filetype handling for x/y).
    n_interp = bit depth of frequency interpolation, n = 2**n_interp. 0 = auto
    
    Et/Ew are envelope fulctions of time/freq electric field
    
    w = 2pi*nu

    S(w) = |Ew|**2
    I(t) = |Et|**2
    '''

    if from_file:
        if type(data) is str:
            if not os.path.exists(data):
                print('File does not exist')
                return -1
        
            imported_data = _import(data,file_type)

            #insert testing for wavelength/intensity location in dataobject
            wavelength = imported_data.wavelength
            intensity = imported_data.intensity
            #get units from dataobject
        else:
            print('invalid filetype')
            return -1

    else:
        wavelength = data[0]
        intensity = data[1]
        imported_data = data

    if n_interp == 0:
        #insert here later - round up to nearest power of two.
        n = 2**12
    else:
        n = 2**12

    #use units to convert wavelength to SI
    wl = wavelength*1E-9
    psd = normalize(intensity)
    nu = c/wl       #nu is SI

    #interpolate psd, linear freq spacing
    nui = np.linspace(min(nu),max(nu),n)
    df = (max(nu)-min(nu))/(n-1)
    Ew = (normalize(np.interp(nui,np.flipud(nu),np.flipud(psd))))**(1/2)
    #i = (np.abs(nui-nu0)).argmin()     #centre freq index

    #perform FT-1, remove centre spike
    t = np.fft.fftshift(np.fft.fftfreq(n,df)[1:-1])
    Et = np.fft.fftshift((np.fft.ifft(np.fft.ifftshift(Ew)))[1:-1])

    output_dict = {'time': t, 'et': Et, 'nu': nui, 'ew': Ew}
    output = DataObj(output_dict)

    return output, imported_data


def ac_x2t(position,aoi=15,config='sym'):
    '''Convert autocorrelation position to time
    Symmetric - stage moves perp to normal.
    Asymmetric - stage moves along incoming optical axis
    '''
    if type(config) is not str:
        print('Unrecognized configuration. Must be symmetric or asymmetric.')
        return position

    if config.lower() in alias_dict['symmetric']:
        time = (1/c)*position*2*np.cos(aoi*pi/180)

    elif config.lower() in alias_dict['asymmetric']:
        time = (1/c)*2*position
        #time = (1/c)*position*(1+np.cos(2*aoi*pi/180))

    else:
        print('Unrecognized configuration. Must be symmetric or asymmetric.')
        #I think this is a general expression, given stage_angle wrt to mirror normal
        #time = (1/c)*2*position*np.cos((aoi-stage_angle)*pi/180)
        return position

    return time


def fit_ac(data, from_file = False, file_type='bt_ac', form='all', bgform = 'constant'):
    '''Fit autocorrelation peak.
    data must be either:
        1. 2 x n array - data[0] = time(delay), data[1] = intensity
        2. datafile name --> from_file must be True

    If there is no 'delay' parameter in data file (only position), the position is
     auto converted to time delay.
    '''
    if from_file:
        if type(data) is str:
            if not os.path.exists(data):
                warn('File does not exist')
                return -1
        
            imported_data = _import(data,file_type)

            #insert testing for power location in dataobject
            position = imported_data.position
            intensity = imported_data.power

            if 'delay' in imported_data.__dict__:
                delay = imported_data.delay
            else:
                delay = ac_x2t(position,aoi=15,config='sym')

            #get units from dataobject
        else:
            warn('invalid filetype')
            return -1

    else:
        imported_data = data
        delay = data[0]
        intensity = data[1]

    x = delay
    y = intensity

    bgpar, bgform = _background(x,y,form = bgform)
    mean = np.average(x,weights = y)       
    stdv = np.sqrt(np.average((x-mean)**2 ,weights = y))
      
#set fitting function (including background)
    if bgform is None:
        def fitfuncGaus(x,sigma,a,x0):
            return gaussian(x,sigma,a,x0)
        def fitfuncSech2(x,sigma,a,x0):
            return sech2(x,sigma,a,x0)

    if bgform.lower() in alias_dict['constant']:
        def fitfuncGaus(x,sigma,a,x0,p0):
            return gaussian(x,sigma,a,x0) + p0
        def fitfuncSech2(x,sigma,a,x0,p0):
            return sech2(x,sigma,a,x0) + p0
        
    elif bgform.lower() in alias_dict['linear']:
        def fitfuncGaus(x,sigma,a,x0,p1,p0):
            return gaussian(x,sigma,a,x0) + p1*x + p0
        def fitfuncSech2(x,sigma,a,x0,p1,p0):
            return sech2(x,sigma,a,x0) + p1*x + p0

    elif bgform.lower() in alias_dict['quadratic']:
        def fitfuncGaus(x,sigma,a,x0,p2,p1,p0):
            return gaussian(x,sigma,a,x0) + p2*x**2 + p1*x + p0
        def fitfuncSech2(x,sigma,a,x0,p2,p1,p0):
            return sech2(x,sigma,a,x0) + p2*x**2 + p1*x + p0
    else:
        def fitfuncGaus(x,sigma,a,x0):
            return gaussian(x,sigma,a,x0)
        def fitfuncSech2(x,sigma,a,x0):
            return sech2(x,sigma,a,x0)

    nFitArgs = len(getargspec(fitfuncGaus).args) - 1

#sets which functions are to be fit... this can be streamlined i think
    if form.lower() in ['both', 'all']:
        fitGaus = True
        fitSech2 = True
        
    elif form.lower() in alias_dict['gaus']:
        fitGaus = True
        fitSech2 = False
        
    elif form.lower() in alias_dict['sech2']:
        fitGaus = False
        fitSech2 = True
        
    else:
        print('Unknown fit form: '+form[0])
        fitGaus = False
        fitSech2 = False
    
    #start fitting 
    popt=[]
    pcov=[]
    fit_results=[]
    
    if type(bgpar) is np.float64:
        p0=[stdv,max(y)-min(y),mean,bgpar]
    elif type(bgpar) is np.ndarray:
        p0=[stdv,max(y)-min(y),mean]+bgpar.tolist()   
    else:
        p0=None
   
    if fitGaus:
        try:
            poptGaus,pcovGaus = curve_fit(fitfuncGaus,x,y,p0) 
        except RuntimeError:
            poptGaus = np.zeros(nFitArgs)
            pcovGaus = np.zeros((nFitArgs,nFitArgs))   
        
        popt.append(poptGaus)
        pcov.append(pcovGaus)
        fit_results.append(FitResult(ffunc=fitfuncGaus, ftype='gaussian',
            popt=poptGaus, pcov=pcovGaus, indep_var='time', bgform=bgform))
        
    if fitSech2:
        try:
            poptSech2,pcovSech2 = curve_fit(fitfuncSech2,x,y,p0)
        except RuntimeError:
            poptSech2 = np.zeros(nFitArgs)
            pcovSech2 = np.zeros((nFitArgs,nFitArgs))
                       
        popt.append(poptSech2)
        pcov.append(pcovSech2)
        fit_results.append(FitResult(ffunc=fitfuncSech2, ftype='sech2',
            popt=poptSech2, pcov=pcovSech2, indep_var='time', bgform=bgform))

    return fit_results, imported_data


def sigma_fwhm(sigma, shape='gaus'):
    '''Convert sigma to full-width half-max
    '''
    if shape.lower() in alias_dict['gaus']:
        A = 2*np.sqrt(2*np.log(2))
    elif shape.lower() in alias_dict['sech2']:
        A = 2*np.arccosh(np.sqrt(2))
    elif shape.lower() in alias_dict['lorentz']:
        A = 1
    else:
        A = 1
        warn('Pulse shape not recognized.')

    return A*sigma


def deconv(sigma, shape='gaus'):
    '''Deconvolution factors
    Go from sigma_ac --> sigma
    '''
    if shape.lower() in alias_dict['gaus']:
        A = 1/np.sqrt(2)
    elif shape.lower() in alias_dict['sech2']:
        A = 0.6482
    elif shape.lower() in alias_dict['lorentz']:
        A = 0.5
    else:
        A = 1
        warn('Pulse shape not recognized.')

    return A*sigma


def _background(x,y,form = 'constant'):
    '''Provides starting values for background parameters.
    Takes x,y data and the desired background form (default to constant)
    returns p, the polynomial coefficients. p is variable in length.
    '''
    if form is None:
        p = np.zeros((3))

    if form.lower() in alias_dict['constant']:
        p = min(y)
        #p = np.hstack((p,[0,0]))
        
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
        warn('Unknown background form')
        p = np.zeros((3))
        
    return p, form
