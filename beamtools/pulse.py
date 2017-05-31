'''
Pulse characterization

Created Fri May 12

@author: cpkmanchee
'''

import numpy as np
import os.path

from beamtools.constants import h,c,pi
from beamtools.common import normalize, gaussian, sech2
from beamtools.import_data_file import import_data_file as _import


__all__ = ['spectrumFT']


def spectrumFT(data,from_file = False, file_type='oo_spec', units_wl='nm', n_interp=0):
    '''Compute transform limited pulse from spectrum.
    data should be wavelength vs. PSD (intensity)
    Units assumed to be nm for wavelength.
    If from_file is set True, data should be filename
    Optional file_format, default is oceanoptics_spectrometer. Currently
    can not change this (filetype handling for x/y).
    n_interp = bit depth of frequency interpolation, n = 2**n_interp. 0 = auto
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
            print(wavelength,intensity)
            #get units from dataobject
        else:
            print('invalid filetype')
            return -1

    else:
        wavelength = data[0]
        intensity = data[1]

    if n_interp == 0:
        #insert here later - round up to nearest power of two.
        n = 2**12
    else:
        n = 2**12

    #use units to conver wavelength to SI
    wl = wavelength*1E-9
    psd = normalize(intensity)
    nu = c/wl       #nu is SI

    #interpolate psd, linear freq spacing
    nui = np.linspace(min(nu),max(nu),n)
    df = (max(nu)-min(nu))/(n-1)
    psdi = normalize(np.interp(nui,np.flipud(nu),np.flipud(psd)))
    #i = (np.abs(nui-nu0)).argmin()     #centre freq index

    #perform FT-1, remove centre spike
    t = np.fft.ifftshift(np.fft.fftfreq(n,df)[1:-1])
    ac =np.fft.ifftshift((np.fft.ifft(np.fft.ifftshift(psdi)))[1:-1])

    return wl, psd, nui, psdi, t, ac
