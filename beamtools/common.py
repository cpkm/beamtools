'''
Common functions and constants for beamtools package

Created Thu Mar 23

@author: cpkmanchee
'''

h = 6.62606957E-34  #J*s
c = 299792458.0     #m/s

def normalize(f):
    return (f-f.min())/(f.max()-f.min())


def gaussian(x,sigma,amp=1,x0=0,const=0):
    return amp*np.exp(-(x-x0)**2/(2*sigma**2)) + const