import numpy as np
import beamtools as bt
import matplotlib.pyplot as plt

def autocorr(x):
    return np.correlate(x, x, mode='same')

l0 = 3E-5
nu0 = bt.c/l0
w0 = 2*bt.pi*nu0

t0 = 0.1E-12
window = 50E-12
Nt = 2**16

t = np.linspace(-window/2,window/2,Nt)
dt = window/(Nt-1)

#Et = bt.gaussian(t,sigma=t0)**(1/2)*np.exp(1j*w0*t)
Et = bt.gaussian(t,sigma=t0)**(1/2)
It = np.abs(Et)**2

Ae = autocorr(Et)
Ai = autocorr(It)

Ew = np.fft.fft(Et)
psd = np.abs(np.fft.fft(Ae))
nu = np.fft.fftfreq(Nt,dt)
Nnu = nu.size
dnu = (nu.max()-nu.min())/(Nnu-1)
Pi = np.fft.fft(It)

EtF = np.fft.ifft(np.abs(Ew))

AeF = autocorr(np.fft.fftshift(EtF))

iFpsd = np.fft.ifft(psd)
it = np.fft.fftfreq(Nnu,dnu)

plt.plot(it,bt.normalize(np.abs(EtF)), '-',c='C0')
plt.plot(t,bt.normalize(np.abs(Et)),'--',c='C1')

plt.plot(t,bt.normalize(np.abs(AeF)),'--',c='C3')
plt.plot(t,bt.normalize(np.abs(Ae)),'-.',c='C2')

plt.figure()
plt.plot(nu,bt.normalize(np.abs(Ew)**2), '-',c='C0')
plt.plot(nu,bt.normalize(psd),'--',c='C1')


plt.show()

#plt.plot(t,np.real(Et),'-',c='C9')
#plt.plot(t,np.abs(Et),'-.',c='C9')
#plt.plot(t,np.real(Ae)/np.abs(Ae[Ae.size//2]),'-',c='C0')
#plt.plot(t,bt.normalize(np.abs(Ae)),'-.', c='C0')
#plt.plot(t,It,'-',c='C5')
#plt.plot(t,bt.normalize(Ai),'-.',c='C9')
#plt.plot(t,bt.normalize(np.abs(Ae)**2),'-', c='C0')
#plt.plot(it,bt.normalize(np.abs(iFpsd)**2),c='C4')



#plt.plot(t,bt.normalize(np.abs(Ae)**2),'-.',c='C9')

#plt.plot(t,It,'-r')
#plt.plot(t,bt.normalize(np.real(Ai)),'--', c='C3')

plt.show()