from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from past.builtins import basestring

import numpy
#import pyfftw

def power( x, dt = 1., nbin = 1, binsize = 0, detrend = False,
           Bartlett = False, Welch = False, Hann = False,
           useRegular = False):
    """
    Take the power spectrum of input x.
    
    dt is the time sample spacing (seconds).

    Can break a long input up into sections or "bins" before transforming.
    Frankly, bin seems like a crap name for this, but it's in the interface,
    so I guess we are stuck with it.  Read "section" wherever you see "bin".
    nbin is the number of sections to take.
    binsize is the number of samples per section.
    Caller can specify either nbin or binsize but not both.
    When you ask for 2 sections, actually average power over the 1st 50%
    of data, the middle 50%, and the last 50%; in general, averages power
    over (2nbin-1) half-offset ranges.
    
    Returns (power, nu, window) where power is the power spectrum in
    units of x^2 per Hz, nu is the frequency sample vector in Hz,
    and window is the window function weight applied to the timestream.
    """

    if binsize and nbin != 1:
        raise ValueError("Please specify either binsize or nbin, but not both")

    nsamp = numpy.shape(x)[-1]

    if binsize == 0:
        binsize = int(nsamp/nbin)             # length of a bin
    else:
        nbin = int(nsamp/binsize)

    if nbin <= 0:
        raise ValueError("You have requested %d bins (len = %d)" % \
            (nbin, nsamp))

    if detrend: detrendData(x, window = 200)

    if Bartlett + Hann + Welch > 1:
        raise ValueError("Please choose at most one type of window")

    if Bartlett:
        window = 1 - abs((numpy.arange(binsize) - binsize/2.0)/(binsize/2.0))
    elif Hann:
        window = 0.5*(1-numpy.cos(2*math.pi*(numpy.arange(binsize))/binsize))
    elif Welch:
        window = 1 - pow((numpy.arange(binsize) - binsize/2.0)/(binsize/2.0), 2)
    else:
        window = 1.0*numpy.ones(binsize)

    
    one_d = x.ndim == 1
    if one_d:
        x.shape = (1,-1)

    if useRegular:
        nt = nextregular(binsize)
    else:
        nt = binsize

    power = 0
    if nbin != 1:
        for b in range(2*nbin - 1):
            y = x[:,b*binsize//2 : b*binsize//2 + binsize].copy()
            detrendData(y, window = 200)
            fx = numpy.fft.rfft(window[numpy.newaxis,:]*y,nt)
            power += (fx.real*fx.real + fx.imag*fx.imag)
            #fx = numpy.fft.fft(window[numpy.newaxis,:]*y)
            #fx = pyfftw.interfaces.numpy_fft.fft(window*y)
            #power += (fx.real*fx.real + fx.imag*fx.imag)[:,:binsize//2+1]
    else:
        y = x.copy()
        detrendData(y, window = 200)
        fx = numpy.fft.rfft(window[numpy.newaxis,:]*y,nt)
        power = (fx.real*fx.real + fx.imag*fx.imag)
        #fx = numpy.fft.fft(window[numpy.newaxis,:]*x)
        #fx = pyfftw.interfaces.numpy_fft.fft(window*x)
        #power = (fx.real*fx.real + fx.imag*fx.imag)[:,:binsize//2+1]
            

    # Normalizations
    power_scale = 2.0                     # 1-sided power
    power_scale /= (2.*nbin - 1.)         # Allow for multiple 'bins'
    power_scale /= pow(sum(window),2)     # Allow for window
    power_scale *= float(nt)*dt      # per unit time (total time in a bin)
    #power_scale *= float(binsize)*dt      # per unit time (total time in a bin)
    power *= power_scale

    # 1-sided power double-counts at 0 and Nyquist frequency.  Correct that.
    #power[0] /= 2.
    #power[-1] /= 2.

    nf = power.shape[-1]
    nu = numpy.arange(nf) / nf / (2*dt)
    nu[0] = 0.5*nu[1]  # Make sure the x=0 isn't killing power

    if one_d:
        x.shape = -1
        power.shape = -1

    return power, nu, window

def detrendData(y, window = 1000):
    """
    @brief Remove the trend and mean from a data vector
    @param y        Data to detrend
    @param window   Number of elements to consider at each end of vector
    """
    n = y.shape[-1]
    one_d = y.ndim == 1
    if one_d: y.shape = (1,-1)
    if window > n//2: window = n//2
    y0 = numpy.mean(y[:,:window],axis=1)
    y1 = numpy.mean(y[:,-window:],axis=1)
    m1 = (y1+y0)/2.0
    m2 = numpy.mean(y,axis=1)
    slope = (y1-y0)/(n-1)
    x = numpy.arange(n)
    y -= (y0 - m1 + m2)[:,numpy.newaxis].repeat(n,1) + slope[:,numpy.newaxis] * x[numpy.newaxis,:]
    if one_d: y.shape = -1


def nextregular(n):
    while not checksize(n): n+=1
    return n

def checksize(n):
    while not (n%16): n//=16
    while not (n%13): n//=13
    while not (n%11): n//=11
    while not (n%9): n//=9
    while not (n%7): n//=7
    while not (n%5): n//=5
    while not (n%3): n//=3
    while not (n%2): n//=2
    return (1 if n == 1 else 0)


