from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
Frequency domain filters!

The functions in this particular sub-module return real or complex
vectors suitable for direct multiplication by an FFT.
"""

import numpy


def genFreqs( ndata, tSample ):
    """
    @brief Generate frequency vector of specified length and corresponding
    to a specified total time period.
    
    @param ndata    Number of elements in vector
    @param tSample  Sample time
    """
    dn = 2      # if you like the central frequency to be negative, change dn to 1
    return 1./(ndata*tSample) * numpy.hstack((numpy.arange(0., (ndata+dn)//2),
                                              numpy.arange(-(ndata+dn)//2+dn, 0)))

def genFreqsTOD(tod):
    n = tod.ctime.shape[0]
    return genFreqs(n, (tod.ctime[-1]-tod.ctime[0]) / (n-1))

def powerLawFilter( tod, power = -2.0, knee = 1):
    """
    """
    #f = genFreqs( tod.ndata, tod.sampleTime )
    f = genFreqsTOD(tod)
    filt = numpy.power(numpy.abs(f), power)
    filt[0] = 0.0
    return filt

def lowFreqWeinerFilter( tod, power = -2.0, fknee = 1.0 ):
    #f = genFreqs( tod.ndata, tod.sampleTime )
    f = genFreqsTOD(tod)
    s = numpy.power(numpy.abs(f)/fknee, power)
    s[0] = s[1]
    n = numpy.ones(len(f))
    filt = s/(s+n)
    return filt

def RCFilter( tod, fc = 2. ):
    #f = genFreqs( tod.ndata, tod.sampleTime )
    f = genFreqsTOD(tod)
    filt = 1/numpy.sqrt(1+numpy.power(f/fc, 2))
    return filt

def sine2highPass( tod=None, fc = 1.0, df = 0.1,
                   nsamps=None, sampleTime=None):
    """
    @brief Sine square high pass filter
    """
    fc = numpy.abs(fc); df = numpy.abs(df)
    if tod is None:
        f = genFreqsTOD(nsamps, sampleTime)
    else:
        f = genFreqsTOD(tod)
    filt = numpy.zeros(len(f))
    filt[numpy.abs(f) > fc + df/2.] = 1.0
    sel = (numpy.abs(f) > fc - df/2.)*(numpy.abs(f) < fc + df/2.)
    filt[sel] = numpy.sin(numpy.pi/2./df*(numpy.abs(f[sel]) - fc + df/2.))**2
    return filt

def sine2lowPass( tod=None, fc = 1.0, df = 0.1,
                  nsamps=None, sampleTime=None ):
    """
    @brief Sine square low pass filter
    @param fc frequency where power is half (Hz)
    @param df width of filter (Hz)
    """
    fc = numpy.abs(fc); df = numpy.abs(df)
    if tod is None:
        f = genFreqs(nsamps, sampleTime)
    else:
        #f = genFreqs( tod.ndata, tod.sampleTime )
        f = genFreqsTOD(tod)
    filt = numpy.zeros(len(f))
    filt[numpy.abs(f) < fc - df/2.] = 1.0
    sel = (numpy.abs(f) > fc - df/2.)*(numpy.abs(f) < fc + df/2.)
    filt[sel] = numpy.sin(numpy.pi/2*(1 - 1/df*(numpy.abs(f[sel]) - fc + df/2.)))**2
    return filt

def highPassButterworth( tod, fc = 1.0, order = 1, gain = 1.0 ):
    """
    @brief Generate a Butterworth high pass filter with specified cutoff frequency, order and gain

    @param fc cutoff frequency
    @param order order of the filter
    @param gain gain of the filter (default = 1.0)
    """
    j = numpy.complex(0,1)
    #f = j * genFreqs( tod.ndata, tod.sampleTime ) / fc
    f = j * genFreqsTOD(tod) / fc
    filt = numpy.ones(len(f))
    for k in range(1,order+1):
        sk = numpy.exp(j * (2*k + order - 1) * numpy.pi / (2*order))
        filt = filt * f / (1 - f*sk)
    filt *= gain

    return numpy.abs(filt)



def lowPassButterworth( tod, fc = 1.0, order = 1, gain = 1.0 ):
    """
    @brief Generate a Butterworth high pass filter with specified cutoff frequency, order and gain

    @param fc cutoff frequency
    @param order order of the filter
    @param gain gain of the filter (default = 1.0)
    """
    j = numpy.complex(0,1)
    #f = j * genFreqs( tod.ndata, tod.sampleTime ) / fc
    f = j * genFreqsTOD(tod) / fc
    filt = numpy.ones(len(f))
    for k in range(1,order+1):
        sk = numpy.exp(j * (2*k + order - 1) * numpy.pi / (2*order))
        filt = filt / (f - sk)
    filt *= gain

    return numpy.abs(filt)


def gaussianFilter( tod, timeSigma = None, frecSigma = None, gain = 1.0, fCenter = 0. ):
    """
    @brief  Generates a Gaussian filter.
    """
    if timeSigma is not None and frecSigma is not None:
        print("WARNING: cannot specify both time and frec sigmas. Using time_sigma.")
    if timeSigma is not None:
        sigma = 1.0 / (2*numpy.pi*timeSigma)
    elif frecSigma is not None:
        sigma = frecSigma
    else:
        sigma = 1.0

    #f = genFreqs( tod.ndata, tod.sampleTime )
    f = genFreqsTOD(tod)
    
    return gain * numpy.exp(-0.5*(numpy.abs(f)-fCenter)**2/sigma**2)

