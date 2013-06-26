from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
from numpy import *

class filter:
    """
    Base class for filter type objects.  Such objects know how
    to construct filters, but do not pre-compute specific data for
    applying the filters to signals.
    """
    def __init__(self, domain='time', isLinear=True):
        self.domain = domain
        self.isLinear = isLinear
        
    def isLinear(self):
        return self.isLinear
    
    def naturalDomain(self):
        return self.domain

    @staticmethod
    def _timeVector(dt, n, positive=True):
        if positive:
            return float(dt)*arange(n)
        else:
            return float(dt)*arange(-n//2,n-n//2)

    @staticmethod
    def _freqVector(df, n, positive=False):
        if positive:
            idx = arange(n)
        else:
            idx = arange(n) - n*(arange(n) // ((n+1)//2))
        return float(df) * idx

    def getExec(self, *args, **kwargs):
        pass


class timeFilter(filter):
    """
    Base class for linear filters that operate naturally on vectors in
    the time domain.

    Such filters should only require the vector length and time-step
    to precompute necessary internal data.
    """
    def __init__(self, execClass=None):
        self.domain = 'time'
        self.isLinear = True
        self.eclass = execClass

    def getExec(self, dt, n):
        return self.eclass(self, dt=dt, n=n)

    def getExecFor(self, tod=None, data=None, t=None, dt=None, df=None, f=None, n=None):
        _n, _dt = None, None
        # Determine vector length, somehow
        for e in [ 'tod.data.shape[-1]', 'data.shape[-1]', 'len(data)', 'n']:
            try:
                __n = eval(e)
                if __n is not None: _n = __n
            except:
                continue
        if _n is None:
            raise RuntimeError("Could not determine vector length from arguments.")
        n = _n
        # Determine time step
        for e in ['tod.dt[1]', 'tod.ctime[1] - tod.ctime[0]', 't/n', '1/df/n',
                  '1/f', 'dt', 'tod.sampleTime']:
            try:
                __dt = eval(e)
                if __dt is not None: _dt = __dt
            except:
                continue
        if _dt is None:
            raise RuntimeError("Could not determine frequency step from arguments.")
        dt = _dt
        return self.getExec(dt, n)


class fourierFilter(filter):
    """
    Base class for filters that are (possibly complex) scalings in
    Fourier space).  Provides _freqVector and getExecFor methods;
    inheritors should extend getFourierFilter according to generate
    their particular filter forms.
    """
    def __init__(self):
        self.domain = 'freq'
        self.isLinear = True

    def getExec(self, df, n):
        data = self.getFourierFilter(df, n)
        return fourierFilterExec(self, df, n, data)

    def getExecFor(self, tod=None, data=None, dt=None, df=None, f=None, n=None):
        _n, _df = None, None
        # Determine vector length, somehow
        for e in [ 'tod.ndata', 'data.shape[-1]', 'len(data)', 'n']:
            try:
                __n = eval(e)
                if __n is not None: _n = __n
            except:
                continue
        if _n is None:
            raise RuntimeError("Could not determine vector length from arguments.")
        n = _n
        # Determine frequency step
        for e in ['1./(tod.ctime[1] - tod.ctime[0])', '1./tod.dt[1]/n', '1/dt/n', 'f/n', 'df']:
            try:
                __df = eval(e)
                if __df is not None: _df = __df
            except:
                continue
        if _df is None:
            raise RuntimeError("Could not determine frequency step from arguments.")
        df = _df
        return self.getExec(df, n)

    # Virtual methods
    def getFourierFilter(df, n):
        """
        Return a discrete representation of the filter.

        Parameters:
        df  frequency spacing.
        n   number of points.

        Returns:
        The fourier space scaling vector.
        """
        pass


"""
A filterExec is created from a filter to operate on a particular
stream of data.  It may precompute and store any information necessary
to filter the particular data that is targeted.

In some cases the filterExec needs to be designed specifically to
match the parent filter class.  In other cases (e.g. standard
Fourier-space filters) a single filterExec class can support a large
set of filters.
"""

class filterExec:
    def __init__(self, fclass):
        self.fclass = fclass

    def naturalDomain(self):
        return self.fclass.naturalDomain()

    def isLinear(self):
        return self.fclass.isLinear()
        
    def _shift(self, z, shift=True, inverse=False):
        if not shift:
            return z
        if inverse:
            pivot = self.n//2
        else:
            pivot = (self.n+1)//2
        return hstack((z[pivot:], z[:pivot]))
            
    def getTimes(self, shift=False, positive=True):
        t = self.fclass._timeVector(self.dt, self.n, positive=positive)
        return self._shift(t, shift=shift, inverse=True)
    
    def getFreqs(self, shift=False, positive=False):
        f = self.fclass._freqVector(self.df, self.n, positive=positive)
        return self._shift(f, shift=shift)

    def transform(self, data, inverse=False):
        if inverse:
            return fft.ifft(data)
        else:
            return fft.fft(data)

    def _applyDomain(self, data, domain, inverse=False, allowTransform=True):
        fdomain = self.fclass.naturalDomain()
        do_tf = (domain != fdomain)
        if do_tf and not allowTransform:
            raise RuntimeError('filterExec has natural domain "%s" but data has domain "%s"' % \
                  (fdomain, domain))
        reals_in = not iscomplexobj(data)
        if do_tf:
            inverse_tf = (domain=='freq')
            data = self.transform(data, inverse=inverse_tf)
        data = self.apply(data, inverse=inverse)
        if do_tf:
            data = self.transform(data, inverse=not inverse_tf)
        if reals_in: # then reals_out
            data = real(data)
        return data

    def applyTime(self, data, inverse=False, allowTransform=True):
        """
        Apply the filter to the (time-space) vector data.
        """
        return self._applyDomain(data, 'time', inverse=inverse, allowTransform=True)

    def applyFreq(self, data, inverse=False):
        """
        Apply the filter to the (frequency-space) vector data.
        """
        return self._applyDomain(data, 'freq', inverse=inverse, allowTransform=True)

    def getImpulseResponse(self, inverse=False, positive=False):
        # Just filter a delta function, remove dt scaling
        delta = (self.getTimes(positive=positive)==0.) \
                .astype('int').astype('float') / self.dt
        return self.applyTime(delta, inverse=inverse)

    def getFourierForm(self, inverse=False, shift=False):
        return self._shift(self.applyFreq(ones(self.n, dtype='complex')), shift)

    # Virtual methods
    def apply(self, data):
        return data


class timeFilterExec(filterExec):
    def __init__(self, fclass, dt, n, data=None):
        filterExec.__init__(self, fclass)
        self.dt, self.n = dt, n
        self.df = 1./(n*dt)
        self.data = data

    
class fourierFilterExec(filterExec):
    def __init__(self, fclass, df, n, data=None):
        filterExec.__init__(self, fclass)
        self.df, self.n = df, n
        self.dt = 1./(n*df)
        self.data = data
        
    def apply(self, data, inverse=False):
        if inverse:
            return data / self.data
        return data * self.data

    
"""
The usable filters.  Get the filter object with required parameters,
then call getExec() to get a filterExec that can be applied to your data.
"""

"""
Frequency-space linear filters.
"""

class highPassButterworth(fourierFilter):
    """
    High-pass Butterworth filter.
    """
    def __init__(self, fc=1., order=1, gain=1.):
        """
        @param fc cutoff frequency
        @param order order of the filter
        @param gain gain of the filter (default = 1.0)
        """
        fourierFilter.__init__(self)
        self.params = (fc, order, gain)

    def getFourierFilter(self, df, n):
        fc, order, gain = self.params
        f = 1j * self._freqVector(df, n) / fc
        filt = ones(len(f))
        for k in range(1,order+1):
            sk = exp(1j * (2*k + order - 1) * pi / (2*order))
            filt = filt * f / (1 - f*sk)
        filt *= gain
        return abs(filt)
        
class lowPassButterworth(fourierFilter):
    """
    Low-pass Butterworth filter.
    """
    def __init__(self, fc=1., order=1, gain=1.):
        """
        @param fc cutoff frequency
        @param order order of the filter
        @param gain gain of the filter (default = 1.0)
        """
        fourierFilter.__init__(self)
        self.params = (fc, order, gain)

    def getFourierFilter(self, df, n):
        fc, order, gain = self.params
        f = 1j * self._freqVector(df, n) / fc
        filt = ones(len(f))
        for k in range(1,order+1):
            sk = exp(1j * (2*k + order - 1) * pi / (2*order))
            filt = filt / (f - sk)
        filt *= gain
        return abs(filt)


class lowPassSine(fourierFilter):
    """
    Low-pass sine (or sine-squared, if you like) fourier-space filter.
    """
    def __init__(self, fc=1., df=1, gain=1.):
        """
        @param fc cut-off frequency (where gain is 1/2)
        @param df frequency half-width of transition region
        @param gain gain of the filter (default = 1.0)
        """
        fourierFilter.__init__(self)
        self.params = (fc, df, gain)

    def getFourierFilter(self, df, n):
        fc, dfc, gain = self.params
        f = (abs(self._freqVector(df, n)) - fc) / dfc  # -1 to 1 in transition
        filt = zeros(len(f))
        filt[f <= -1.] = 1.
        sel = (f > -1)*(f < 1.)
        filt[sel] = 0.5*(1.-sin(pi/2*f[sel]))
        return gain * filt

class highPassSine(fourierFilter):
    """
    High-pass sine (or sine-squared, if you like) fourier-space filter.
    """
    def __init__(self, fc=1., df=1, gain=1.):
        """
        @param fc cut-off frequency (where gain is 1/2)
        @param df frequency half-width of transition region
        @param gain gain of the filter (default = 1.0)
        """
        fourierFilter.__init__(self)
        self.params = (fc, df, gain)

    def getFourierFilter(self, df, n):
        fc, dfc, gain = self.params
        f = (abs(self._freqVector(df, n)) - fc) / dfc  # -1 to 1 in transition
        filt = zeros(len(f))
        filt[f >= 1.] = 1.
        sel = (f > -1)*(f < 1.)
        filt[sel] = 0.5*(1.+sin(pi/2*f[sel]))
        return gain * filt


class gaussianBand(fourierFilter):
    """
    Gaussian band-pass filter.
    """
    def __init__(self, fc=1., df=1, fwhm=None, gain=1.):
        """
        Call with parameters 'fc' and one of 'df' or 'fwhm'.

        @param fc central frequency
        @param df bandwidth (as standard deviation)
        @param fwhm bandwidth (as FWHM)
        """
        fourierFilter.__init__(self)
        if fwhm is not None:
            df = fwhm / (8.*log(2))**0.5
        self.params = (fc, df, gain)

    def getFourierFilter(self, df, n, ):
        fc, dfc, gain = self.params
        f = (abs(self._freqVector(df, n)) - fc) / dfc
        return gain*exp(-f**2 / 2)


class fourierFilterSet(fourierFilter):
    """
    Compose a set of fourierFilters.
    """
    def __init__(self, *args):
        """
        @params fourierFilter objects.
        """
        fourierFilter.__init__(self)
        self.filters = args

    def getFourierFilter(self, df, n):
        filt = self.filters[0].getFourierFilter(df, n)
        for f in self.filters[1:]:
            filt *= f.getFourierFilter(df, n)
        return filt


"""
Time-space linear filters.
"""

class timeConvolutionFilter(timeFilter):
    def __init__(self):
        #Set timeConvolutionExec as exec constructor
        timeFilter.__init__(self, timeConvolutionExec)

    # Virtual methods
    def getKernel(self, dt, n):
        return ones(1)

class timeConvolutionExec(timeFilterExec):
    def __init__(self, fclass, dt, n, mode='symmetric'):
        timeFilterExec.__init__(self, fclass, dt, n)
        self.kernel = self.fclass.getKernel(dt, n)
        self.mode = mode

    def apply(self, data, inverse=False):
        if inverse:
            raise RuntimeError('Convolution doesn\'t have simple inverse; use fourier filter.')
        # Keep the whole convolution
        z = convolve(data, self.kernel, mode='full')
        # Assume periodicity and combine wings
        n = self.kernel.shape[-1] - 1
        z[:n] += z[-n:]
        # Mode determines causality
        if self.mode == 'symmetric':
            return z[n//2:-n+n//2]
        elif self.mode == 'lag':
            return z[n:]
        elif self.mode == 'lead':
            return z[:-n]
        else:
            raise ValueError('Unknown phase mode "%s"' % self.mode)


class boxcar(timeConvolutionFilter):
    def __init__(self, width=None, width_samples=None):
        timeConvolutionFilter.__init__(self)
        self.width, self.width_samples = width, width_samples
        
    def getKernel(self, dt, n):
        w, w0 = self.width, self.width_samples
        if w0 is None:
            w0 = int(round(w / dt))
        return ones(w0, dtype='float') / w0

if __name__ == '__main__':
    import biggles as bg
    # testing, testing...
    filters = [
        ('hi-pass butter', highPassButterworth(fc=50., order=4)),
        ('lo-pass butter', lowPassButterworth(fc=50., order=4)),
        ('hi-pass sine', highPassSine(fc = 50., df=50.)),
        ('lo-pass sine', lowPassSine(fc = 50., df=50.)),
        ('gaussian band', gaussianBand(fc = 100., df=20.)),
        ('box-car', boxcar(width=0.02)),
        ]
    n, f_samp = 1000, 400.
    # Paneling
    n_cols = 3
    n_rows = (len(filters)+n_cols-1)//n_cols
    # Spectra
    table = bg.Table(n_rows, n_cols)
    for i, (name, filt) in enumerate(filters):
        fexec = filt.getExecFor(n=n, f=f_samp)
        x, y = fexec.getFreqs(shift=True), fexec.getFourierForm(shift=True)
        fg = bg.FramedPlot()
        fg.add(bg.Curve(x, abs(y)))
        fg.title = name
        table[i%n_rows, i//n_rows] = fg
    table.show()
    # Green functions -- shorter time interval, same density
    n, dt = n, 1e-4
    table = bg.Table(n_rows, n_cols)
    for i, (name, filt) in enumerate(filters):
        fexec = filt.getExecFor(n=n, dt=dt)
        x, y = fexec.getTimes(positive=False), fexec.getImpulseResponse(positive=False)
        fg = bg.FramedPlot()
        fg.add(bg.Curve(x, real(y)))
        fg.title = name
        table[i%n_rows, i//n_rows] = fg
    table.show()
    
    
