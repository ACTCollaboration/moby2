from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy
from numpy.fft import fft2,ifft2,fftfreq

import moby2.util.log as psLib

#
# If you add preconditioners, create a key for them in the
# preconditioner_registry at the bottom of this file.
#

class preconditionerList(list):
    """
    @brief a class for combining multiple preconditioners
    """
    def applyPrecond( self, arr ):
        """
        Apply, in order, the preconditioners
        """
        for p in self:
            p.applyPrecond( arr )

class preconditioner(object):
    """
    @brief preconditioner class for use with pcg.pcg object

    preconditioners act on numpy arrays which are double precision proxies for maps.
    """
    def __init__( self, map):
        """
        @param mp a fitsMap that this preconditioner will operate on
        """
        pass

    def applyPrecond(self, arr):
        """
        @brief This is the function used by pcg to apply the preconditioner
        @param arr numpy array to precondition
        """
        pass


class meanRemover( preconditioner ):
    """
    @brief Remove the mean of the map.
    """
    def __init__( self, mp ):
        """
        """
        self.mask = mp.weight != 0

    def applyPrecond( self, arr ):
        """
        """
        m = arr[self.mask].mean()
        arr[self.mask] -= m

        psLib.trace( "moby", 3, "MeanRemover: Mean = %f" % m )

class trendRemover( preconditioner ):
    """
    @brief Remove the trend of the map.
    """
    def __init__( self, mp ):
        """
        """
        self.mask = mp.weight != 0
        #self.mask = self.mask.ravel()
        if self.mask.sum() == 0:
            psLib.trace( "moby", 2, "trendRemover: no weight in map!" )
            return
        x, y = numpy.indices(mp.weight.shape)
        #x = x.ravel()
        #y = y.ravel()
        self.A = numpy.vstack((x[self.mask], y[self.mask], x[self.mask]*0+1))
        self.invAA = numpy.linalg.inv(numpy.dot(self.A, self.A.transpose()))

    def applyPrecond( self, arr ):
        """
        @brief Remove trend from map
        @param arr   map data from where to subtract the trend
        """
        if self.mask.sum() == 0:
            return
        mz = numpy.dot(self.A.transpose(),
                       numpy.dot(self.invAA, numpy.dot(self.A, arr[self.mask])))
        arr[self.mask] -= mz


# Not updated
#
class trendRemoverPlanet( preconditioner ):
    """
    @brief Remove the trend of the map, given a giant source in the map.
    """
    def __init__( self, mp, lon = 0., lat = 0., cutBoxWidth = -1, weightLowerLimit = 0, power = 1 ):
        """
        @param mp map
        @param lon longitude of source
        @param lat latitude of source
        @param cutBoxRadius radius of mask around source in arcminutes
        @param weightLowerLimit mask out data with weight lower than this for the background fit
        @param power the power of the polynomial fit to remove. ie: 0 = mean, 1 = linear, 2 = quadratic, 3 = cubic
        """
        self.lon = lon
        self.lat = lat
        self.cutBoxWidth = cutBoxWidth
        self.weightLowerLimit = weightLowerLimit
        self.power = power
        self.mask = numpy.ones(numpy.shape(mp.weight), dtype = 'bool')
        self.mask[numpy.where(mp.weight <= weightLowerLimit)] = False
        if cutBoxWidth > 0:
            r0, c0 = mp.sky2pix( lon, lat )
            cbrRow = cutBoxWidth/60./2 / mp.dLatPerPix()
            cbrCol = cutBoxWidth/60./2 / mp.dLonPerPix()
            self.mask[r0-cbrRow:r0+cbrRow, c0-cbrCol:c0+cbrCol] = False
        self.mask = self.mask.ravel()
        x, y = numpy.indices(numpy.shape(mp.weight))
        x = x.ravel() - mp.ncol/2 #//?
        y = y.ravel() - mp.nrow/2 #//?
        #x = x[self.mask]
        #y = y[self.mask]
        x_k = numpy.zeros( (power, len(x)) )
        y_k = numpy.zeros( (power, len(x)) )
        k = 1
        while k <= power:
            x_k[(k-1)] = x**k
            y_k[(k-1)] = y**k
            k += 1
                                            
        self.A = numpy.vstack((x_k[:,self.mask], y_k[:,self.mask], numpy.ones(len(x[self.mask]))))
        #self.A = numpy.vstack((x[self.mask], y[self.mask], numpy.ones(len(x[self.mask]))))
        self.invAA = numpy.linalg.inv(numpy.dot(self.A, self.A.transpose()))

        self._mask = numpy.ones(numpy.shape(mp.weight), dtype = 'bool')
        self._mask[numpy.where(mp.weight == 0)] = False
        self._mask = self._mask.ravel()
        self._A = numpy.vstack((x_k[:,self._mask], y_k[:,self._mask], numpy.ones(len(x[self._mask]))))
        #self._A = numpy.vstack((x[self._mask], y[self._mask], numpy.ones(len(x[self._mask]))))


    def applyPrecond( self, arr ):
        """
        @brief Remove trend from map
        @param arr   map data from where to subtract the trend
        """
        z = arr.ravel()
        c = numpy.dot(self.invAA, numpy.dot(self.A, z[self.mask]))
        z[self._mask] -= numpy.dot( self._A.transpose() , c)
        arr[:] = numpy.reshape(z, numpy.shape(arr))


epsilon = 0.000001 # a small number

class divideByWeights( preconditioner ):
    """
    @brief divide each map pixel by the number of observations that fall in that pixel (a.k.a. the weight)
    """

    def __init__( self, mp ):
        """
        @param mp the Map.Map from which we'll derive the preconditioner
        """
        weight = numpy.array(mp.weight, 'float')
        weight[weight==0] = -1.
        self.invWeight = 1/weight
        self.invWeight[weight<0] = 0

    def applyPrecond( self, arr ):
        """
        @brief divides arr by self.weight
        @param arr numpy array to precondition
        """
        arr *= self.invWeight
        psLib.trace( "moby", 3, "Divide by Weight: mean = %f" % arr.mean())

class undoCommonMode( preconditioner ):
    """
    @brief  Correct for the modes lost to common mode subtraction by filtering the map
    The form of the filter in l-space is: F_l = A (l/l0)^alpha / (1+l/l0)^alpha
    The Map is FT'ed to l space and divided by F_l and FT'ed back
    """
    def __init__( self, mp, A = 1.82, alpha  =1.67, l0 = 892. ):
        """
        @param mp a Map.Map for the preconditioner
        @param A scale factor for the filter
        @param alpha The power of l in the filter
        @param l0 The pivot l for the filter

        """
        self.ra0 = mp.lon0
        self.ra1 = mp.lon1
        self.dec0 = mp.lat0
        self.dec1 = mp.lat1
        self.A = A
        self.alpha = alpha
        self.l0 = l0
 
    def applyPrecond( self, arr): 

        ftMap = fft2( arr )
                
        Ny = (arr.shape)[0]
        Nx = (arr.shape)[1]

        pixScaleX = numpy.abs(self.ra1 - self.ra0)/Nx*numpy.pi/180.\
                    *numpy.cos(numpy.pi/180.*0.5*(self.dec0+self.dec1))
        pixScaleY = numpy.abs(self.dec1-self.dec0)/Ny*numpy.pi/180. 

        lx =  2*numpy.pi  * fftfreq( Nx, d = pixScaleX ) 
        ly =  2*numpy.pi  * fftfreq( Ny, d = pixScaleY )

        ix = numpy.mod(numpy.arange(Nx*Ny),Nx) 
        iy = numpy.arange(Nx*Ny)//Nx

        lMap = ftMap*0.0
        
        lMap[iy,ix] = numpy.sqrt(lx[ix]**2 + ly[iy]**2)
        
        lOverl0 = lMap/self.l0
        
        Fl = self.A*((lOverl0)**self.alpha)/(1.+lOverl0**self.alpha)
        
        filteredFtMap = ftMap/Fl
        
        filteredFtMap[0,0] = 0.0    # avoids nan and makes mean 0
        
        arr[:,:] = (numpy.real(ifft2( filteredFtMap )))[:,:]
        
        
preconditioner_registry = {
    'meanRemover': meanRemover,
    'trendRemover': trendRemover,
    'divideByWeights': divideByWeights,
}
