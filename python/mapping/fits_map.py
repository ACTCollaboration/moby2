from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
fitsMap and friends.

numpy-based handling of FITS maps and some WCS.
"""

import os,sys

try:
    from astropy.io import fits
except ImportError:
    import pyfits as fits       # requires pyfits or astropy

import numpy as np

# Assistance

def _to_deg(x):
    return (180./np.pi)*x

def _to_rad(x):
    return (np.pi/180.)*x

def _as_arrays(*args, **kwargs):
    return [np.asarray(x, **kwargs) for x in args]

def euler(alpha, axis):
    o, z = np.ones(np.shape(alpha)), np.zeros(np.shape(alpha))
    c, s = np.cos(alpha), np.sin(alpha)
    M = np.array([[o,  z,  z],
               [z,  c, -s],
               [z,  s,  c]])
    idx = [(3-axis)%3,(4-axis)%3,(5-axis)%3]
    M = M[idx][:,idx]
    return M


MAP_DEFAULT_DTYPE = 'float64'
MAP_DEFAULT_WTYPE = 'int64'
MAP_DEFAULT_FOURIERNORMSCHEME = 'integral'

    
"""
FITS coordinate systems

To go from a sky coordinate to a particular map:

1. In: celestial spherical coordinates

2. Rotate to native spherical coordinates using CRVAL, LONPOLE, LATPOLE.

3. Project to intermediate world coordinates using CTYPE, PV.

4. Transform to pixel coordinates (relative to the reference pixel),
   using CDELT, PC.

5. Transform to map pixel coordinates using CRPIX.


Transformations 1-4 are called a pixelization.  Transformation 5
should be described along with the bounds of the map (NAXISi).

Transformations 2 and 3 are closely linked.  The form of
transformation 4 is fixed.  Transformation 5 can only be done with a
particular map in mind.

Our base classes to handle this are thus

fitsProjection - takes celestial spherical coordinates to intermediate
  world coordinates.

fitsGrid - takes intermediate world coordinates to pixel coordinates.

fitsMapGrid - takes pixel coordinates to map pixel coordinates.

fitsPixelization - groups a fitsProjection, fitsGrid, and possibly a
  fitsMapGrid together.


Our coordinate systems are abbreviated:
   'sky'    celestial spherical coordinates
   'xy'     intermediate world coordinates
   'ij'     pixel coordinates (column and row, with reference pixel at 0,0)
   'IJ'     map pixel coordinates (column, row), zero-indexed.

Note that the data arrays in fitsMap (.data, .weight) are still
indexed by (row, column) -- [J,I].
   
"""

class fitsHeaderData(dict):
    def __init__(self, keys=[], header={}):
        dict.__init__(self)
        self.key_order = keys
        self.values = {}
        for k in keys:
            self[k] = header.get(k, None)
    def copy(self):
        return self.__class__(header=self)
    def items(self):
        all = list(self.keys())
        out = []
        for k in self.key_order:
            out.append((k, self[k]))
            if k in all: all.remove(k)
        for k in all:
            out.append((k, self[k]))
        return out

class fitsProjection(fitsHeaderData):
    def __init__(self, header={}, keys=[]):
        keys = ['CRVAL1','CRVAL2','LONPOLE0','LATPOLE0','CTYPE1','CTYPE2'] + keys
        fitsHeaderData.__init__(self, header=header, keys=keys)

    @staticmethod
    def getFrom(header=None):
        ctype = header['CTYPE1']
        if ctype.endswith('-CEA'):
            return CEAProjection(header)
        if ctype.endswith('-TAN'):
            return TANProjection(header)
        if ctype.endswith('-CAR'):
            return CARProjection(header)
        if ctype.endswith('-GLS'):
            print('Warning: GLS coordinates not properly supported...')
            return TANProjection(header)
        if ctype == 'X' or ctype == 'RA---CAR':
            return LinearProjection(header)
        raise ValueError('unknown projection type "%s"' % ctype)

    # virtual methods.
    def xy_from_sky(self, phi, theta):
        pass
    def sky_from_xy(self, x, y):
        pass

    def arcScales(self, x=None, y=None):
        """
        Return arclength derivatives with respect to x and y at x and y.
        """
        return 1., 1.

class CEAProjection(fitsProjection):
    """
    Provides the conversions between native spherical and intermediate
    world coordinates for the Cylindrical Equal Area projection.
    """
    def __init__(self, header={}, crval=(0., 0.), lamb=1., confLat=None):
        """
        Create a CEA projection. Pass in conformal latitude (confLat)
        in degrees, or lambda parameter (lamb).
        """
        if confLat is not None:
            lamb = np.cos(_to_rad(confLat))**2
        _header = {
            'CRVAL1': crval[0],
            'CRVAL2': crval[1],
            'PV2_1': lamb,
            'CTYPE1': 'RA---CEA',
            'CTYPE2': 'DEC--CEA',
            'LONPOLE0': 0.,
            'LATPOLE0': 0.,
            }
        if header is not None:
            _header.update(dict(list(header.items())))
        fitsProjection.__init__(self, header=_header, keys=['PV2_1'])
        # Warn of our limitations...
        if self['CRVAL2'] != 0.:
            raise RuntimeError("%s is not set up to handle non-zero CRVAL2!" % \
                  str(self.__class__))

    def xy_from_sky(self, phi, theta):
        phi, theta = _as_arrays(phi, theta)
        return phi - self['CRVAL1'], _to_deg(np.sin(_to_rad(theta))) / self['PV2_1']

    def sky_from_xy(self, x, y):
        x, y = _as_arrays(x, y)
        return x + self['CRVAL1'], _to_deg(np.arcsin(_to_rad(y)*self['PV2_1']))

    def cdelt(self, pixelSize):
        """
        Suggest FITS CDELT values for this projection, based on the
        desired pixel size.  (Pixels will be square at the conformal
        latitude.)

        Returns tuple (CDELT1, CDELT2).
        """
        p = pixelSize / self['PV2_1']**0.5
        return -p, p

    def arcScales(self, x=None, y=None):
        if x is None or y is None:
            phi, theta = self.sky_from_xy(0., 0.)
        else:
            phi, theta = self.sky_from_xy(x, y)
        # Analytic shifts -- partial derivatives
        ##  d(phi cos theta) = cos(theta) dx + (...)
        ##  d(theta)         = lambda / cos(theta) * dy
        return np.cos(_to_rad(theta)), self['PV2_1'] / np.cos(_to_rad(theta))
        

class CARProjection(fitsProjection):
    """
    Provides the conversions between native spherical and intermediate
    world coordinates for the Plate Carree projection.
    """
    def __init__(self, header={}, crval=(0., 0.)):
        """
        Create a CAR projection.  The reference pixel must be on the
        celestial equator.
        """
        _header = {
            'CRVAL1': crval[0],
            'CRVAL2': crval[1],
            'CTYPE1': 'RA---CAR',
            'CTYPE2': 'DEC--CAR',
            }
        if header is not None:
            _header.update(dict(list(header.items())))
        fitsProjection.__init__(self, header=_header)
        # Warn of our limitations...
        if self['CRVAL2'] != 0.:
            raise RuntimeError("%s is not set up to handle non-zero CRVAL2!" % \
                  str(self.__class__))

    def xy_from_sky(self, phi, theta):
        phi, theta = _as_arrays(phi, theta)
        return phi - self['CRVAL1'], theta - self['CRVAL2']

    def sky_from_xy(self, x, y):
        x, y = _as_arrays(x, y)
        return x + self['CRVAL1'], y + self['CRVAL2']

    def cdelt(self, pixelSize):
        """
        Suggest FITS CDELT values for this projection, based on the
        desired pixel size.  (Pixels will be square at the reference
        latitude.)

        Returns tuple (CDELT1, CDELT2).
        """
        p1 = -pixelSize / np.cos(_to_rad(self['CRVAL2']))
        p2 = pixelSize 
        return p1, p2

    def arcScales(self, x=None, y=None):
        if x is None or y is None:
            phi, theta = self.sky_from_xy(0., 0.)
        else:
            phi, theta = self.sky_from_xy(x, y)
        # Analytic shifts -- partial derivatives
        ##  d(phi cos theta) = cos(theta) dx + (...)
        ##  d(theta)         = lambda / cos(theta) * dy
        return np.cos(_to_rad(theta)), 1.
        

class TANProjection(fitsProjection):
    """
    Provides the conversions between native spherical and intermediate
    world coordinates for the gnomic (tangent plane) projection.
    """
    def __init__(self, header={}, crval=(0., 0.)):
        """
        Create a gnomic (TAN) projection.
        """
        _header = {
            'CRVAL1': crval[0],
            'CRVAL2': crval[1],
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
#            'LONPOLE0': 0.,
#            'LATPOLE0': 0.,
            }
        if header is not None:
            _header.update(dict(list(header.items())))
        fitsProjection.__init__(self, header=_header)
        self.compute()

    def compute(self):
        phi0, theta0 = _to_rad(self['CRVAL1']), _to_rad(self['CRVAL2'])
        # Let R be the rotation that takes [0,0,1] to (phi,theta) = CRVAL,
        # and takes (for small e)          [0,e,1] to (phi,theta+e)
        # i.e. it takes orientation along +y to orientation along +theta.
        # Then R consists of first a rotation about the x axis, then about
        # the z axis.  The first orientation takes [0,0,1] to theta by
        # rotating by pi/2-theta.  This takes the zenith into the yz
        # plane, so the next rotation must include an extra pi/2.
        self.R = np.dot(euler(np.pi/2+phi0, 2), euler(np.pi/2-theta0, 0))
        self.Rinv = np.dot(euler(-np.pi/2+theta0, 0), euler(-np.pi/2-phi0, 2))

    def xy_from_sky(self, phi, theta):
        # Generate x,y,z on celestial sphere.
        phi, theta = _to_rad(np.asarray(phi)), _to_rad(np.asarray(theta))
        x, y, z = np.cos(phi)*np.cos(theta), np.sin(phi)*np.cos(theta), np.sin(theta)
        # Rotate reference pixel to zenith
        x1, y1, z1 = np.dot(self.Rinv, np.array([x,y,z]))
        # Project to tangent plane and rescale
        return _to_deg(x1/z1), _to_deg(y1/z1)

    def sky_from_xy(self, x, y):
        # Move tangent plane x,y onto unit sphere near zenith.
        x, y = _to_rad(np.asarray(x)), _to_rad(np.asarray(y))
        z = (1.+x**2+y**2)**.5   # this is actually the distance...
        x, y, z = x/z, y/z, 1/z  # on shell
        # Rotate zenith to reference pixel
        x1, y1, z1 = np.dot(self.R, (x,y,z))
        # Convert to polar coordinates
        phi = np.arctan2(y1, x1)
        theta = np.arcsin(z1)
        return _to_deg(phi), _to_deg(theta)
        
    def cdelt(self, pixelSize):
        """
        Suggest FITS CDELT values for this projection, based on the
        desired pixel size.  (Pixels will be square at the conformal
        latitude.)

        Returns tuple (CDELT1, CDELT2).
        """
        p = pixelSize
        return -p, p

class LinearProjection(fitsProjection):
    """
    Provide access to linear coordinates.
    """
    def __init__(self, header={}, crval=(0., 0.)):
        """
        Create a projection suitable for use in linear coordinates.
        """
        _header = {
            'CRVAL1': crval[0],
            'CRVAL2': crval[1],
            'CTYPE1': 'X',
            'CTYPE2': 'Y',
            'LONPOLE0': 0.,
            'LATPOLE0': 0.,
            }
        if header is not None:
            _header.update(dict(list(header.items())))
        fitsProjection.__init__(self, header=_header)

    def xy_from_sky(self, phi, theta):
        phi, theta = _as_arrays(phi, theta)
        return phi - self['CRVAL1'], theta - self['CRVAL2']

    def sky_from_xy(self, x, y):
        x, y = _as_arrays(x, y)
        return x + self['CRVAL1'], y + self['CRVAL2']

    def cdelt(self, pixelSize):
        """
        Suggest FITS CDELT values for this projection.  Trivial.

        Returns tuple (CDELT1, CDELT2).
        """
        return pixelSize, pixelSize

class fitsGrid(fitsHeaderData):
    def __init__(self, header=None, keys=None, cdelt=(.01, .01), pc=None):
        """
        Create a fitsGrid object.  Initializes from header (a
        fitsHeader or a dict), or set cdelt and (optionally) the
        rotation matrix PCi_j as pc.
        """
        if keys is None:
            keys = ['PC1_1','PC1_2','PC2_1','PC2_2','CDELT1','CDELT2']
        if header is None:
            header = {'CDELT1': cdelt[0],
                      'CDELT2': cdelt[1]}
            if pc is not None:
                header.update(_matrix_to_cards(pc, 'PC%i_%i'))
        else:
            header = dict(list(header.items()))
        keys = [k for k in keys if k in header]
        fitsHeaderData.__init__(self, header=header, keys=keys)
        self.compute()

    def compute(self):
        if not 'PC1_1' in self:
            self.pc, self.pc_inv = None, None
        else:
            self.pc = _cards_to_matrix(self, 'PC%i_%i')
            self.pc_inv = np.linalg.inv(self.pc)
        self.delt = np.array([self['CDELT1'], self['CDELT2']])
        return self
    
    def ij_from_xy(self, x, y):
        if self.pc_inv is not None:
            i, j = np.tensordot(self.pc_inv, np.array((x,y)), axes=[1,0])
        else:
            i, j = x, y
        i, j = i / self.delt[0], j / self.delt[1]
        return i, j

    def xy_from_ij(self, i, j):
        if self.pc is not None:
            x, y = np.tensordot(self.pc, np.array((i, j),'float'), axes=[1,0])
        else:
            x, y = np.array(i, 'float'), np.array(j, 'float')
        x, y = self.delt[0]*x, self.delt[1]*y
        return x, y

    def differentials(self, i=0, j=0):
        return self.xy_from_ij(i+1, j+1)

    def getMapGrid(self, x, y):
        """
        Given a set of points xy coordinates, return a fitsMapGri
        that includes all those points (and no other, unneccessary
        points). This is meant to work well with pixel centers.

        Note that if you want to use sky coordinates instead of xy,
        you have to convert first (using a fitsProjection).
        """
        if len(x) == 0:
            return fitsMapGrid()

        # Convert intermediate world coordinates to pixel coordinates
        i, j = self.ij_from_xy(x, y)
        i_lo, i_hi = np.floor(np.min(i)), np.ceil(np.max(i))
        j_lo, j_hi = np.floor(np.min(j)), np.ceil(np.max(j))

        # Count
        n1 = max(1, int(i_hi-i_lo))
        n2 = max(1, int(j_hi-j_lo))

        # I think the crpix computation assumes CRVAL=(0,0)...
        return fitsMapGrid(naxis=(n1, n2), crpix=(-i_lo, -j_lo))

    def decimated(self, step=(2,2), offset=(0,0)):
        """
        Return a coarsened fitsGrid.
        """
        fg = self.copy()
        fg['CDELT1'] *= step[0]
        fg['CDELT2'] *= step[1]
        fg.compute()
        return fg

    def pixelArea(self):
        """
        Returns the area, in square degrees, of a "typical" pixel.

        The meaning of "typical" varies depending on the projection.
        """
        return self['CDELT1'] * self['CDELT2']
    
    @staticmethod
    def getFrom(header=None):
        return fitsGrid(header=header)

class fitsMapGrid(fitsHeaderData):
    def __init__(self, header={}, keys=None, naxis=(1,1), crpix=(1,1)):
        _header = {
            'NAXIS': 2,
            'NAXIS1': naxis[0],
            'NAXIS2': naxis[1],
            'CRPIX1': crpix[0],
            'CRPIX2': crpix[1],
            }
        if header is not None:
            _header.update(list(header.items()))
        if keys is None:
            keys = ['NAXIS','NAXIS1','NAXIS2','CRPIX1','CRPIX2']
        fitsHeaderData.__init__(self, header=_header, keys=keys)
        
    def IJ_from_ij(self, i, j):
        i, j = _as_arrays(i,j)
        return i + self['CRPIX1'] - 1, j + self['CRPIX2'] - 1

    def ij_from_IJ(self, I, J):
        I, J = _as_arrays(I,J)
        return I + 1 - self['CRPIX1'], J + 1 - self['CRPIX2']

    def bounds(self, centers=False):
        d = 0.5
        if centers:
            d = 0.
        return np.array([-d, self['NAXIS1'] - 1 + d]), \
               np.array([-d, self['NAXIS2'] - 1 + d])

    def contains(self, I=None, J=None):
        s = True
        if I is not None:
            s = s * (I >= -0.5)*(I<self['NAXIS1'])
        if J is not None:
            s = s * (J >= -0.5)*(J<self['NAXIS2'])
        return s

    def submap(self, I, J):
        """
        Return a fitsMapGrid instance appropriate for use with new
        limits in I=[I_min,I_max) and J=[J_min,J_max+1).
        Note that this can be used for supermaps as well as submaps.
        """
        d = dict(list(self.items()))
        # Use floor and ceiling to get correct limit behaviour in the
        # case of fractional I and J
        I, J = np.sort((I,J))
        I = [np.floor(I[0] + .5), np.ceil(I[1] - .5)]
        J = [np.floor(J[0] + .5), np.ceil(J[1] - .5)]
        d['NAXIS1'] = I[1] - I[0]
        d['NAXIS2'] = J[1] - J[0]
        d['CRPIX1'] = self['CRPIX1'] - I[0]
        d['CRPIX2'] = self['CRPIX2'] - J[0]
        return fitsMapGrid(header=d)

    def decimated(self, step=(2,2), offset=(0,0)):
        """
        Return a coarsened fitsMapGrid.
        """
        fm = self.copy()
        for i in [0,1]:
            # Take, e.g., I to (I-offset)/step
            k1, k2 = 'NAXIS%i'%(i+1), 'CRPIX%i'%(i+1)
            n, cr = self[k1], self[k2]
            fm[k1] = (n - offset[i]) / step[i]  # //?
            fm[k2] = float(self[k2] - 1 - offset[i])/step[i] + 1
        return fm

    @staticmethod
    def getFrom(header=None):
        return fitsMapGrid(header=header)

    
class _converter:
    def __init__(self, parent, source, dest):
        self.data = (parent, source, dest)
    def __call__(self, c1, c2):
        p, source, dest = self.data
        fs, man = p.forward_sequence, p.managers
        forward = (fs.index(dest) >= fs.index(source))
        while source != dest:
            si = fs.index(source)
            if forward:
                att = man[si]
                _source = fs[si + 1]
            else:
                att = man[si-1]
                _source = fs[si - 1]
            if not hasattr(p, att):
                raise RuntimeError('object does not have %s attribute required to ' \
                      'convert from %s to %s coordinates' % (att, source, _source))
            att = getattr(p, att)
            c1, c2 = getattr(att, '%s_from_%s' % (_source, source))(c1, c2)
            source = _source
        return c1, c2
        

class fitsPixelization:
    """
    Contains some combination of
       projection  a fitsProjection
       grid        a fitsGrid
       mapgrid     a fitsMapGrid

    The presence of these attributes allows a fitsPixelization to
    convert between various coordinate systems.  The "convert" method
    can be used, or the special attributes

       <dest>_from_<source>(c1, c2)

    where dest and source are the two coordinate systems ('sky', 'xy',
    'ij', 'IJ').
    """

    forward_sequence = ['sky', 'xy', 'ij', 'IJ']
    managers = ['projection', 'grid', 'mapgrid']

    def __init__(self, projection=None, grid=None, mapgrid=None):
        self.projection = projection
        self.grid = grid
        self.mapgrid = mapgrid

    def header(self, h=None):
        if h is None:
            h = {}
        for a in ['projection', 'grid', 'mapgrid']:
            if hasattr(self, a):
                h.update(getattr(self,a))
        return h

    def copy(self):
        return fitsPixelization(self.projection.copy(), self.grid.copy(),
                                self.mapgrid.copy())

    def submap(self, c1, c2, coords='IJ'):
        # Convert to map coordinates.
        c1, c2 = np.sort(self.convert('IJ', coords, c1, c2))
        return fitsPixelization(self.projection.copy(), self.grid.copy(),
                                self.mapgrid.submap(c1, c2))

    def __getattr__(self, name):
        words = name.split('_')
        if len(words) == 3:
            if words[1] == 'from':
                return _converter(self, words[2], words[0])
            if words[1] == 'to':
                return _converter(self, words[0], words[2])
        raise AttributeError("bad attribute request '%s'" % name)

    def convert(self, dest, source, c1, c2):
        return _converter(self, source, dest)(c1, c2)

    def bounds(self, coords='sky', centers=False):
        """
        Get boundaries of map, in specified coordinate system.  Note
        that the map "boundaries" are defined by the pixel edges, not
        than the pixel centers.  (Pass centers=True if you want.)

        If you are getting the bounds in order to do an extraction on
        another map with the same pixelization, you'd be wise to use
        'ij' coordinates in the calls to "bounds" and "extract".

        Returns arrays (min_c1, max_c1), (min_c2, max_c2).
        """
        c1, c2 = self.mapgrid.bounds(centers=centers)
        c1, c2 = self.convert(coords, 'IJ', c1, c2)
        return np.sort(c1), np.sort(c2)

    def contains(self, c1, c2, coords='sky'):
        if coords is not None:
            c1, c2 = self.convert('IJ', coords, c1, c2)
        return self.mapgrid.contains(c1, c2)

    @classmethod
    def getFrom(cls, header=None):
        self = cls()
        self.projection = fitsProjection.getFrom(header)
        self.grid = fitsGrid.getFrom(header)
        self.mapgrid = fitsMapGrid.getFrom(header)
        return self
        
    def conjugate(self, fheader, inverse=None):
        """
        Computes fitsPixelization appropriate for the fourier
        transform of a map.  Makes a flat sky approximation
        and computes a k-coordinate system for the map.
        """
        #New header
        h = fheader.copy()
        if inverse is None:
            inverse = ('_CDELT1' in h)
        arc_keys = ['CDELT','CRVAL','CRPIX','CTYPE']
        if inverse:
            # restore
            for k in arc_keys:
                for i in ['1','2']:
                    h[k+i] = h['_'+k+i]
                    del h['_'+k+i]
        else:
            # archive
            for k in arc_keys:
                for i in ['1','2']:
                    h['_'+k+i] = h[k+i]
            # We will work with xy coordinates, since they're square.
            # Conversion to sky degrees (at some fiducial point):
            deg_x, deg_y = self.projection.cdelt(1.)
            # Bounds of map -- use edges not enters
            (x0,x1), (y0,y1) = self.bounds('xy', centers=False)
            # k-spacing
            dkx, dky = 2.*np.pi/(x1-x0)*abs(deg_x), 2.*np.pi/(y1-y0)*abs(deg_y)
            # Put k=0 in the center of the map; should match _shift_fft.
            n1, n2 = h['NAXIS1'], h['NAXIS2']
            h['CRPIX1'], h['CRPIX2'] = int(n1)//2 + 1, int(n2)//2 + 1
            h['CRVAL1'], h['CRVAL2'] = 0., 0.
            h['CDELT1'], h['CDELT2'] = dkx, dky
            h['CTYPE1'], h['CTYPE2'] = 'X', 'Y'
        return fitsPixelization.getFrom(h), h
                          
    def flattened(self, center=None, rescale=(1.,1.), recenter=None):
        """
        Returns a linearized pixelization.  This can be used to
        convert (inexactly) a given projection into a tangent plane
        equivalent.

        If center is not passed, it is taken from the center of the
        map.  The linear spacing of the pixels is estimated using the
        whole map.

        If rescale is a pair of floats, they are used to rescale the
        x- and y-axis pixel spacing.  If recenter is a pair of floats,
        this position will map (approximately) to the origin in the
        new pixelization.  If recenter==True then 'center' will map to
        the origin.

        The 'rescale' and 'recenter' parameters can be used to convert
        an equatorial coordinates map into a pseudo-tangent planet map
        centered on some important feature.
        """
        if center is None:
            center = self.bounds('sky')
            center = [0.5*(x[0]+x[1]) for x in center]
        if recenter == True:
            recenter = center
        elif recenter is None:
            recenter = (0., 0.)
        # Reference pixel
        crval = center[0]-recenter[0], center[1]-recenter[1]
        crpix = self.IJ_from_sky(*center)
        # Linearize coordinates using map edges
        I, J = self.bounds('IJ')
        x, _ = self.sky_from_IJ(I, center[1] + I*0)
        _, y = self.sky_from_IJ(center[0] + J*0, J)
        dx = (x[-1] - x[0]) / (I[-1] - I[0])
        dy = (y[-1] - y[0]) / (J[-1] - J[0])
        # Adjusted mapgrid reference pixel
        mg = self.mapgrid.copy()
        mg['CRPIX1'], mg['CRPIX2'] = crpix[0]+1, crpix[1]+1
        # New pixelization with requested reference pixel
        pixn = fitsPixelization(
            projection=LinearProjection(crval=crval),
            grid=fitsGrid(cdelt=(dx * rescale[0], dy * rescale[1])),
            mapgrid=mg)
        return pixn
    
    def decimated(self, step=(2,2), offset=(0,0)):
        grid, mapgrid = self.grid, self.mapgrid
        if grid is not None:
            grid = grid.decimated(step, offset)
        if mapgrid is not None:
            mapgrid = mapgrid.decimated(step, offset)
        return fitsPixelization(self.projection.copy(), grid, mapgrid)

    def truePitch(self, x=None, y=None, coords='sky'):
        if x is None or y is None:
            # Take map center as reference point
            x0, y0 = [np.mean(z) for z in self.bounds(coords='xy')]
        else:
            # Ensure xy coords
            x0, y0 = self.convert('xy', coords, x, y)
        # How big is a pixel?
        i0, j0 = self.convert('ij', 'xy', x0, y0)
        x1, y1 = self.convert('xy', 'ij', i0+1, j0+1)
        dx, dy = x1-x0, y1-y0
        # What is the arc scale in this projection at this point?
        dadx, dady = self.projection.arcScales(x0, y0)
        return dx * dadx, dy*dady

    def getPixelIndex(self, xy, radians=False):
        """
        Use optimized routines to get map pixel indices from input
        xy coordinates.

        xy must be an array of floats with first dimension 2.

        If radians=False, the conversion expects xy to be standard
        fitsMap xy coordinates in degrees.  If radians=True it assumes
        the input data is in radians.

        Returns indices into the map (with -1 representing out of
        bounds).
        """
        info = []
        for k in [0,1]:
            delt = self.grid['CDELT%i'%(k+1)]
            rpix = self.mapgrid['CRPIX%i'%(k+1)]
            rval = self.projection['CRVAL%i'%(k+1)]
            n = self.mapgrid['NAXIS%i'%(k+1)]
            if radians:
                delt = _to_rad(delt)
            scale = 1./delt
            origin = rval + delt*(1-rpix)
            info.append((n, scale, origin))
        s = xy.shape
        xy.shape = (2,-1)
        out = np.empty(xy.shape[1], 'int32')
        mobyLib.convertToPixelIndex(out, xy[0], xy[1],
                                    info[0][0], info[0][1], info[0][2],
                                    info[1][0], info[1][1], info[1][2])
        xy.shape = s
        out.shape = s[1:]
        return out

def _matrix_to_cards(m, keyformat):
    """
    Returns dictionary containing elements of matrix m, for use in FITS headers.
    """
    output = {}
    m = np.array(m)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            output[keyformat %(i+1,j+1)] = m[i,j]
    return output

def _cards_to_matrix(cards, keyformat):
    """
    Parse FITS header cards dictionary to extract a matrix.
    """
    n1, n2 = 0, 0
    while cards.get(keyformat %(n1+1,n2+1), None) is not None:
        n1 += 1
    while cards.get(keyformat %(n1,n2+1), None) is not None:
        n2 += 1
    m = np.array([[ cards[keyformat % (i+1, j+1)] for j in range(n2)] for i in range(n1)])
    return m

scaleParams = {
    'arcmin': 60.,
    'rad': np.pi/180.,
    'arcsec': 3600.,
    'deg': 1.,
    None: 1.,
    'ell': 180./np.pi,
    }


class fitsMap:
    """
    Base class for working with simple FITS maps.
    """

    def __init__(self, filename=None, weight=None, pixn=None,
                 dtype=MAP_DEFAULT_DTYPE,
                 wtype=MAP_DEFAULT_WTYPE,
                 fourierNormScheme=MAP_DEFAULT_FOURIERNORMSCHEME):
        self.data, self.weight, self.fheader = None, None, None
        if pixn is not None:
            self._load_pixn(pixn, dtype=dtype, wtype=wtype)
        if filename is not None:
            self.read(filename, weight=weight, dtype=dtype, wtype=wtype)
        self._fourierNormScheme = fourierNormScheme

    def _load_pixn(self, pixn, keepData=False,
                   dtype=MAP_DEFAULT_DTYPE,
                   wtype=MAP_DEFAULT_WTYPE):
        self.pixn = pixn
        n1, n2 = self.pixn.mapgrid['NAXIS1'], self.pixn.mapgrid['NAXIS2']
        self.fheader = fits.Header()
        for k, v in list(pixn.header().items()):
            self.fheader[k] = v
        if not keepData:
            self.data = np.zeros((n2,n1), dtype)
            self.weight = np.zeros((n2,n1), wtype)
        self.setXY()
        
    def _load_header(self, header=None):
        if header is not None:
            self.fheader = header
        self.pixn = fitsPixelization.getFrom(self.fheader)
        self.setXY()

    def clear(self):
        if self.data is not None:
            self.data[:] = 0.
        if self.weight is not None:
            self.weight[:] = 0.
            
    def copy(self, target=None, copyData=True):
        def get(x):
            if copyData:
                return x.copy()
            else:
                return np.zeros(x.shape, x.dtype)
        if target is None:
            target = self.__class__()
        target.pixn = self.pixn.copy()
        target.data = get(self.data)
        if self.weight is None:
            target.weight = None
        else:
            target.weight = get(self.weight)
        target._load_header(self.fheader.copy())
        target.x = self.x.copy()
        target.y = self.y.copy()
        target._fourierNormScheme = self._fourierNormScheme
        return target

    def read(self, filename, weight=None, dtype=None, wtype=None):
        hdu = fits.open(filename)
        self.fheader = hdu[0].header
        self.data = hdu[0].data
        if self.data.ndim == 3:
            # Work-around, for abscal runs on enki planet maps.
            from pixell import enmap     # You might need pixell to load 3-d maps.
            m = enmap.read_fits(filename)[0]
            self.fheader = m.wcs.to_header()
            self.fheader['NAXIS2'], self.fheader['NAXIS1'] = m.data.shape
            self.data = np.array(m, dtype=dtype)
            self._load_header()
            return

        if dtype is not None:
            self.data = self.data.astype(dtype)
        self._load_header()
        if weight in [True, None]:
            # Convert True/None to False, a valid weights filename, or an exception
            _rw = None
            if filename.endswith('.fits'):
                _rw = filename[:-4] + 'weight.fits'
                if not os.path.exists(_rw):
                    _rw = None
            if _rw is None:
                if weight is None:
                    # Fine, no weights.
                    weight = False
                else:
                    raise RuntimeError("Failed to guess weights filename for %s" % filename)
            else:
                weight = _rw
        if weight != False:
            if not os.path.exists(weight):
                raise RuntimeError("Weights file %s not found" % weight)
            hdu = fits.open(weight)
            self.weight = hdu[0].data
            if wtype is not None:
                self.weight = self.weight.astype(wtype)
        else:
            self.weight = None

    def setXY(self):
        """
        Update self.x, self.y, self.x_sign, self.y_sign from self.pixn information.
        """
        n2, n1 = self.data.shape
        c, r = np.arange(n1), np.arange(n2)
        self.x, _ = self.pixn.sky_from_IJ(c, c*0)
        _, self.y = self.pixn.sky_from_IJ(r*0, r)
        self.x = self.x.reshape(1, -1)
        self.y = self.y.reshape(-1, 1)
        dx, dy = self.pixn.grid.differentials()
        self.x_sign, self.y_sign = np.sign(dx), np.sign(dy)

    def write(self, filename, weight=None, force=False):
        #print 'Updating header keys from pixelization information:'
        header = self.fheader
        for o in [self.pixn.projection, self.pixn.grid, self.pixn.mapgrid]:
            if o is None: continue
            for k, v in list(o.items()):
                if k[0] == '_': continue
                if v is None: continue
                # if header.has_key(k):
                #    print 'Replacing %s' % k
                header[k] = v
        fits.writeto(filename, data=self.data, header=header,
                       clobber=force)
        if weight is None:
            weight = self.weight is not None
        if weight == True:
            if filename[-5:] == '.fits':
                weight = filename[:-5] + '.weight.fits'
            else:
                weight = filename + 'weight.fits'
        if weight != False:
            fits.writeto(weight, data=self.weight, header=header,
                           clobber=force)

    def _to_units(self, x, units=None):
        if hasattr(units, '__float__'):
            scale = units
        else:
            scale = scaleParams[units]
        if hasattr(x, 'float'):
            return x * scale
        else:
            return [xx * scale for xx in x]

    # Plotting

    def _checkRegularity(self):
        # Check squareness
        msg = (0, None)
        for src in self.x[0], self.y[:,0]:
            d = src[1:] - src[:-1]
            dist = abs(d.std() / d.mean())
            if dist > 0.01:
                msg = max(msg, (2,'significant distortion'))
            elif dist > 0.0001:
                msg = max(msg, (1, 'some distortion'))
        return msg

    def _meanAspect(self):
        """
        Return the mean aspect ratio of the image (typical dy / dx).
        This can be passed to imshow to keep each pixel roughly square.
        """
        dx = self.x[0,1:] - self.x[0,:-1]
        dy = self.y[1:,0] - self.y[:-1,0]
        return abs(dx.mean() / dy.mean())

    def _grid_plot_args(self, units=None, coords=None):
        # Common processing for .imshow and .contour to determine map boundaries.
        ## Extend bounds by half a pixel
        dx, dy = np.diff(self.x[0,:]), np.diff(self.y[:,0])
        x = (self.x[0,0] - dx[0]/2, self.x[-1,-1] + dx[-1]/2)
        y = (self.y[0,0] - dy[0]/2, self.y[-1,-1] + dy[-1]/2)
        x1, x2, y1, y2 = (min(x), max(x),
                          min(y), max(y))
        if self.x_sign < 0:
            x1, x2 = x2, x1
        if self.y_sign < 0:
            y1, y2 = y2, y1
        if coords is not None:
            (x1,x2),(y1,y2) = self.pixn.convert(coords,'sky',(x1,x2),(y1,y2))
        return {'extent': self._to_units([x1, x2, y1, y2], units=units)}

    def imshow(self, data=None, units=None, sig='amp', sigma=None,
               coords=None, quiet=False, ax=None, **args):
        import pylab
        err, msg = self._checkRegularity()
        if err > 0 and not quiet:
            print('Warning: %s.imshow -- %s due to curved coordinates (%s)' % \
                  (str(self.__class__), msg, str(self.pixn.projection.__class__)))
        args.update(self._grid_plot_args(coords=coords, units=units))
        if data is None:
            data = self.data
        if sig == 'amp':
            d = data
        elif sig == 'abs':
            d = abs(data)
        elif sig == 'phase' or sig == 'angle':
            d = angle(data)
        if sigma is not None:
            d0, ds = d[d!=0].mean(), d[d!=0].std()
            args['vmin'] = d0 - sigma*ds
            args['vmax'] = d0 + sigma*ds
        if 'aspect' not in args:
            args['aspect'] = self._meanAspect()
        if 'interpolation' not in args:
            args['interpolation'] = 'nearest'
        if ax is None:
            return pylab.imshow(d, origin='lower', **args)
        else:
            return ax.imshow(d, origin='lower', **args)

    def contour(self, levels, data=None, units=None,
                coords=None, quiet=False, **args):
        import pylab
        err, msg = self._checkRegularity()
        if err > 0 and not quiet:
            print('Warning: %s.contour -- %s due to curved coordinates (%s)' % \
                  (str(self.__class__), msg, str(self.pixn.projection.__class__)))
        args.update(self._grid_plot_args(coords=coords, units=units))
        if data is None:
            data = self.data
        #if not args.has_key('aspect'):
        #    args['aspect'] = self._meanAspect()
        return pylab.contour(data, levels, **args)

    def zoom(self, r, x0=0, y0=0, units=None):
        import pylab
        lims = (x0-r*self.x_sign, x0+r*self.x_sign) + \
               (y0-r*self.y_sign, y0+r*self.y_sign)
        lims = self._to_units(lims, units=units)
        pylab.xlim(*lims[0:2]) # x0-r*self.x_sign, x0+r*self.x_sign)
        pylab.ylim(*lims[2:4]) # y0-r*self.y_sign, y0+r*self.y_sign)

    def circleMask(self, r, x0=0., y0=0.):
        return (self.x-x0)**2 + (self.y-y0)**2 <= r**2

    def diagMap(self, direction='UL'):
        nx, ny = self.data.shape
        ri, ci = indices((nx,ny))
        if direction=='UL':
            return ri+ci
        elif direction=='UR':
            return ri-ci

    def extract_rc(self, row_lims, col_lims, index=False, data=None):
        """
        Creates a submap bounded by 2-element tuples row_lims and
        col_lims.  These should be non-negative begin:end indices; the
        first but not the last index will be included.  For standard
        coordinate systems, use the 'extract' method.

        If index=True, returns row and column limits for the user to
        apply to any other map-shaped arrays.
        """
        def insidify(x, lims):
            return [ max(lims[0], min(lims[1], _x)) for _x in x]

        cl = insidify([int(x) for x in col_lims], (0,self.data.shape[1]))
        rl = insidify([int(x) for x in row_lims], (0,self.data.shape[0]))
        if index:
            return rl, cl
        #map = self.__class__()
        map = self.copy(copyData=False)
        map.data = self.data[rl[0]:rl[1],cl[0]:cl[1]].copy()
        if self.weight is None:
            map.weight = None
        else:
            map.weight = self.weight[rl[0]:rl[1],cl[0]:cl[1]].copy()
        #map.fheader = self.fheader.copy()
        # Get new pixelization
        map.pixn = self.pixn.submap(cl, rl)
        # Update CRPIX, NAXIS from new mapgrid
        for k, v in list(map.pixn.mapgrid.items()):
            map.fheader[k] = v
        map.setXY()
        return map

    def extractXY(self, xlims, ylims):
        print('deprecated, use extract(...)')
        return self.extract(xlims, ylims, 'sky')

    def extract(self, xlims, ylims, coords='sky', index=False):
        clims, rlims = self.pixn.convert('IJ',coords,xlims,ylims)
        clims, rlims = np.sort(np.array((clims, rlims)))
        clims[1] += 1
        rlims[1] += 1
        return self.extract_rc(rlims, clims, index=index)

    def extract_center(self, x, y, size, coords='sky'):
        return self.extract((x-size/2,x+size/2), (y-size/2,y+size/2),
                            coords=coords)
        
    def grow(self, xlims, ylims, coords='sky', index=False):
        """
        Returns a map padded to include the limits specified.  The
        limits will be rounded to the nearest pixel center.
        """
        def insert_data(dest, src, offset):
            n1, n2 = src.shape
            dest[offset[0]:offset[0]+n1, offset[1]:offset[1]+n2] = src
        clims, rlims = self.pixn.convert('IJ',coords,xlims,ylims)
        clims, rlims = np.sort(np.array((clims, rlims)).round().astype('int'))
        clims, rlims = clims + [0,1], rlims + [0,1]
        r0, c0 = rlims[0], clims[0]
        nr, nc = rlims[1] - r0, clims[1] - c0
        map = self.__class__()
        map.data = np.zeros((nr,nc), self.data.dtype)
        insert_data(map.data, self.data, (-r0, -c0))
        if self.weight is not None:
            map.weight = np.zeros((nr,nc), self.weight.dtype)
            insert_data(map.weight, self.weight, (-r0, -c0))
        #Coordinate system
        pixn = self.pixn.submap(clims, rlims, 'IJ')
        map._load_pixn(pixn, keepData=True)
        return map

    def flattened(self, **kwargs):
        """
        Flatten projection so it looks like raw degrees.  This
        introduces a lot of distortion to most projections, but it's
        convenient for going to fourier space or making centered
        postage stamps.

        Arguments are passed on to self.pixn.flattened.
        """
        map = self.copy()
        map.pixn = self.pixn.flattened(**kwargs)
        map.fheader['CTYPE1'] = 'X'
        map.fheader['CTYPE2'] = 'Y'
        map.setXY()
        return map

    def decimated(self, step=(2,2), offset=(0,0)):
        def decimate_data(x):
            return x[offset[1]::step[1],offset[0]::step[0]].copy()
        pixn = self.pixn.decimated(step, offset)
        map = self.__class__(pixn=pixn)
        map.data = decimate_data(self.data)
        if self.weight is not None:
            map.weight = decimate_data(self.weight)
        return map

    def peak(self, mask=None, units=None):
        if mask is None:
            mask = np.ones(self.data.shape, dtype='bool')
        mm = np.argmax(self.data[mask].reshape(-1))
        J, I = [idx[mm] for idx in mask.nonzero()]
        x0, y0 = self.x[0,I], self.y[J,0]
        return tuple(self._to_units((x0, y0), units=units))

    def correctSign(self):
        mx, mn, mid = np.amax(self.data), np.amin(self.data), np.mean(self.data)
        if mid-mn > mx-mid:
            self.data *= -1.
            return -1.
        return 1.
    

    @classmethod
    def simpleMap(cls, x_lims, y_lims, spacing, **kwargs):
        """
        Return a map in linear projection that contains x_lims,
        y_lims, and has requested pixel spacing.
        """
        pixn = fitsPixelization(
            projection = LinearProjection(),
            grid = fitsGrid(cdelt=spacing))
        pixn.mapgrid = pixn.grid.getMapGrid(x_lims, y_lims)
        return cls(pixn=pixn, **kwargs)

    def radii(self, center=(0.,0.)):
        return ((self.x-center[0])**2+(self.y-center[1])**2)**0.5


class spaceMap(fitsMap):
    def fftMap(self, shift=True, rephase=False):
        f = freqMap()
        f.data = np.fft.fft2(self.data)
        f.weight = self.weight
        if f.weight is not None:
            f.weight = f.weight.copy()
        if shift:
            f.data = _shift_fft(f.data)
        f.pixn, f.fheader = self.pixn.conjugate(self.fheader)
        f.setXY()
        if rephase:
            f.data *= np.exp(-1.j*(f.x*self.x[0,0] + f.y*self.y[0,0]))
        # Apply scaling
        fns = self._fourierNormScheme
        dx, dy = self.pixn.truePitch()
        dkx, dky = f.pixn.truePitch()
        nx, ny = f.data.shape
        if fns == 'bizarre':
            f.data /= (dkx*dky)**0.5*nx*ny
        elif fns == 'integral':
            # The k=0 term is equal to the spatial integral of the map
            f.data *= abs(dx*dy)
        else:
            raise ValueError('unknown fourierNormScheme "%s"'%self.fourierNormScheme)
        f._fourierNormScheme = fns
        return f

    def pretty_plot(self, filename=None, title=None, clf=None,
                    zoom=None, units=None, imshow_args={}, **kwargs):
        import pylab as pl
        self.imshow(units=units, **imshow_args)
        pl.colorbar()
        if zoom is not None:
            self.zoom(zoom)
        if units is None:
            units = 'deg'
        pl.xlabel('X (%s)' % units)
        pl.ylabel('Y (%s)' % units)
        if title:
            pl.title(title,)
        if filename:
            pl.savefig(filename)
            if clf is None:
                clf = True
        if clf:
            pl.clf()
        
    @classmethod
    def from_enmap(cls, emap):
        from pixell import enmap
        assert(emap.ndim == 2)   # sry, kill those dims for me.
        self = cls()
        self.fheader = emap.wcs.to_header()
        self.fheader['NAXIS2'], self.fheader['NAXIS1'] = emap.data.shape
        self.data = np.array(emap)
        self._load_header()
        return self


class freqMap(fitsMap):
    def ifftMap(self, shift=True, rephase=False):
        m = spaceMap()
        m.pixn, m.fheader = self.pixn.conjugate(self.fheader)
        m.data = self.data.copy()
        m.setXY()
        # Normalization correction
        fns = self._fourierNormScheme
        dx, dy = m.pixn.truePitch()
        dkx, dky = self.pixn.truePitch()
        ny, nx = self.data.shape
        if fns == 'bizarre':
            # I'd rather not even explain this
            m.data *= (dkx*dky)**0.5*nx*ny # /= 2.*np.pi*(dx * dy)**0.5
        elif fns == 'integral':
            # The k=0 term is equal to the spatial integral of the map
            m.data /= abs(dx*dy)
        else:
            raise ValueError('unknown fourierNormScheme "%s"'%fns)
        m._fourierNormScheme = fns
        # Phasing and shifting
        if rephase:
            m.data = m.data * np.exp(1.j*(self.x*m.x[0,0] + self.y*m.y[0,0]))
        if shift:
            m.data = _shift_fft(m.data, inverse=True)
        m.data = np.fft.ifft2(m.data).real
        m.weight = self.weight
        if m.weight is not None:
            m.weight = m.weight.copy()
        return m

    def pad(self, real_tweak=-.5):
        """
        Double the maximum fourier frequency, padding with 0s, as if the real space
        map had twice its resolution.

        This will mess up the real space coordinate system.
        """
        ny, nx = self.data.shape
        #dx = self.x[0,1] - self.x[0,0]
        #dy = self.y[1,0] - self.y[0,0]
        #ri, ci = indices((ny*2,nx*2))
        #self.x = (ci-nx) * dx
        #self.y = (ri-ny) * dy
        # I don't know the cute way to do this
        new_data = np.zeros((ny*2, nx*2), self.data.dtype)
        for i in range(ny):
            new_data[i+(ny+1)//2,(nx+1)//2:(3*nx+1)//2] = self.data[i]
        self.data = new_data
        # Mess with the pixelization
        ## The mapgrid: NAXIS and CRPIX
        for target in [self.pixn.mapgrid, self.fheader]:
            target['NAXIS1'] *= 2
            target['NAXIS2'] *= 2
            target['CRPIX1'] = target['CRPIX1'] * 2
            target['CRPIX2'] = target['CRPIX2'] * 2
        ## the grid: CDELT invariant in Fourier space.
        ## the projection: does not change in any case
        # Fix archived values
        if '_CDELT1' in self.fheader:
            for i in ['1', '2']:
                self.fheader['_CDELT'+i] /= 2
                self.fheader['_CRPIX'+i] = self.fheader['_CRPIX'+i] * 2 + real_tweak
        # Refresh coordinates
        self.setXY()

    def _meanAspect(self):
        return 1.

def _shift_fft(x,inverse=False):
    n0, n1 = x.shape
    m0, m1 = (n0+1)//2, (n1+1)//2
    if inverse:
        m0, m1 = n0-m0, n1-m1
    y = x[(np.arange(n0)+m0) % n0]
    y = y[:,(np.arange(n1)+m1) % n1]
    return y

def _shift(x):
    n0, n1 = x.shape
    m0, m1 = n0//2, n1//2
    y = x*0.
    idxA = np.arange(n0)
    idxB = np.arange(m0,m0+n0) % n0
    y[idxA,:] = x[idxB,:]
    y[:,idxA] = y[:,idxB]
    return y

def shiftImage(image, dx, dy, asFFT=False):
    if not asFFT:
        imageF = image.fftMap()
    else:
        imageF = image
    rephase = np.exp(-1.j*(imageF.x*dx + imageF.y*dy))
    imageF.data *= rephase
    if not asFFT:
        return imageF.ifftMap()
    else:
        return imageF

def simpleMap(x_lims, y_lims, spacing, **kwargs):
    pixn = fitsPixelization(
        projection = LinearProjection(),
        grid = fitsGrid(cdelt=spacing))
    pixn.mapgrid = pixn.grid.getMapGrid(x_lims, y_lims)
    return spaceMap(pixn=pixn, **kwargs)

if __name__ == '__main__':
    # Test 1, load pixelization from FITS file.
    import sys, os
    filenames = sys.argv[1:]
    if len(filenames) == 0:
        filenames = [
            '/u/mhasse/store/work_100930/2008/saturn_ar1/links/1232013057.1232013075.ar1.fits',
            '/u/sievers/maps/test_act_ar1_2008_equatorial_riseset/ACT_148_2008_equatorial_rising_downweight_0.25_find_10_modes_detau_noise_dedark_debutter_noprior_no_carpets_gainly_200.fits' ]
    for filename in filenames:
        if not os.path.exists(filename): continue
        print(filename)
        fm = spaceMap(filename)
        mp = fitsPixelization.getFrom(header=fm.coords.header)
        print(' Map boundaries (pixel):  ', mp.bounds('IJ'))
        print(' Map boundaries (sky):    ', mp.bounds('sky'))
        print(' Inversion test, map center:')
        I, J = mp.bounds('IJ')
        I = I.mean()
        J = J.mean()
        print('  IJ:   ', I, J)
        r,d = mp.sky_from_IJ(I, J)
        print('  sky:  ', r, d)
        I,J = mp.IJ_from_sky(r, d)
        print('  IJ:   ', I, J)
    print()
    
    # Test 2, create new pixelization.
    # Create CEA projection at standard latitude of -53.5
    fp = CEAProjection(confLat=-53.5)
    # Fits grid with 30" pixels
    fg = fitsGrid(cdelt=fp.cdelt(30./3600))
    # A map big enough to hold given (ra, dec) corners
    x, y = list(zip(*[fp.xy_from_sky(r, d) for r,d in [(-100., -55.), (30., -50)]]))
    mg = fg.getMapGrid(x, y)
    # Examine map extent
    px = fitsPixelization(projection=fp, grid=fg, mapgrid=mg)
    print('Map corners:')
    print('      I,J              i,j                 x,y                   RA,DEC')
    for I in [0, mg['NAXIS1']-1]:
        for J in [0, mg['NAXIS2']-1]:
            i, j = px.ij_from_IJ(I, J)
            x, y = px.xy_from_IJ(I, J)
            r, d = px.sky_from_IJ(I, J)
            print('%7.1f,%7.1f  %7.1f,%7.1f  %10.4f,%10.4f  %10.4f,%10.4f' % \
                  (I,J, i,j, x,y, r,d))
            
    
