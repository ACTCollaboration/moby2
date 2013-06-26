from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
Provides some of the convenience of spaceMap, but for multiple images
sharing a grid; e.g. IQU.
"""
import numpy as np

try:
    from astropy.io import fits
except ImportError:
    import pyfits as fits       # requires pyfits or astropy

from . import fits_map

def read_fits(filename, hdu=0):
    hdu_list = fits.open(filename.format(submap_code='I'))
    hdu = hdu_list[hdu]
    if hdu.header['NAXIS'] == 3:
        m = IQUMap()
        m.read_fits_iqu(filename)
    else:
        m = IQUMap()
        m.read_fits_set(filename)
    return m

class MultiMap:
    def __init__(self):
        self.data = {}
        self.images = []
        self.weights = []
        self.pixn = None
    def set_image(self, name, fits_map, copy=True):
        if len(self.images) == 0:
            self.pixn = fits_map.pixn.copy()
        if not name in self.images:
            self.images.append(name)
        if copy:
            fits_map = fits_map.copy()
        self.data[name] = fits_map
    def write_fits_set(self, filename_pattern, force=False):
        for k in self.images:
            self.data[k].write(filename_pattern.format(submap_code=k), force=force)
    def read_fits_set(self, filename_pattern, images):
        for k in images:
            m = fits_map.spaceMap(filename_pattern.format(submap_code=k))
            self.set_image(k, m, copy=False)
    def copy(self, copyData=True):
        out = self.__class__()
        for k, v in list(self.data.items()):
            out.set_image(k, v.copy(copyData=copyData))
        out.pixn = self.pixn.copy()
        return out
    def fftMap(self):
        fmap = self.__class__()
        for k in self.images:
            fmap.set_image(k, self.data[k].fftMap())
        fmap.pixn = fmap.data[self.images[0]].pixn.copy()
        return fmap
    def ifftMap(self):
        omap = self.__class__()
        for k in self.images:
            omap.set_image(k, self.data[k].ifftMap())
        omap.pixn = omap.data[self.images[0]].pixn.copy()
        return omap


class IQUMap(MultiMap):

    def read_fits_set(self, filename_pattern):
        return MultiMap.read_fits_set(self, filename_pattern, ['I','Q','U'])

    read = read_fits_set
    write = MultiMap.write_fits_set

    def extract(self, xlims, ylims, coords='sky', index=False):
        map_out = IQUMap()
        for i in self.images:
            map_out.set_image(i, self.data[i].extract(
                    xlims, ylims, coords=coords, index=index))
        return map_out

    def pretty_plot(self, filename=None, title=None, clf=None,
                    zoom=None, units=None, sticks=False, **kwargs):
        import pylab as pl

        rescale = 1.

        imshowkw = kwargs.get('imshow_args', {})

        fig, sps = pl.subplots(2,2)
        def spi(i):
            pl.sca(sps[i//2, i%2])

        codes = 'IQU'
        for i in range(3):
            spi(i)
            m = self.data[codes[i]]
            m.imshow(data=m.data*rescale, units=units, **imshowkw)
            if zoom is not None:
                m.zoom(zoom)
            pl.text(0.95, 0.95, codes[i], color='white',
                    ha='right', va='top',
                    transform=pl.gca().transAxes)
            pl.colorbar()

        # Compute polarized amplitude and direction
        Q, U = self.data['Q'].data, self.data['U'].data
        ang = np.pi/2 - np.arctan2(U, Q) / 2
        P = (Q**2+U**2)

        m0 = self.data['I']

        spi(3)
        if sticks in [True, False]:
            sticks = {'show': sticks}
        if sticks.get('show_pol'):
            m0.imshow(data=P*rescale, units=units, **imshowkw)
        else:
            m0.imshow(data=P*0, units=units, vmin=0, vmax=1, **imshowkw)

        if zoom is not None:
            self.data['I'].zoom(zoom)
        lims = pl.xlim(), pl.ylim()
        def is_in(x, lims):
            return (min(lims)<=x) and (x<=max(lims))
        pl.colorbar()
        if sticks.get('show', True):
            stick_step = sticks.get('step', 1)
            mag = P / P.max() * stick_step * 1
            dscale = fits_map.scaleParams[units]
            pitch = m0.pixn.truePitch()
            x0, y0 = m0.x[0]*dscale, m0.y[:,0]*dscale
            mx = mag*pitch[0]*dscale
            my = mag*pitch[1]*dscale
            dx = mx * np.cos(ang)
            dy = my * np.sin(ang)
            for j in range(stick_step//2, Q.shape[0], stick_step):
                if not is_in(y0[j], lims[1]): continue
                for i in range(stick_step//2, Q.shape[1], stick_step):
                    if not is_in(x0[i], lims[0]): continue
                    if mag[j,i] / mag.max() < .02: continue
                    pl.plot([x0[i]-dx[j,i],x0[i]+dx[j,i]],
                            [y0[j]-dy[j,i],y0[j]+dy[j,i]], c='white')

        pl.xlim(*lims[0]), pl.ylim(*lims[1])

        if units is None:
            units = 'deg'
        for i in [0,1]:
            sps[1,i].set_xlabel('X (%s)' % units)
            sps[i,0].set_ylabel('Y (%s)' % units)

        if title:
            pl.figtext(0.5, 0.95, title, ha='center')
        if filename:
            pl.savefig(filename)
            if clf is None:
                clf = True
        if clf:
            pl.clf()

    def write_fits_iqu(self, filename, force=False):
        """
        Write IQU map to a FITS file as a 3d HDU.
        """
        fheader = fits.Header()
        m0 = self.data[self.images[0]]

        for k,v in list(m0.pixn.header().items()):
            fheader[k] = v

        ## see ftp://ftp.aoc.nrao.edu/pub/software/aips/TEXT/PUBL/AIPSMEM117.PDF
        ## CRPIX is 1 by fiat; CRVAL and CDELT then identify the components
        ## being stepped through.  [1,2,3] corresponds to [I,Q,U].
        fheader['NAXIS'] = 3
        fheader['CTYPE3'] = 'STOKES'
        fheader['NAXIS3'] = len(self.images)
        fheader['CRPIX3'] = 1.0
        fheader['CDELT3'] = 1.0
        fheader['CRVAL3'] = 1.0

        data = np.array([self.data[k].data for k in self.images])
        fits.writeto(filename, data=data, header=fheader,
                       clobber=force)
        
    def read_fits_iqu(self, filename):
        hdu = fits.open(filename)
        fheader = hdu[0].header
        data = hdu[0].data
        assert (data.ndim == 3)

        # steal stokes axis params
        stokes_fields = ['CTYPE3', 'NAXIS3', 'CRPIX3', 'CDELT3', 'CRVAL3']
        stype, snaxis, scrpix, scdelt, scrval = [fheader[k] for k in stokes_fields]
        # remove them...
        for k in stokes_fields:
            del fheader[k]
        fheader['NAXIS'] = 2

        components = ['X', 'I', 'Q', 'U', 'V']
        for i in range(1, fheader['NAXIS3']+1):
            code = components[int((i - scrpix)*scdelt + scrval)]
            m = fits_map.fitsMap()
            m.fheader = fheader
            m.data = data[i]
            m._load_header()
            self.set_image(code, m)
