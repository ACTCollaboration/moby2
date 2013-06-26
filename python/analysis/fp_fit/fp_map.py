from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import numpy as np

class FPMap:
    def __init__(self, x=None, y=None, sig=None, n=100, dx=None):
        self.alt_rate = 0.
        if x is not None:
            self.setup_xy(x, y, n=n, dx=dx)
            if  sig is not None:
                self.make_map(x, y, sig)
        else:
            if sig is not None:
                print('Pass x,y with sig to force instant map-making.')

    def setup_tod(self, tod, n=100, dx=None):
        self.n = [n, n]
        self.alt0 = np.mean(tod.alt)
        self.alt_rate = 15.*np.cos(self.alt0) * np.pi/180 / 3600
        x = tod.az
        y = tod.alt + tod.dt*self.alt_rate
        self.setup_xy(x, y, n, dx)
        self.clear()
        return y

    def setup_xy(self, x, y, n=100, dx=None):
        if dx is None:
            self.n = [n, n]
            self.dx = [(max(xx)-min(xx)) / self.n[i]
                       for i,xx in enumerate([x,y])]
        else:
            self.dx = dx
            self.n = [round((max(xx)-min(xx)) / self.dx[i] + 0.5)
                      for i,xx in enumerate([x,y])]
        self.x0 = [min(x), min(y)]
        self.clear()

    def clear(self):
        self.map = np.zeros(tuple(self.n))
        self.wts = np.zeros(tuple(self.n))

    def bin(self, xi, yi, sig):
        bins = (np.arange(-0.5, self.n[0]+0.4), np.arange(-0.5, self.n[1]+0.4))
        new_map, _,_ = np.histogram2d(xi, yi, bins=bins, weights=sig)
        new_wts, _,_ = np.histogram2d(xi, yi, bins=bins)
        self.map = self.map*self.wts + new_map
        self.wts += new_wts
        w = (self.wts > 0)
        self.map[w] /= self.wts[w]
        return new_map
        
    def make_map_t(self, az, alt, dt, sig):
        y = alt + dt*self.alt_rate
        xi = ((az - self.x0[0]) / self.dx[0]) #.astype('int')
        yi = ((y  - self.x0[1]) / self.dx[1]) #.astype('int')
        self.bin(xi, yi, sig)

    def make_map(self, x, y, sig):
        xi = ((x - self.x0[0]) / self.dx[0]) #.astype('int')
        yi = ((y - self.x0[1]) / self.dx[1]) #.astype('int')
        return self.bin(xi, yi, sig)

    def peak_position(self):
        mx = np.unravel_index(np.argmax(self.map), self.map.shape)
        c = [self.x0[i]+mx[i]*self.dx[i] for i in range(2)]
        if self.alt_rate != 0:
            #Convert y to time.
            c[1] = (c[1] - self.alt0) / self.alt_rate
        return tuple(c) + (self.map[mx],)
        
    def imshow(self, *args, **kwargs):
        """
        See fitsMap.imshow.
        """
        import moby2
        lims = []
        for i in [0,1]:
            lims.append((self.x0[i], self.x0[i] + (self.n[i]-.9)*self.dx[i]))
        map1 = moby2.mapping.fits_map.simpleMap(lims[0], lims[1], self.dx)
        map1.data[:self.n[1],:self.n[0]] = self.map.transpose()
        im = map1.imshow(*args, **kwargs)
        return im


def find_source(x, y, sig, res=0.001):
    """
    Makes a rough guess at the position and brightness of a single source by
    creating a coarse map.
    """
    # Guess planet location based on coarse alt-az map
    s = fpMap(x=x, y=y, sig=sig, dx=[res,res])
    x0, y0, h = s.peak_position()
    return x0, y0, h


if __name__ == '__main__':
    # Through coarse projection, find the position of a source in a TOD.
    from optparse import OptionParser
    import moby2
    
    o = OptionParser()
    o.add_option('-d','--det', action='append', type='int', nargs=2)
    opts, args = o.parse_args()

    r_list, c_list = [r for r,_ in opts.det], [c for _,c in opts.det]

    for filename in args:
        tod = moby2.TOD.from_dirfile(filename, rows=r_list, cols=c_list)
        for d in tod.data:
            print('Actually, let us not.')
            
