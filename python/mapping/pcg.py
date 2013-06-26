from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import os
#import psLib, numpy, actUtils
import numpy as np
import moby2.util.log as psLib

#import mobyLib
#import Mapset, preconditioner, prior
#from moby.todUtils import TOD, cuts
from moby2.tod import TOD, TODCuts
#from moby.pointing import pointingOffset
#from moby.database import result

from .preconditioner import preconditionerList, preconditioner_registry

try:
    if os.environ['MB_USE_MPI']:
        psLib.trace("moby", 3, "###################### pcg initializing MPI #########################")
        from mpi4py import MPI
        from moby.utilities import mbMPI
except:
    pass

def map2numpy( m, n, accum = False ):
    if accum:
        n[:] += m.data[:]
    else:
        n[:] =  m.data[:]

def numpy2map( n, m, accum = False ):
    if accum:
        m.data[:] += n[:]
    else:
        m.data[:] =  n[:]

def ravel_dot(a, b):
    # Where a and b have the same, multiple, dimensions.
    return (a*b).sum()


class PCG( object ):
    """
    @brief class for implementing the preconditioned conjugate gradient solution

    See (e.g.) http://en.wikipedia.org/wiki/Conjugate_gradient_method for an overview
    """

    def __init__( self, filterClass=None, filterParams=None,
                  root=None, tmp_out="./", reportFile=None,
                  clear_tod=False):
        """
        @brief Setup the pcg problem

        @param filterClass a class from scripts/pcgFilters.py
        @param filterParams dictionary of keyword arguments for filter initialization
        @param root if using MPI, the rank of the root node
        """
        self.maps = {}
        self.weights = {}
        self.nmap = 0
        self.tods = {}
        self.ntod = 0
        self.proj = {}
        self.filterClass = filterClass
        self.filterParams = filterParams
        self.clear_tod = clear_tod
        self.tolerance = 1e-4
        self.iter = 0
        self.filters    = {}   # for noise weighting in Fourier space
                                                    # have "applyFilter" methods which take a tod
        self.preconditioners = {}   # objects with "applyPrecond(map)" methods
        self.priors = {}            # objects with "applyPrior(arrIn, arrOut)" methods
        self.pcParams        = {}   # parameters for preconditioners
        self.priorParams     = {}   # parameters for priors
        self.r               = {}   # r = b - Ax
        self.b               = {}   # b = M^T N^{-1} d
        self.p               = {}   # the step in map space
        self.q               = {}   # Ap
        self.x               = {}   # The latest map


        # Temporal Outputs
        self.tmp_out = tmp_out
        self.todPlot_count = 0

        # Reporting
        if reportFile is None: self.reportFile = reportFile
        else: self.reportFile = "%s/%s"%(tmp_out, reportFile)

        # Mode-solving parameters
        self.modes = {}

        #Additional MPI parameters
        self.root = root
        if self.root is not None:
            self.rootPort = MPI.WORLD[self.root]
            self.doMPI = True
        else:
            self.doMPI = False

    def trace(self, unit, level, msg):
        psLib.trace(unit, level, msg)
        

    def addMap( self, mp, name, pcParams=None, priorParams=None ):
        """
        @brief add a map to the pcg
        @param mp Map.Map object to add
        @param name string identifying the map, must be different from other maps in the pcg
        @param pcParams optional list of preconditioner dictionaries (see mobyDefault.par)

        This should be done before adding tods which project to the added map.
        """

        self.nmap += 1
        self.maps[name] = mp
        shape = mp.data.shape
        self.nmap += 1
        self.weights[name] = np.zeros(shape)
        self.r[name] = np.zeros(shape)
        self.b[name] = np.zeros(shape)
        self.p[name] = np.zeros(shape)
        self.q[name] = np.zeros(shape)
        self.x[name] = np.array(mp.data).copy() # initial map
        self.preconditioners[name] = None
        self.priors[name] = None
        self.priorParams[name] = priorParams
        self.trace("moby", 3, "pcg: setting up preconditioners.")
        self.pcParams[name] = pcParams


    def addTOD( self, tod, projDict,
                tod_name=None,
                modesObject = None,
                setInitConditions = False, exportProjections = False, outPath = ".",
                modesLevel = 0.0001):
        """
        @brief add a tod and associated projectors to the pcg
        @param tod the TOD.TOD to  add to the pcg
        @param projDict a dictionary with map name keys and projector values
        """
        if tod_name is None:
            try:
                tod_name = tod.info.name
                if tod_name is None:
                    raise
            except:
                tod_name = str(id(tod))
                self.trace("moby", 2, "Assigning name '%s' to unknown TOD" % tod_name)
        self.ntod += 1
        self.tods[tod_name] = tod

        for mn in list(projDict.keys()):
            self.proj[ self._createProjKey( tod_name, mn ) ] = projDict[mn]
        if self.filterClass is not None:
            self.trace("moby", 3, "Estimating noise for %s." % tod_name)
            data = tod.data.copy()
            mapsZero = True
            for mapName in self.getMapNamesForTOD(tod_name):
                if self.x[mapName].max() > 0. or self.x[mapName].min() < 0.:
                    mapsZero = False
            if not mapsZero:
                self.trace("moby", 3, "Subtracting initial maps from %s." % tod_name)
                # Remove maps from data before estimating noise filters
                tod.data[:] *= -1
                self.loadMaps( weight=False )
                self.projectMapsToTOD( tod )
                tod.data[:] *= -1
            #now estimate filters
            self.filters[tod_name] = self.filterClass( tod, self.filterParams )
            tod.data[:] = data[:]
            del data
            self.filters[tod_name].setTOD(tod)
            self.filters[tod_name].applyFilter()
            self.filters[tod_name].setTOD( None )
        else:
            self.filters[tod_name] = None

        psLib.trace("moby", 3, "Initial projection of %s." % tod_name)
        self.clearMaps()
        self.projectTODToMaps(tod_name)
        mapNames = self.getMapNamesForTOD(tod_name)
        for mapName in mapNames:
            map2numpy( self.maps[mapName], self.b[mapName], accum = True )
            self.weights[mapName][:] += np.array(self.maps[mapName].weight[:])
            if exportProjections:
                filename = '%s/%s_proj.fits' % (outPath, tod_name)
                m = self.maps[mapName].trim()
                m.write(filename, weight = True)

        # If solving for correlated modes, add modesObject
        if modesObject is not None:
            # Cut sections of the TOD that fall outside the maps area.
            print("Datacopy2")
            data = tod.data.copy()
            self.cutTODexcess( tod )
            tod.data[:] = data[:]
            del data
            mo = self.modes[tod_name] = modesObject
            mo.loadModes( tod )
            nPix = 0
            for mapName in mapNames:
                mo.addMap(mapName)
                mo.b[mapName] = mo.projectModes(initialize = True, divideByWeights = False)
                nPix += self.maps[mapName].npix
            mo.setWeights(nPix = nPix)
            for mapName in mapNames:
                mo.b[mapName] /= mo.weight
                # Set initial conditions
                if setInitConditions:
                    mo.a[mapName] = mo.b[mapName].copy()
                    mo.p[mapName] = mo.a[mapName].copy() # For use in applyInverseCovariance
                    #print "Datacopy3"
                    #data = tod.data.copy()
                    mo.modes.coeff[:] = mo.b[mapName][:] * mo.weight
                    mo.modes.removeModeSet(forceFit = False, dets = None, noise = None)
                    self.clearMaps()
                    self.projectTODToMaps(tod_name)
                    map2numpy( self.maps[mapName], self.x[mapName], accum = True )
                    #tod.data[:] = data[:]
                    #del data
            mo.unloadModes()

    def setup( self, setInitMap = False ):
        """
        @brief initialize all PCG vectors and preconditioners
        Execute after you've added your maps and tods.
        """
        # Initialize report
        if self.reportFile is not None:
            self.reportInit()
        
        # If we're on a cluster, accumulate the b map
        if self.doMPI:
            self.trace("moby", 3, "pcg: Reducing b and weights")
            for mapName in list(self.maps.keys()):
                mbMPI.reduceArray( self.b[mapName], self.root )
                mbMPI.reduceArray( self.weights[mapName], self.root )
                if self.root == MPI.WORLD.rank:
                    self.maps[mapName].weight[:][:] = self.weights[mapName][:][:]
        
        # setup and apply preconditioners
        if not self.doMPI or self.root == MPI.WORLD.rank: 
            self.trace("moby", 3, "pcg: masking initial maps with the weight maps")
            for mapName in list(self.maps.keys()):
                if np.any(self.x[mapName] != 0.):
                    if setInitMap:
                        self.trace("moby", 3, "pcg: Dividing initMaps by weights.")
                        self.x[mapName] /= self.weights[mapName]
                    else:
                        self.x[mapName][:] = 0
                    self.x[mapName][np.where(self.weights[mapName] == 0)] = 0.
            self.loadMaps()
            self.trace("moby", 3, "pcg: Setting up and applying preconditioners to b.")
            for mapName in list(self.maps.keys()):
                if self.pcParams[mapName] is None:
                    continue
                pcList = preconditionerList()
                for pcDict in self.pcParams[mapName]:
                    pc_constructor = preconditioner_registry.get(
                        pcDict.get('name', None), None)
                    pc_params = pcDict.get('keywords', {})
                    if pc_constructor is None:
                        self.trace("moby", 0, "pcg: Failed to decode preconditioner spec %s" % \
                                       pcDict)
                    try:
                        pcList.append(pc_constructor(
                                self.maps[mapName], **pcDict.get('keywords',{})))
                    except:
                        self.trace("moby", 0, "pcg: Invalid preconditioner specification for %s" %\
                                mapName)
                        raise
                pcList.applyPrecond(self.b[mapName])
                self.preconditioners[mapName] = pcList

        # setup priors
        if not self.doMPI or self.root == MPI.WORLD.rank: 
            self.trace( "moby", 3, "pcg: Creating map priors" )
            for mapName in list(self.maps.keys()):
                if self.priorParams[mapName] is not None:
                    self.priors[mapName] = prior.priorList()
                    for priorDict in self.priorParams[mapName]:
                        try:
                            self.priors[mapName].append( eval(priorDict['class'])(*[self.maps[mapName]], **pcDict['keywords']) )
                        except:
                            self.trace("moby", 0, "pcg: Invalid prior specification for %s" %\
                                    mapName)
                            raise

        # get Ax  (stor in q)
        self.trace("moby", 3, "pcg: Applying inverse covariance to initial map")
        self.applyInverseCovariance( self.x, self.q )

        if self.doMPI:
            mbMPI.reduceArray( self.q[mapName], self.root )

        if not self.doMPI or self.root == MPI.WORLD.rank: 

            for mapName in list(self.maps.keys()):                    
                if self.preconditioners[mapName] is not None:   # Q M^T N^-1 M x + Q K^-1 x
                    if np.any(self.q[mapName]):
                        self.preconditioners[mapName].applyPrecond(self.q[mapName])

            for mapName in list(self.maps.keys()):
                self.r[mapName] = self.b[mapName] - self.q[mapName]
                self.p[mapName] = self.r[mapName].copy()
                if len(self.modes) > 0:
                    for tod_name in list(self.tods.keys()):
                        if mapName in self.getMapNamesForTOD(tod_name):
                            mo = self.modes[tod_name]
                            mo.r[mapName] = mo.b[mapName] - mo.q[mapName]
                            mo.p[mapName] = mo.r[mapName].copy()

        if self.doMPI:
            for mapName in list(self.maps.keys()):
                mbMPI.broadcastArray( self.p[mapName], self.root )


    def step( self ):
        """
        Iterate the solver. 

        @return True if converged
        """
        #First check convergence
        finished = False
        if not self.doMPI or self.root == MPI.WORLD.rank:
            conv1 = self.convergence()
            if self.tolerance > conv1:
                if self.root is None:
                    return True
                finished = True

        if self.doMPI:
            finished = self.rootPort.Bcast( finished )
            if finished:
                return True

        # Increment iteration counter
        self.iter += 1

        # q = Ap
        self.trace("moby", 3, "pcg step: Applying inverse covariance.")
        self.applyInverseCovariance( self.p, self.q )

        if self.doMPI:
            for mapName in list(self.maps.keys()):
                mbMPI.reduceArray( self.q[mapName], self.root )

        if not self.doMPI or self.root == MPI.WORLD.rank: 
            for mapName in list(self.maps.keys()): # Q M^T N^-1 M x + Q K^-1 x
                if self.preconditioners[mapName] is not None:
                    self.preconditioners[mapName].applyPrecond(self.q[mapName])
            
            #compute alpha and r^T r
            rTr1 = 0; pTAp1 = 0; rTr2 = 0; pTAp2 = 0
            for mapName in list(self.maps.keys()):
                rTr1  += ravel_dot(self.r[mapName], self.r[mapName])
                pTAp1 += ravel_dot(self.p[mapName], self.q[mapName])  # q = Ap
                if len(self.modes) > 0:
                    for tod_name in list(self.tods.keys()):
                        if mapName in self.getMapNamesForTOD(tod_name):
                            mo = self.modes[tod_name]
                            rTr2 += ravel_dot(mo.r[mapName], mo.r[mapName])
                            pTAp2 += ravel_dot(mo.p[mapName], mo.q[mapName])
            rTr = rTr1 + rTr2
            pTAp = pTAp1 + pTAp2
            alpha = rTr/pTAp
            self.trace("moby", 4, "pcg: rTr1 = %f\n     rTr2 = %f"%(rTr1, rTr2))
            self.trace("moby", 3, "pcg: rTr = %f"%rTr)
            self.trace("moby", 4, "pcg: pTAp1 = %f\n     pTAp2 = %f"%(pTAp1, pTAp2))
            self.trace("moby", 3, "pcg: pTAp = %f"%pTAp)
            self.trace("moby", 2, "pcg: alpha = %f"%alpha)
    
            #compute new x and r
            for mapName in list(self.maps.keys()):
                self.x[mapName] += alpha*self.p[mapName]
                self.r[mapName] -= alpha*self.q[mapName]
                if len(self.modes) > 0:
                    for tod_name in list(self.tods.keys()):
                        if mapName in self.getMapNamesForTOD(tod_name):
                            mo = self.modes[tod_name]
                            mo.a[mapName] += alpha*mo.p[mapName]
                            mo.r[mapName] -= alpha*mo.q[mapName]
    
            #compute beta
            rTrNew1 = 0; rTrNew2 = 0
            for mapName in list(self.maps.keys()):
                rTrNew1 += ravel_dot(self.r[mapName], self.r[mapName])
                if len(self.modes) > 0:
                    for tod_name in list(self.tods.keys()):
                        if mapName in self.getMapNamesForTOD(tod_name):
                            mo = self.modes[tod_name]
                            rTrNew2 += ravel_dot(mo.r[mapName], mo.r[mapName])
            rTrNew = rTrNew1 + rTrNew2
            beta = rTrNew/rTr
            self.trace("moby", 4, "pcg: rTrNew1 = %f\n     rTrNew2 = %f"%(rTrNew1, rTrNew2))
            self.trace("moby", 3, "pcg: rTrNew = %f"%rTrNew)
            self.trace("moby", 2, "pcg: beta = %f"%beta)
    
            #compute the next step
            for mapName in list(self.maps.keys()):
                self.p[mapName] = self.r[mapName] + beta * self.p[mapName]
                if len(self.modes) > 0:
                    for ts in list(self.tods.keys()):
                        if mapName in self.getMapNamesForTOD(tod_name):
                            mo = self.modes[tod_name]
                            mo.p[mapName] = mo.r[mapName] + beta * mo.p[mapName]

        if  self.doMPI:
            for mapName in list(self.maps.keys()):
                mbMPI.broadcastArray( self.p[mapName], self.root )

        if self.reportFile is not None:
            conv2 = self.convergence()
            dat = [conv1, rTr1, rTr2, rTr, pTAp1, pTAp2, pTAp, alpha, 
                   rTrNew1, rTrNew2, rTrNew, beta, conv2]
            self.reportStep(dat)

        return False
        

    def applyInverseCovariance(self, ps, qs):
        """
        @brief compute qs = A ps = M^T N^-1 M ps

        Uses TODs and maps as scratch.

        @param ps a list of pcg vectors to apply the inverse covariance to
        @param qs a list of pcg vectors to accumulate the output in
        """
        # Initalize accumator to 0
        for q in list(qs.values()):
            q[:] = 0.

        for mapName in list(self.maps.keys()):
            if self.priors[mapName] is not None:  # q = K^-1 x
                self.priors[mapName].applyPrior(ps[mapName], qs[mapName])

        for tod_name in list(self.tods.keys()):
            tod = self.getTOD(tod_name)
            self.loadMaps( arr = ps, weight = False )
            self.projectMapsToTOD(tod_name)                       # M x + m a
            if self.filters[tod_name] is not None:
                self.filters[tod_name].setTOD(tod)
                self.filters[tod_name].applyFilter()              # N^-1 M x 
                self.filters[tod_name].setTOD(None)
            self.clearMaps()
            self.projectTODToMaps(tod_name)                       # M^T N^-1 M x + m^T...
            self.unloadMaps( arr = qs, accum = True )             # q = M^T N^-1 M x + K^-1 x
            del tod
        

    def projectMapsToTOD(self, tod_name):
        """
        @brief project all maps associated with a tod into that tod
               Also projects the modes solution "p" onto the TOD
        @param tod TOD.TOD associated with this pcg
        """
        mapNames = self.getMapNamesForTOD(tod_name)
        tod = self.getTOD(tod_name)
        for mapName in mapNames:
            projName = self._createProjKey(tod_name, mapName)
            m = self.maps[mapName]
            if not np.any(m.data):
                continue
            proj, cuts1, cuts2 = self.proj[projName]
            proj.projectMapToTOD(tod, m, cuts=cuts1)
            # Does this work?
            if tod_name in self.modes:
                mo = self.modes[tod_name]
                mo.loadModes(tod)
                mo.deProjectModes(mo.p[mapName])
                mo.unloadModes()
            
    def projectTODToMaps(self, tod_name):
        """
        @brief Project a TOD into its associated maps
        @param tod TOD.TOD associated with this pcg
        """
        tod = self.tods[tod_name]
        mapNames = self.getMapNamesForTOD(tod_name)
        for mapName in mapNames:
            projName = self._createProjKey(tod_name, mapName)
            m = self.maps[mapName]
            proj, cuts1, cuts2 = self.proj[projName]
            proj.projectTODToMap(m, tod, cuts=cuts2)
            if self.clear_tod:
                del tod.data

    def convergence( self ):
        """
        @brief Compute the ratio of the residual over the projected data

        @return |r|/|b|
        """
        num = 0.0
        den = 0.0
        for mapName in list(self.maps.keys()):
            r = self.r[mapName]
            b = self.b[mapName]
            num += ravel_dot(r,r)
            den += ravel_dot(b,b)
            if len(self.modes) > 0:
                for tn in list(self.modes.keys()):
                    r = self.modes[tn].r[mapName]
                    b = self.modes[tn].b[mapName]
                    num += ravel_dot(r,r)
                    den += ravel_dot(b,b)
        return (num/den)**0.5


    def setWeights( self ):
        """
        @brief sets the weights for all maps
        """
        self.trace("moby", 3, "pcg: setting weights")
        self.clearMaps()
        for tod_name in list(self.tods.keys()):
            tod = self.getTOD(tod_name)
            self.projectTODToMaps(tod)
            del tod
        for mapName in list(self.maps.keys()):
            self.weights[mapName] = np.array(self.maps[mapName].weight).copy()

    def getTOD(self, tod_name):
        tod = self.tods[tod_name]
        if tod.data is None:
            tod.data = np.zeros((len(tod.det_uid), tod.nsamps),
                                    'float32')
        else:
            tod.data[:] = 0.
        return tod

    def loadMaps( self, arr = None, weight = True ):
        """
        @brief load current map (self.x[i]) into the mapset map (self.mapsets[i].map)
        @param arr if specified, load different arrays. for instance loadMaps(pcg.b) will give you 
          raw the filtered maps
        @param weight bool: load weights too
        In general, the maps are used as scratch space.
        """  
        if arr is None:
            arr = self.x
        for mapName in list(self.maps.keys()):
            if weight:
                if self.weights[mapName] is None:
                    self.setWeights()
                self.maps[mapName].weight[:][:] = self.weights[mapName][:][:]
            numpy2map( arr[mapName], self.maps[mapName] )


    def unloadMaps( self, arr = None, accum = False ):
        """
        @brief Unload maps into an array
        @param dictionary of arrays to unload into (defaults to self.x)
        @param accum bool: add into array instead of using assignment
        """
        if arr is None:
            arr = self.x

        for mapName in list(self.maps.keys()):
            map2numpy( self.maps[mapName], arr[mapName], accum = accum )


    def clearMaps(self):
        for m in list(self.maps.values()):
            m.clear()


    def _createProjKey( self, todName, mapName ):
        """
        @brief construct a key for the projector dictionary out of the associated map and tod names
        @param todName name of tod
        @param mapName name of map
        @return projector key string
        """
        return "%s.%s" % (todName, mapName)


    def getMapNamesForTOD( self, todName ):
        """
        @brief Get all map names associated with a particular TOD
        """
        todMapNames = []
        allMapNames  = list(self.maps.keys())
        projNames = list(self.proj.keys())
        for mapName in allMapNames:
            projName = self._createProjKey(todName, mapName)
            if projName in projNames:
                todMapNames += [mapName]
        return todMapNames       


    def cutTODexcess( self, tod ):
        """
        @brief Cut sections of the TOD that fall outside the maps area. 
        This destroys the data in the tod and in the maps.
        """
        print('cutTODexcess is in use.')
        mapNames = self.getMapNamesForTOD( tod.name )
        for mapName in mapNames:
            arr = self.x[mapName].copy()
            arr[:][:] = 1.0
            self.loadMaps(arr = {mapName: arr}, weight = False)
            projName = self._createProjKey(tod.name, mapName)
            m = self.maps[mapName]
            p = self.proj[projName]
            tod.data[:][:] = 0
            mobyLib.projectMapToTOD( tod.ctod, m.cmap, tod.cuts.cCuts, p )
            tod.cuts.cutTODZeros(tod)
        self.clearMaps()


    def exportTODplot( self, tod, tag = "tod_deproj"):
        """
        @brief  Plots a single detector TOD with cut regions in red
        """
        d,r,c = tod.listUncut()
        p = tod.plotWithCuts(r[0],c[0], show = False)
        p.save_as_eps("%s/%s_%d.eps"%(self.tmp_out, tag, self.todPlot_count))
        self.todPlot_count += 1


    def exportFits( self, dat, outPath, tag = "export", trimMap = False):
        """
        @brief Export a map (or any intermadiate map) as a fits file 
        """
        self.loadMaps(arr = dat, weight = True)
        for mapName in list(self.maps.keys()):
            filename = '%s/%s_%s.fits' % (outPath, mapName, tag)
            if trimMap:
                m = self.maps[mapName].trim()
            else:
                m = self.maps[mapName]
            m.write(filename, weight = True)


    def reportStep( self, dat ):
        """
        @brief Write internal step variables to file
        """
        keys = ["Conv. before", "rTr1", "rTr2", "rTr", "pTAp1", "pTAp2", "pTAp", "Alpha",
                "rTrNew1", "rTrNew2", "rTrNew", "Beta", "Conv. after"]
        f = open(self.reportFile, "a")
        f.write("\nStep %d:\n"%self.iter)
        for i in range(len(keys)):
            f.write("%s: %e\n"%(keys[i], dat[i]))
        f.close()

    def reportInit( self ):
        """
        @brief  Initialize stepping report
        """
        f = open(self.reportFile, "w")
        f.write("### Stepping report for PCG ###\n\n")
        f.write("Map names:\n")
        for name in list(self.maps.keys()):
            f.write("    %s\n"%name)
        f.write("\nNumber of TODs: %d\n\n"%self.ntod)
        f.close()

        
class PCGFilter( object ):
    """
    @brief A class to describe our pcg time domain filters

    Implements (1-C)N^{-1}(1-C) where 1-C is the correlation mode removal and
    N^{-1} is the independent noise filtering. Both operations must be linear
    symmetric and non-singular.
    """

    def __init__( self, noiseTOD, filterNoise = True ):
        """
        @param noiseTOD TOD.TOD with best estimate of noise
        @param filterNoise bool if False just perform decorrelation
        """
        pass

    def setTOD( self, tod ):
        """
        @brief Set the tod associated with filters 
        @param tod TOD.TOD to point the filters to
        """
        pass

    def removeCorrelations( self ):
        """
        @brief remove correlation modes (1-C)
        """
        pass

    def noiseFilter( self ):
        """
        @brief filter individual detector noises
        """
        pass

    def applyFilter(self, data=None):
        """
        @brief Filter the TOD. Performs the pcg (1-C)^TN^{-1}(1-C) step.
        """
        psLib.trace( "moby", 2, "Filtering")
        self.removeCorrelations(data) # don't need to repeat 1-C b/c (1-C)^T(1-C) = 1-C
        if self.params['filterNoise']:
            psLib.trace( "moby", 3, "Noise Filtering")
            self.noiseFilter(data)
        if not(self.params['useNoiseInCorr']):
            self.removeCorrelations(data)


class modesObject( object ):
    """
    @brief  class for implementing modes solution in the context of the PCG.
    """

    def __init__( self, depot, tag, level ):
        """
        @param  modes  object containing the modes to be solved.
        """
        pass

    def loadModes( self, tod ):
        """
        @brief Loads a modeSet for a given tod
        """
        pass

    def unloadModes( self ):
        """
        @brief Clears modes to reduce memory usage
        """
        pass

    def projectModes( self ):
        """
        @brief  Excecute the "m^T" operation (dot product between modes and TOD).
        """
        pass

    def deProjectModes( self, coeff = None ):
        """
        @brief  Excecute the "m" operation (generate a virtual TOD from the modes weighted by
        the solution coefficients)
        @param  coeff  coefficients to use for de-projection
        """
        pass
