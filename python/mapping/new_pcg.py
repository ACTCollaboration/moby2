from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
import numpy as np

from .pcg_filters import PCG_filter_list
from .pol_map import IQUMap

def ravel_dot(a, b):
    # Where a and b have the same, multiple, dimensions.
    return (a*b).sum()

class PCGMap:
    """
    Base class for any portion of the PCG solution vector.
    """
    tods = None
    r = None
    p = None
    x = None

    def __init__(self, name=None):
        self.name = name
        self.tods = []

    def auto_dot(self, a, b):
        if isinstance(a, basestring):
            a = self.get_ref(a)
        if isinstance(b, basestring):
            b = self.get_ref(b)
        return ravel_dot(a, b)

    def add_TOD(self, key, tod):
        pass

    def project_from_p(self, key, data, cuts):
        pass

    def deproject_to_q(self, key, data, cuts):
        pass

    def get_qp(self):
        return self.auto_dot('q', 'p')

    def step_x(self, alpha):
        pass

    def get_rz(self):
        return self.auto_dot('r', 'z')

    def step_p(self, beta):
        pass

    def get_convergence(self):
        pass

    def get_ref(self, code):
        if code in ['x', 'r', 'q', 'z', 'p', 'weight']:
            return getattr(self, code, None)

    def set(self, code, value):
        ref = self.get_ref(code)
        if ref is None:
            setattr(self, code, value.copy())
        else:
            ref[:] = value
        
    def clear(self, code):
        vect = self.get_ref(code)
        vect[:] = 0

    def accumulate(self, code, vals):
        vect = self.get_ref(code)
        vect[:] += vals


default_pcg_params = {
    'tolerance': 1e-6,
}

class PCG:
    maps = None
    tods = None
    cuts = None
    def __init__(self, params={}):
        self.maps = []
        self.tods = []
        self.cuts = {}
        self.logger = moby2.util.log.get_logger()
        self.params = default_pcg_params.copy()
        self.params.update(params)
        self.filters = {}

    def trace(self, level, msg):
        self.logger.trace(level, msg)

    def add_map(self, map_):
        self.maps.append(map_)

    def add_TOD(self, key, tod, cuts):
        self.tods.append((key, (tod.data.shape, tod.data.dtype)))
        # Make a copy of the cuts and kill the sample index.  Where
        # we're going, we won't need the TOD's old reference point.
        assert (tod.info.sample_index == cuts.sample_offset)
        cuts = cuts.copy()
        cuts.sample_offset = 0
        self.cuts[key] = cuts
        # Deproject initial maps
        for m in self.maps:
            m.add_TOD_step1(key, tod, cuts)
        # Fill cuts before filtering.
        moby2.tod.fill_cuts(data=tod.data, cuts=self.cuts[key])
        # Filter TOD...
        filter_class = self.params.get('filter_class', 'identity')
        cls = PCG_filter_list[filter_class]
        self.filters[key] = cls(tod, self.params['filter_params'])
        self.filters[key].applyFilter(tod.data)
        # Init each map
        for m in self.maps:
            m.add_TOD_step2(key, tod.data, cuts)

    def prepare_iteration(self):
        for m in self.maps:
            m.prepare_iteration()

    def step(self):
        # 1. The slow, filtering step; compute
        #        q = pT A p
        #    For each map.
        for key, (shape, dtype) in self.tods:
            data = np.zeros(shape, dtype)
            # Project p into time domain
            for m in self.maps:
                m.project_from_p(key, data, self.cuts[key])
            # Fill cuts before filtering
            moby2.tod.fill_cuts(data=data, cuts=self.cuts[key])
            # Filter...
            self.filters[key].applyFilter(data)
            # Deproject to get q = Ap
            for m in self.maps:
                m.deproject_to_q(key, data, self.cuts[key])
            del data

        # 2. Compute alpha.
        rz, qp = 0., 0.
        for m in self.maps:
            rz += m.get_rz()
            qp += m.get_qp()
        alpha = rz / qp

        # 3. Pass alpha to all maps so they can update their solutions.
        for m in self.maps:
            m.step_x(alpha)

        # 4. Compute beta.
        new_rz = 0.
        for m in self.maps:
            new_rz += m.get_rz()
        beta = new_rz / rz

        # 5. Update conjugate vector or whatever
        for m in self.maps:
            m.step_p(beta)

        # 6. Test
        return self.get_convergence() < self.params['tolerance']
        
    def get_convergence(self):
        conv = 0
        for m in self.maps:
            conv = max(conv, m.get_convergence())
        return conv

    def write(self, filename, **kwargs):
        pass

    def plot(self, filename, **kwargs):
        pass


class PCG_MPI(PCG):

    def __init__(self, params={}):
        PCG.__init__(self, params)
        from mpi4py import MPI
        self.comm = MPI.COMM_WORLD
        self.root = (self.comm.Get_rank() == 0)

    def prepare_iteration(self):
        # Get map.r and map.weight from all peers.
        rw = [(m.get_ref('r'), m.get_ref('weight')) for m in self.maps]
        rw_all = self.comm.gather(rw, root=0)
        # Accumulate r and weight in root node, and do preconditioning.
        if self.root:
            for rw_set in rw_all[1:]:
                for (r, w), m in zip(rw_set, self.maps):
                    m.accumulate('r', r)
                    m.accumulate('weight', w)
            # Run map preconditioners on root node copy.
            PCG.prepare_iteration(self)
        # Synchronize
        self.comm.Barrier()

    def step(self):
        # 0. Copy p to nodes
        p = None
        if self.root:
            p = [m.get_ref('p') for m in self.maps]
        p_set = self.comm.bcast(p, root=0)
        if not self.root:
            # Set p, clear q.
            for m,p in zip(self.maps, p_set):
                m.set('p', p)
                m.clear('q')

        # 1. The slow, filtering step; compute
        #        q = pT A p
        #    For each map.
        for key, (shape, dtype) in self.tods:
            data = np.zeros(shape, dtype)
            # Project p into time domain
            for m in self.maps:
                m.project_from_p(key, data, self.cuts[key])
            # Fill cuts before filtering
            moby2.tod.fill_cuts(data=data, cuts=self.cuts[key])
            # Filter...
            self.filters[key].applyFilter(data)
            # Deproject to get q = Ap
            for m in self.maps:
                m.deproject_to_q(key, data, self.cuts[key])
            del data

        # 1.4 Retrieve q from nodes
        q = [m.get_ref('q') for m in self.maps]
        q_all = self.comm.gather(q, root=0)
        if self.root:
            for mi, m in enumerate(self.maps):
                for i in range(1, len(q_all)):
                    m.accumulate('q', q_all[i][mi])
            del q_all

            # 2. Compute alpha.
            rz, qp = 0., 0.
            for m in self.maps:
                rz += m.get_rz()
                qp += m.get_qp()
            alpha = rz / qp

            # 3. Pass alpha to all maps so they can update their solutions.
            for m in self.maps:
                m.step_x(alpha)

            # 4. Compute beta.
            new_rz = 0.
            for m in self.maps:
                new_rz += m.get_rz()
            beta = new_rz / rz

            # 5. Update conjugate vector or whatever
            for m in self.maps:
                m.step_p(beta)

            # 6. Test
            converged = self.get_convergence() < self.params['tolerance']
        else:
            converged = None

        # Implicit synchronization
        return self.comm.bcast(converged, root=0)
        

class Precond_DivideByWeights:
    def prepare_iteration(self, weight=None, **kwargs):
        self.weight = weight
    def apply(self, r, in_place=True):
        s = self.weight != 0
        if not in_place:
            r = r.copy()
        r[s] /= self.weight[s]
        return s

class Precond_RemoveMean:
    def prepare_iteration(self, weight=None, **kwargs):
        self.weight = weight
    def apply(self, r, in_place=True):
        s = self.weight != 0
        if not in_place:
            out = r.copy()
        else:
            out = r
        out[s] -= out[s].mean()
        return out

class _FitsMapOutputter:
    """
    Provide write/plot interface for fitsMap compatible data.
    """
    def _get_map_kw(self, kwargs):
        map0 = self.map
        kw = {}
        for k,v in list(kwargs.items()):
            if k == 'extract':
                map0 = self.map.extract(*v)
                continue
            kw[k] = v
        return map0, kw

    def write(self, filename, outputm=None, info=None, **kwargs):
        if self.is_null():
            return
        map0, kwargs = self._get_map_kw(kwargs)
        if outputm is not None:
            ofile = outputm.get_filename(filename, info)
        else:
            ofile = filename.format(**info)
        map0.write(ofile, **kwargs)

    def plot(self, filename, outputm=None, info=None, **kwargs):
        if self.is_null():
            return
        map0, kwargs = self._get_map_kw(kwargs)
        if outputm is not None:
            ofile = outputm.get_filename(filename, info)
        else:
            ofile = filename.format(**info)
        map0.pretty_plot(ofile, **kwargs)


class PCGIntensityMap(PCGMap, _FitsMapOutputter):
    def __init__(self, name=None, fits_map=None, proj_gen=None):
        PCGMap.__init__(self, name=name)
        self.map = fits_map
        self.proj_gen = proj_gen
        self.x = self.map.data
        s, d = self.x.shape, self.x.dtype
        self.r = np.zeros(s, d)
        self.q = np.zeros(s, d)
        self.weight = self.map.weight
        self.projs = {}
        self.preconds = [Precond_DivideByWeights(), Precond_RemoveMean()]
    def add_TOD_step1(self, key, tod, cuts):
        self.tods.append(key)
        self.projs[key] = self.proj_gen.get(
            map_name=self.name, tod_name=key, tod=tod,
            fits_map=self.map)
        # Project an initial solution from the TOD?
        # - This requires two full passes through all TODs, with TOD data.
        #
        # Deproject starting solution from TOD!
        if np.any(self.x):
            self.projs[key].deproject_from_map(-self.x, tod.data, cuts=cuts)
    def add_TOD_step2(self, key, data, cuts):
        # Initialize r from the TOD
        self.projs[key].project_to_map(data, output=self.r, cuts=cuts,
                                       weights_output=self.weight)
    def prepare_iteration(self):
        for precond in self.preconds:
            precond.prepare_iteration(weight=self.weight)
        self.z = self.r.copy()
        for precond in self.preconds:
            precond.apply(self.z)
        self.p = self.z.copy()
    def project_from_p(self, key, data, cuts):
        self.projs[key].deproject_from_map(self.p, output=data, cuts=cuts)
    def deproject_to_q(self, key, data, cuts):
        self.projs[key].project_to_map(data, output=self.q, cuts=cuts)
    def step_x(self, alpha):
        self.x += alpha * self.p
        self.r -= alpha * self.q
        self.z[:] = self.r
        for precond in self.preconds:
            precond.apply(self.z)
    def step_p(self, beta):
        self.p[:] = self.p*beta + self.z
        self.q[:] = 0
    def get_convergence(self):
        D = (self.x**2).sum()
        if D == 0:
            return 0.
        return (self.r**2).sum() / D
    def is_null(self):
        return not np.any(self.r)

class PCGPolarizationMap(PCGMap, _FitsMapOutputter):
    def __init__(self, name=None, fits_map=None, proj_gen=None):
        PCGMap.__init__(self, name=name)
        
        self.map = IQUMap()
        for code in 'IQU':
            self.map.set_image(code, fits_map)
            self.map.data[code].data *= 0
        self.proj_gen = proj_gen
        self.projs = {}
        self.pcg_maps = [PCGIntensityMap(fits_map=self.map.data[code])
                         for code in 'IQU']

    def get_ref(self, code):
        if code in ['x', 'r', 'q', 'z', 'p', 'weight']:
            return [getattr(m, code) for m in self.pcg_maps]

    def set(self, code, value):
        setattr(self, code, value)

    def clear(self, code):
        vects = self.get_ref(code)
        for vect in vects:
            vect[:] = 0

    def accumulate(self, code, vals):
        vects = self.get_ref(code)
        for vect, v in zip(vects, vals):
            vect[:] += v

    def get_qp(self):
        return sum([ravel_dot(m.q,m.p) for m in self.pcg_maps])

    def get_rz(self):
        return sum([ravel_dot(m.r,m.z) for m in self.pcg_maps])

    def add_TOD_step1(self, key, tod, cuts):
        self.tods.append(key)
        self.projs[key] = self.proj_gen.get(
            map_name=self.name, tod_name=key, tod=tod,
            fits_map=self.pcg_maps[0].map)
        # Project an initial solution from the TOD?
        # - This requires two full passes through all TODs, with TOD data.
        #
        # Deproject starting solution from TOD!
        if np.any([np.any(m.x) for m in self.pcg_maps]):
            iqu = [-m.x for m in self.pcg_maps]
            self.projs[key].deproject_from_map(iqu, tod.data, cuts=cuts)
    def add_TOD_step2(self, key, data, cuts):
        # Initialize r from the TOD
        iqu = [m.r for m in self.pcg_maps]
        self.projs[key].project_to_map(data, iqu_output=iqu, cuts=cuts,
                                       weights_output=self.pcg_maps[0].map.weight)
    def prepare_iteration(self):
        [m.prepare_iteration() for m in self.pcg_maps]
    def project_from_p(self, key, data, cuts):
        iqu = [m.p for m in self.pcg_maps]
        self.projs[key].deproject_from_map(iqu, output=data, cuts=cuts)
    def deproject_to_q(self, key, data, cuts):
        iqu = [m.q for m in self.pcg_maps]
        self.projs[key].project_to_map(data, iqu_output=iqu, cuts=cuts)
    def step_x(self, alpha):
        [m.step_x(alpha) for m in self.pcg_maps]
    def step_p(self, beta):
        [m.step_p(beta) for m in self.pcg_maps]
    def get_convergence(self):
        return max([m.get_convergence() for m in self.pcg_maps])
    def is_null(self):
        return all([m.is_null() for m in self.pcg_maps])


class PCGMapContainer(PCGMap):
    """
    Holder for multiple PCGMaps, that acts as a PCGMap.  This can be
    used as a base class for actual useful containers.
    """
    def __init__(self, name=None, child_maps=[], child_map_info=[]):
        PCGMap.__init__(self, name=name)
        self.sub_maps = [m for m in child_maps]
        self.sub_info = [c for c in child_map_info]
        assert len(child_maps) == len(child_map_info)
    def add_TOD_step1(self, key, tod, cuts):
        for im, m in enumerate(self.sub_maps):
            m.add_TOD_step1(key, tod, cuts)
    def add_TOD_step2(self, key, data, cuts):
        for im, m in enumerate(self.sub_maps):
            m.add_TOD_step2(key, data, cuts)
    def prepare_iteration(self):
        for i, m in enumerate(self.sub_maps):
            m.prepare_iteration()
    def project_from_p(self, key, data, cuts):
        for im, m in enumerate(self.sub_maps):
            m.project_from_p(key, data, cuts)
    def deproject_to_q(self, key, data, cuts):
        for im, m in enumerate(self.sub_maps):
            m.deproject_to_q(key, data, cuts)
    def step_x(self, alpha):
        [m.step_x(alpha) for m in self.sub_maps]
    def step_p(self, beta):
        [m.step_p(beta) for m in self.sub_maps]
    def get_convergence(self):
        return max([m.get_convergence() for m in self.sub_maps])
    def get_qp(self):
        return sum([m.get_qp() for m in self.sub_maps])
    def get_rz(self):
        return sum([m.get_rz() for m in self.sub_maps])


class _FitsMapContainer:
    def write(self, filename, outputm=None, info=None, **kwargs):
        assert filename.endswith('.fits')
        info_copy = {}
        if info is not None:
            info_copy.update(info)
        for m, sub_info in zip(self.sub_maps, self.sub_info):
            info_copy.update(sub_info)
            m.write(filename, outputm=outputm, info=info_copy, **kwargs)

    def plot(self, filename, outputm=None, info=None, **kwargs):
        assert filename.endswith('.png')
        info_copy = {}
        if info is not None:
            info_copy.update(info)
        for m, sub_info in zip(self.sub_maps, self.sub_info):
            info_copy.update(sub_info)
            m.plot(filename, outputm=outputm, info=info_copy, **kwargs)

class PCGDetectorMapSplitter(PCGMapContainer, _FitsMapContainer):
    """
    This is a container for mapping the same region with different
    subsets of array detectors.  The pointing wand is re-used and the
    cuts are automatically divided.
    """
    def __init__(self, name=None, fits_map=None, proj_gen=None,
                 split_params=None, map_class=None):
        """
        name, fits_map, proj_gen have their usual meaning.

        split_params is a list of instructions for selecting detectors
        based on tod.info.array_data.  map_class should be either
        PCGIntensityMap or PCGPolarizationMap.
        """
        # split_params will be used to project different detectors
        # into different maps.  Translate the requests into
        # "select_inner" format.
        ## This is a list, where each element represents a subset of
        ## the detectors.  Each element of the list is a tuple of two
        ## dictionaries.  The first will be passed to
        ## array_data.select_inner.  The second will be used to format
        ## strings that describe the subset.

        split_list = [({},{})]
#        self.split_params = [({},{})]
        for sp_par in split_params:
            sp_type, sp_par = sp_par[0], sp_par[1:]
            if sp_type == 'array_data_parameter':
                ad_par, vals = sp_par[0], sp_par[1]
                to_merge = [({ad_par: [v]}, {ad_par: v}) for v in vals]
            elif sp_type == 'det_uid_groups':
                gr_key, val_lists = sp_par[0], sp_par[1]
                to_merge = [({'det_uid': vals},{gr_key: gr_val})
                            for (gr_val, vals) in val_lists]
            else:
                raise ValueError("Unknown split_param spec '%s'" % sp_type)
            new_split_list = []
            for (sel0, inf0) in split_list:
                for (sel1, inf1) in to_merge:
                    sel0.update(sel1)
                    inf0.update(inf1)
                    new_split_list.append((sel0.copy(), inf0.copy()))
            split_list = new_split_list
        # Split that into the list of selection dicts and list of info dicts.
        self.split_params, split_info = list(zip(*split_list))
        # Tack on a basic numbering scheme...
        for spar_i, inf0 in enumerate(split_info):
            inf0['split'] = spar_i

        # With those detector splits... we create multiple maps and a
        # special projection generator object.
        self.projs = {}
        self.proj_gen = proj_gen
        class MiniProjGenlet:
            def get(_self, map_name=None, tod_name=None, tod=None,
                    fits_map=None):
                if not tod_name in self.projs:
                    self.projs[tod_name] = self.proj_gen.get(
                        map_name=self.name, tod_name=tod_name,
                        tod=tod, fits_map=fits_map)
                return self.projs[tod_name]
        mpg = MiniProjGenlet()
        child_maps = [map_class(name, fits_map.copy(), mpg)
                      for i in range(len(self.split_params))]
        PCGMapContainer.__init__(self, name=name, child_maps=child_maps,
                                 child_map_info=split_info)
        self.child_cuts = {}

    def get_cuts(self, cuts, split_idx, key, tod=None):
        tod_cuts = self.child_cuts.get(key, {})
        if not split_idx in tod_cuts:
            # Base on TOD's det_uid
            cv = moby2.tod.TODCuts(nsamps=cuts.nsamps, det_uid=cuts.det_uid,
                                   sample_offset=cuts.sample_offset)
            cut_det = moby2.tod.CutsVector.new_always_cut(cv.nsamps)
            # Throw out dets not listed.
            keepers = list(tod.info.array_data.select_outer(
                self.split_params[split_idx]))
            for i, u in enumerate(cv.det_uid):
                if not u in keepers:
                    cv.cuts[i] = cut_det
            tod_cuts[split_idx] = cv
            self.child_cuts[key] = tod_cuts
        cv = tod_cuts[split_idx].copy()
        cv.merge_tod_cuts(cuts)
        return cv
    def add_TOD_step1(self, key, tod, cuts):
        for im, m in enumerate(self.sub_maps):
            m.add_TOD_step1(key, tod, self.get_cuts(cuts, im, key, tod=tod))
    def add_TOD_step2(self, key, data, cuts):
        for im, m in enumerate(self.sub_maps):
            m.add_TOD_step2(key, data, self.get_cuts(cuts, im, key))
    def prepare_iteration(self):
        for m in self.sub_maps:
            m.prepare_iteration()
    def project_from_p(self, key, data, cuts):
        for im, m in enumerate(self.sub_maps):
            m.project_from_p(key, data, self.get_cuts(cuts, im, key))
    def deproject_to_q(self, key, data, cuts):
        for im, m in enumerate(self.sub_maps):
            m.deproject_to_q(key, data, self.get_cuts(cuts, im, key))

