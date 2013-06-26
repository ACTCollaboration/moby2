from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
"""
Wrappers for scipy.optimize routines that make (relatively) easy to
hold some parameters fixed while letting others vary.

It would be cool to get an MCMC fitter in here somehow.
"""

from scipy.optimize import leastsq, fmin
import numpy as np

class _MFExpander:
    def __init__(self, func, free, x0):
        self.func, self.free, self.x0 = func, np.array(free), np.array(x0)

    def __call__(self, x, *args):
        x_all = self.get_all_x(x)
        return self.func(x_all, *args)
        
    def get_all_x(self, x):
        x_all = self.x0.copy()
        x_all[self.free] = x
        return x_all
        

def multi_leastsq(func, x0, args=(), free=None, Dfun=None, diag=None, **kwargs):
    """
    This is intended as an almost drop-in replacement for
    scipy.optimize.leastsq.  Most arguments (including those buried in
    kwargs) are described in the scipy documentation.  The "free"
    argument should be a list of integers or booleans, indicating
    which parameters in x0 should be varied.  The other parameters
    remain fixed to the values given in x0.

    The function returns the full vector of optimized values,
    including any that were held fixed.
    """
    if free is None:
        free = np.ones(len(x0), 'bool')
    else:
        free = np.asarray(free)

    mf_func = _MFExpander(func, free, x0)
    if diag is not None:
        diag = diag[free]
    if Dfun is not None:
        Dfun = _MFExpander(func, free, x0)

    # This one always returns a tuple.  Thanks.
    output = leastsq(mf_func, np.asarray(x0)[free], args=args, Dfun=Dfun, diag=diag, **kwargs)

    return (mf_func.get_all_x(output[0]),) + output[1:]

def multi_fmin(func, x0, args=(), free=None, **kwargs):
    """
    This is intended as an almost drop-in replacement for
    scipy.optimize.fmin.  Most arguments (including those buried in
    kwargs) are described in the scipy documentation.  The "free"
    argument should be a list of integers or booleans, indicating
    which parameters in x0 should be varied.  The other parameters
    remain fixed to the values given in x0.

    The function returns the full vector of optimized values,
    including any that were held fixed.

    Output of retall=True is not expanded.  It probably should be.
    """
    if free is None:
        free = np.ones(len(x0), 'bool')
    else:
        free = np.asarray(free)
    mf_func = _MFExpander(func, free, x0)

    out_is_tuple = bool(kwargs.get('full_output', False)) or \
        bool(kwargs.get('retall', False))

    # Get output from fmin, which might be a tuple
    output = fmin(mf_func, np.asarray(x0)[free], args, **kwargs)

    if out_is_tuple:
        return (mf_func.get_all_x(output[0]),) +  output[1:]

    return mf_func.get_all_x(output)

def multi_mc(func, x0, args=(), free=None, x0_step=None,
             n_samp=1000, n_burn=300, n_restep=200,
             verbosity=1):
    """
    Run a Markov Chain Monte Carlo using the likelihood function func.
    func will be called with arguments

        func(p, *args)

    where args is optionally provided by the user.  x0 is the starting
    point for the solution.  x0_step is the starting step size, and
    can be a scalar, a vector of same length as x0 (assumed diagonal)
    or 2d (assumed to be covariance matrix).  To only step certain
    parameters, list their indices in free, e.g. free=[0,1,4].
    """

    def debug_print(level, msg, iteration=None):
        if verbosity >= level:
            if iteration is not None:
                print('Iter %5i' % iteration, end=' ')
            print(msg)
    x0 = np.asarray(x0, float)

    # Identify parameters we will manipulate
    if free is None:
        free = np.ones(len(x0), 'bool')
    else:
        free = np.asarray(free)
    if free.dtype == 'bool':
        free = free.nonzero()[0]

    last_sample = x0
    last_like = func(last_sample, *args)

    if x0_step is None:
        debug_print(2, 'Estimating step sizes')
        # Measure sensitivity to parameter changes, I guess.
        # Get diagonal part of covariance matrix
        step_vects = np.zeros((len(free), len(free)))
        for fi, i in enumerate(free):
            x1 = np.array(x0)
            dx = max(x1[i]*1e-3, 1e-3)
            x1[i] = x0[i] + dx
            dL = func(x1, *args) - last_like
            if abs(dL/dx) < 1e-7:
                debug_print(2, 'Parameter %i is insensitive.' % i)
                step_vects[fi,fi] = dx
            else:
                step_vects[fi,fi] = 0.2 * dx/dL
            debug_print(2, 'Parameter %i step size = %e' % (i, step_vects[fi,fi]))
        step_sizes = np.ones(len(free))
        #w, v = np.diag([1]*len(free)), x0_stdev[free]
    else:
        # Up convert input step to a covariance matrix.
        x0_step = np.asarray(x0_step)
        if len(x0_step.shape) == 0:
            x0_step = x0_step + x0*0
        if x0_step.ndim == 1:
            x0_step = np.diag(x0_step)**2
        # Factor 
        step_sizes, step_vects = np.linalg.eig(x0_step[free][:,free])
        #print w, v
        # Kill negatives
        step_sizes[step_sizes<=0] = 0.
        step_sizes = step_sizes**.5 * .4
        for s, v in zip(step_sizes, step_vects):
            debug_print(2, 'step vector: %e %s' % (s, str(v)))


    #x0_stdev = np.asarray(x0_stdev)
    #if len(x0_stdev.shape) == 0:
    #    x0_stdev = x0_stdev + x0*0
        
    samples = np.zeros((n_burn + n_samp, len(x0)))
    likes = np.zeros(n_burn + n_samp)

    acceptance = np.zeros(n_restep, bool)
    acceptance_i = 0

    debug_print(2, 'Initial sample:    ' + str(last_sample))
    debug_print(2, 'Initial -log_like: %.4e' % last_like)

    use_covar = True

    dx = np.zeros(len(x0))
    n_accept, n_reject = 0,0
    n_next_restep = 0

    alert_interval = n_samp // 10

    # The sample loop
    for i in range(n_burn+n_samp):
        # Next sample
        if not use_covar:
            dx[free] = np.random.normal(size=len(free)) * step[free]
        else:
            dx[free] = np.dot(step_vects,
                              step_sizes*np.random.normal(size=len(free)))

        trial_sample = last_sample + dx
        trial_like = func(trial_sample, *args)
        accept = (trial_like < last_like or 
                np.random.uniform() < np.exp(last_like - trial_like))

        debug_print(3, 'Step: ' + str(dx), i)
        debug_print(3, 'Sample: ' + str(trial_sample), i)
        debug_print(3, 'Like: %.4e (%s)' % \
                        (trial_like, {True: 'accepted',
                                      False:'rejected'}[bool(accept)]),i)

        acceptance[n_next_restep] = accept
        if accept:
            n_accept += 1
            last_sample = trial_sample
            last_like = trial_like

        samples[i] = last_sample
        likes[i] = last_like

        n_next_restep += 1
        if n_next_restep >= n_restep:
            debug_print(1, 'Estimating new steps.', i)
            acceptance_f = 1.*acceptance.sum() / n_next_restep
            debug_print(2, 'Acceptance level is %.4f' % acceptance_f, i)
            if acceptance_f < .1:
                debug_print(1, 'rejecting too much; shrinking step size.', i)
                #step[free] /= 4.
                step_sizes /= 3
            elif acceptance_f > .9:
                debug_print(1, 'accepting too much; increasing step size.', i)
                #step[free] *= 2
                step_sizes *= 2
            elif not use_covar:
                new_step = samples[i-n_restep:i].std(axis=0) * 0.4
                new_step[new_step==0] = step[new_step==0] / 10
                step[free] = new_step[free]
            else:
                # Check slope of likelihood...
                like_y = likes[max(0,i-n_restep):i]
                like_x = np.linspace(-.5, .5, len(like_y))
                like_p = np.polyfit(like_x, like_y, 1)
                conv_param = (like_y - np.polyval(like_p, like_x)).std() / like_y.std()
                if conv_param < .7: # seeking!
                    step_mul = 4.
                else:
                    step_mul = .8
                debug_print(2, 'like_slope param is %.4f' % conv_param, i)

                # Diagonalize covariance
                C = np.cov(samples[max(0,i-n_restep):i,free].transpose())
                step_sizes, step_vects = np.linalg.eig(C)
                # Kill negatives
                step_sizes[step_sizes<=0] = 0.
                step_sizes = step_sizes**.5 * step_mul
            debug_print(2, 'New step sizes and vectors:', i)
            for s, v in zip(step_sizes, step_vects):
                debug_print(2, '  %e %s' % (s, str(v)))
            
            n_next_restep = 0
                
        if n_samp % alert_interval == alert_interval-1:
            acc_rate = float(n_accept) / (n_accept+n_reject)
            n_accept, n_reject = 0, 0
            debug_print(1, 'acceptance rate of last %i samples was %.3f' %
                        (alert_interval, acc_rate), i)

    return likes[n_burn:], samples[n_burn:]
    
