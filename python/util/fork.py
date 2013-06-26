from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import os
import shutil
import sys
import time

class OutputWrangler:
    """
    This is used in place of an output file when forking to trivially
    parallelize computations.

    For example, if you have 100 jobs to do and you want to fork into
    10 processes, and have the output written to output_file.txt, call:

      fork_out = OutputWrangler('output_file.txt, n_forks=10, n_jobs=100)
    
    To then cause a fork, call

      my_fork, my_jobs = fork_out.fork()
    
    Then process the jobs, writing output with fork_out.write(...).
    The output from each process will be written to a unique file.
    The files will be combined once all forks have completed their
    jobs.
    """

    def __init__(self, filename=None, n_forks=1, n_jobs=1, force=False):
        self.dummy = (filename is None)
        if self.dummy:
            filename = 'mobyfork'

        if n_forks > 1 and sys.flags.interactive:
            print('Interactive mode detected, suppressing fork.')
            n_forks = 1

        if n_jobs < n_forks:
            n_forks = n_jobs

        self.tmp_files = [filename + '_%i.tmp' % i for i in range(n_forks)]
        self.done_files = [filename + '_%i.done' % i for i in range(n_forks)]
        self.filename = filename
        self.children = None

        ok = True
        for f in self.done_files:
            if os.path.exists(f):
                if force:
                    print('Removing %s' % f)
                    try:
                        os.remove(f)
                    except:
                        print('Could not remove %s' % f)
                        ok = False
                else:
                    print('File %s exists, remove to procede.' % f)
                    ok = False
        if not ok:
            raise RuntimeError
    
        self.n_forks = n_forks
        self.n_jobs = n_jobs

    def fork(self, fork_index=None):
        """
        Fork into multiple child processes.

        If you do not want to fork, but want to use the object within
        an existing forked process, pass fork_index=my_fork (where
        my_fork is an index from 0 to n_jobs-1).
        """
        if fork_index is None:
            self.children = []
            for i in range(1, self.n_forks):
                pid = os.fork()
                if pid == 0:
                    print('Spawning child process %i' % i)
                    self._set_index(i)
                    return i, job_indices(self.n_jobs, self.n_forks, i)
                self.children.append(pid)
            fork_index = 0

        self._set_index(fork_index)
        return fork_index, job_indices(self.n_jobs, self.n_forks, fork_index)

    def _set_index(self, index):
        self.index = index
        if len(self.tmp_files) > 0:
            self.fout = open(self.tmp_files[index], 'w')
        
    def write(self, data):
        return self.fout.write(data)

    def flush(self):
        return self.fout.flush()

    def close(self):
        # Close temporary file and rename it
        del self.fout
        os.rename(self.tmp_files[self.index], self.done_files[self.index])
        # If not main process, exit.
        if self.index != 0:
            return
        # Main process waits for all forks to complete.
        first_time = True, 
        not_done = [i for i in self.done_files if not os.path.exists(i)]
        if len(not_done) > 0:
            print('Blocking for all threads to complete...')
            sleepage = 1
            while len(not_done) > 0:
                time.sleep(sleepage)
                sleepage = min(sleepage+1, 10)
                not_done = [i for i in self.done_files if not os.path.exists(i)]

        # Zombie avoidance protocol
        if self.children is not None:
            for pid in self.children:
                os.waitpid(pid,0)

        if not self.dummy:
            print('Assembling output to %s' % self.filename)
            fout = open(self.filename, 'w')
            for infile in self.done_files:
                shutil.copyfileobj(open(infile), fout)
            del fout
        for infile in self.done_files:
            os.remove(infile)

    def __del__(self):
        if hasattr(self, 'fout'):
            self.close()
            
    def cleanup(self):
        for d in self.done_files + self.tmp_files:
            if os.path.exists(d):
                os.remove(d)
        

def job_indices(n_jobs, n_forks, fork_index):
    if n_jobs == 0:
        return []
    return range((fork_index)*n_jobs//n_forks,
                  (fork_index+1)*n_jobs//n_forks)

