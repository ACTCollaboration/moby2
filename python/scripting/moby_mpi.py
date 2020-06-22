from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring

# Set to True or False once import of mpi4py has been attempted.
mpi_ok = None


class MobyMPI:
    def __init__(self, use=True, force=False):
        global mpi_ok
        try:
            from mpi4py import MPI
            mpi_ok = True
        except:
            mpi_ok = False

        use = use and (mpi_ok or force)
        if use:
            if not mpi_ok:
                raise RuntimeError("MPI system not found, but MobyMPI "\
                    "called with force=True.")
            comm = MPI.COMM_WORLD
            size, rank = comm.Get_size(), comm.Get_rank()
        else:
            comm, size, rank = None, 1, 0
        self.comm, self.size, self.rank = comm, size, rank

    def join(self):
        if self.size > 1:
            self.comm.bcast(None, root=0)

    def bcast(self, data, root=0):
        if self.size == 1:
            assert(root==0)
            return data
        return self.comm.bcast(data, root=root)

    def gather(self, data, root=0):
        if self.size == 1:
            assert(root==0)
            return [data]
        return self.comm.gather(data, root=root)

    def get_my_jobs(self, jobs):
        if self.size <= 1:
            return jobs
        i,n,nj = self.rank, self.size, len(jobs)
        return jobs[i*nj//n:(i+1)*nj//n]
