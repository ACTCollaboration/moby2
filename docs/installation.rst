Installation
============

When moby2 is properly installed and configured, the following should
run without error:

.. code-block:: shell

  python -c "import moby2"

moby2 depends on a bunch of different packages that are needed to
support high-precision pointing and weird file i/o.

Setup scripts
-------------

The steps are:

* Setup environment
* Install libactpol dependencies
* Install libactpol
* Install moby2 dependencies
* Install moby2

Environment variables
---------------------

We will isolate the software from the main system.  Choose a prefix
where the libraries, include files, python packages will live.  I like
$HOME/build.  Set the MOBY2_PREFIX:

.. code-block:: shell

  export MOBY2_PREFIX=$HOME/build

For building, execute that line and the following lines in the shell
where you will be compiling stuff:

.. code-block:: shell

  export LDFLAGS="-L$MOBY2_PREFIX/lib $LDFLAGS"
  export CFLAGS="-I$MOBY2_PREFIX/include $CFLAGS"
  export CPPFLAGS="-I$MOBY2_PREFIX/include $CPPFLAGS"

For running, expose the packages to python by defining MOBY2_PREFIX
and adding the following lines to .bashrc (note these might need a bit
of modification if you're on a 32-bit system or running older python
versions):

.. code-block:: shell

  export PATH=$MOBY2_PREFIX/bin:${PATH}
  export PYTHONPATH=$PYTHONPATH:$MOBY2_PREFIX/lib64/python2.7/site-packages:$MOBY2_PREFIX/lib/python2.7/site-packages
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MOBY2_PREFIX/lib64:$MOBY2_PREFIX/lib


Build libactpol dependencies
----------------------------

The patched version of libactpol that we use depends on the following
dev packages, which may be installed / installable without much effort:

* libcfitsio3-dev -- sudo apt-get install libcfitsio3-dev
* libzzip-dev -- sudo apt-get install libzzip-dev
* wcslib -- sudo apt-get install wcslib-dev
* libslim -- build *after* installing libzzip-dev
* sofa -- custom Makefile
* slarefro -- a single slalib routine; custom package (requires gfortran)

You may have already installed some of these.  For ones you haven't:

* Go here: http://www.astro.princeton.edu/~mhasse/moby2_install/libactpol_deps
* Descend into a recent dated folder.  E.g. 140124.
* Get the libactpol_deps.tar file, or just the zips you need and the
  build.bash script.  Put them all in a directory somewhere where you
  can cook.
* Edit build.bash and comment out the things you don't want to
  install. E.g. to compile everything except cfitsio and zziplib:

.. code-block:: shell

  CFIT= #cfitsio3340.tar.gz
  SLAR=sla_refro-moby2-1.tar.gz
  SOFA=sofa-moby2-1.tar.gz
  WCSL=wcslib-4.17.tar.bz2
  ZZIP= #zziplib-0.13.59.tar.bz2
  SLIM=slim_2.6.5.tgz

* You might also set "UNZIP" at this point, the folder where the
  source code will be configured and build.  Make sure the folder exists.
* Then run `bash build.bash`.  It will try to make and install
  everything.  At the end it tells you whether it thinks it worked or
  not.  Successful output from the above would be:

.. code-block:: shell

  Summary:
    cfit_ok = 
    slar_ok = 1
    sofa_ok = 1
    wcsl_ok = 1
    zzip_ok = 
    slim_ok = 1


Build patched libactpol
-----------------------

Get `libactpol-1.2.0-moby2-4.tar.gz`_.  Unzip it.  Enter the
directory.  Run:

.. _libactpol-1.2.0-moby2-4.tar.gz:  http://www.astro.princeton.edu/~mhasse/moby2_install/libactpol-moby2/libactpol-1.2.0-moby2-4.tar.gz

.. code-block:: shell

  ./configure --enable-shared --prefix=$MOBY2_PREFIX
  make
  make install

If this fails, make sure you've defined all the variables.


Get moby2 dependencies
----------------------

The -dev packages are needed for building; the python stuff is only
needed at run time.  The moby2 dependencies can all be found by the
Ubuntu package manager, or easy_install.

.. code-block:: shell

  sudo apt-get -y install libfftw3-dev \
                          liblapack-dev \
                          libgsl0-dev \
                          python-dev \
                          python-tz \
                          python-numpy \
                          python-matplotlib \
                          python-scipy \
			  python-mysqldb \
                          python-setuptools
  sudo easy_install pyephem pyfits

If you have to compile your own libfftw3, make sure to enable shared
library and float32 support:

.. code-block:: shell

  ./configure --prefix=$MOBY2_PREFIX --enable-shared --with-pic --enable-single

If you have to compile your own pyephem, do it like this:

.. code-block:: shell

  python setup.py build
  python setup.py install --prefix=$MOBY2_PREFIX



Get moby2
---------

Use git to clone the moby2 repository.  Our main copy is a private
repo on github.com:

.. code-block:: shell

  git clone ssh://git@github.com/ACTCollaboration/moby2.git moby2

In a pinch you can get a less current version from this http mirror:

.. code-block:: shell

  git clone http://www.astro.princeton.edu/~mhasse/repos/moby2.git moby2


Compile and install moby2
-------------------------

In the moby2 source directory:

#. Make sure ``MOBY2_PREFIX`` is set properly.
#. Run ``make``.  Pause for laughter.
#. Run ``make install``.
#. Test it: ``python -c 'import moby2'``.
#. Add the necessary paths to your ``.bashrc``, or whatever, so that
   the system can find ``moby2`` next time you log in.  There's a
   template in ``python/data/configs/moby2_env``; you can copy it
   somewhere, update the ``MOBY2_PREFIX`` variable, and source the
   resulting file from your ``.bashrc``.
#. Create a ``~/.moby2`` file for your user.  Copy the template from, e.g.
   ``python/data/configs/dot_moby2_actpol``.


Installation on feynman
-----------------------

**Initialize .moby2**

The template copy of .moby2 points to the locations of TOD data, APEX
weather, IOP parameters, etc.  Before trying to run moby2 on feynman,
initialize your .moby2 file from the template copy:

.. code-block:: shell

  cp /mnt/act2/mhasse/shared/dot_moby2_feynman $HOME/.moby2


**Use system install with python2.6 for cluster nodes**

If you want to run on the cluster nodes, you must use python2.6.  So
ensure the "python" module is not loaded, and initialize your
environment as follows.  Note these lines can be simply added to
.bashrc:

.. code-block:: shell

  module unload python
  source /mnt/act2/mhasse/shared/moby2_env2.6.bash


**Use python2.7 on head node or on node030**

The python2.7 version of the shared code is likely to be slightly more
up-to-date.  It's not available on the cluster nodes, but it is
available on feynman headnode and on node030.  To activate python2.7
and the shared moby2, add to .bashrc:

.. code-block:: shell

  module load python
  source /mnt/act2/mhasse/shared/moby2_env.bash

If you are building your own copy of moby2 and patched libactpol on
feynman, please note the following:


**Building from scratch**

Some system environment variables need to be unset for the builds to
work cleanly:

.. code-block:: shell

  unset FLAGS
  unset U

The python dependencies can be installed through one of (choose your
python version...):

.. code-block:: shell

  # python2.7
  easy_install --prefix=$MOBY2_PREFIX/lib64/python2.7/site-packages/ \
    pyephem pyfits
  # or
  easy_install --prefix=$MOBY2_PREFIX/lib64/python2.6/site-packages/ \
    pyephem pyfits


For database access on the cluster nodes, MySQLdb is needed.
easy_install refuses to install this on the head node, because it is
already installed (though in a place not accessible from the cluster
nodes).  A work-around is to launch the easy_install command (similar
to above but with package "mysql-python") within a PBS job.

