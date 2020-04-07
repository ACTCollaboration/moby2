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

MOBY2_PREFIX HAS BEEN REMOVED.  THESE DOCS NEED UPDATE:

- use an anonymous PREFIX for installing libactpol_deps and libcatpol
- provide standard environment modifications for making those deps
  visible to moby2 build system (env script and modulefile)


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


libactpol dependencies
----------------------

As of mid-2018, all dependencies are either easily obtained through
system package managers or pip, or are stored in an ACTCollaboration
github repository.

libactpol dependencies: easily obtained packages
------------------------------------------------

The patched version of libactpol that we use depends on some less
common dev packages.  You should be able to find the following
packages in your distribution's package manager:

* libcfitsio3-dev -- ``sudo apt-get install libcfitsio3-dev``
* libzzip-dev -- ``sudo apt-get install libzzip-dev``
* wcslib -- ``sudo apt-get install wcslib-dev``

(For Redhat you want something like: ``sudo yum install zziplib-devel
cfitsio-devel wcslib``.)

If you do not have root access on your machine, see if the system
administrator has or can make them available.  Alternately, install
them just for your user account.

libactpol dependencies: special modules
---------------------------------------

There are three rather specialized packages required by libactpol:
libslim (compression), sofa (astrometric conversions), slarefro
(refraction).  These are most easily obtained through the
libactpol_deps repository.

You can access the repository by cloning::

  git clone ssh://git@github.com/ACTCollaboration/libactpol_deps.git

This repository contains 3 installable modules.  **See the README file
for the latest instructions.**

**libslim**: In order to support uint8, we may be using a patched
version of libslim.  This may become unneccessary in the future.

**sofa**: This is the Standards of Fundamental Astronomy library from
the IAU, http://www.iausofa.org/ .  At this writing, we use a recent
libsofa, unaltered except to include a Makefile.  This may change in
the future to support leap seconds more flexibly.

**sla_refro**: This is a very simple Fortran -> C wrapping of a single
function from slalibf that is used by libactpol to compute atmospheric
refraction.

Once all three of these packages have been installed, it should be
possible to compile libactpol.


Build patched libactpol
-----------------------

You can access the repository by cloning::

  git clone ssh://git@github.com/ACTCollaboration/libactpol.git
  cd libactpol.git

As of this writing **moby2 does not work with the** ``master`` **branch of**
``libactpol``!  Instead you should switch to the ``moby2_mods`` branch::

  git checkout moby2_mods

Then proceed with::

  autoreconf -i
  ./configure --enable-shared --disable-oldact --disable-slalib --prefix=$MOBY2_PREFIX
  make
  make install


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

  cp /mnt/act3/users/mhasse/shared/dot_moby2_feynman $HOME/.moby2


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

