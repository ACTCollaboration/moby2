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
* Point to a working .moby2

Environment variables
---------------------

Often it is necessary or desirable to install ``moby2`` and/or its
dependencies into a special area.  You may define the MOBY2_PREFIX to
help with that.  For example, to install stuff into /usr/local, set:

.. code-block:: shell

  export MOBY2_PREFIX=/usr/local

For building, execute that line and the following lines in the shell
where you will be compiling stuff:

.. code-block:: shell

  export LDFLAGS="-L$MOBY2_PREFIX/lib $LDFLAGS"
  export CFLAGS="-I$MOBY2_PREFIX/include $CFLAGS"
  export CPPFLAGS="-I$MOBY2_PREFIX/include $CPPFLAGS"

Later, Python must be able to find the package, the executable
scripts, and the shared libraries.  If you have chosen a non-standard
installation location, the environment must be modified so Python and
friends can find stuff.  **In that case**, modify the lines below for
your Python version, and make sure they, too are loaded into the shell
(perhaps through your .bashrc, along with the lines MOBY2_PREFIX
definition you have chosen):

.. code-block:: shell

  export PATH=$MOBY2_PREFIX/bin:${PATH}
  export PYTHONPATH=$PYTHONPATH:$MOBY2_PREFIX/lib/python3.6/site-packages/
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MOBY2_PREFIX/lib

You might prefer to set these things up in a modulefile.  That would
look like this::

  #%Module 1.0
  #
  # moby and dependencies
  #

  ## Modify these definitions and leave the other stuff...
  set mroot /path/to/moby2_deps
  set pystr python3.6
  set dot_moby2 /path/to/dot_moby2

  setenv                  MOBY2_PREFIX    $mroot
  setenv                  DOT_MOBY2       $dot_moby2
  prepend-path            PATH            $mroot/bin
  prepend-path            LIBRARY_PATH    $mroot/lib
  prepend-path            LD_LIBRARY_PATH $mroot/lib
  prepend-path            MANPATH         $mroot/share/man
  prepend-path            CPATH           $mroot/include
  prepend-path            FPATH           $mroot/include
  prepend-path            PYTHONPATH      $mroot/lib/$pystr/site-packages


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

The -dev packages are needed for compiling and linking the C
extensions. The Python stuff is only needed at run time.  The moby2
dependencies can all be found by the Ubuntu package manager, or pip,
or easy_install.

The -dev packages::

  sudo apt-get -y install libfftw3-dev \
                          liblapack-dev \
                          libgsl0-dev

Depending on your system and the Python version, you will need these
(Python 3)::

                          python3-dev \
                          python3-tz \
                          python3-numpy \
                          python3-matplotlib \
                          python3-scipy \
			  python3-mysqldb \
                          python3-setuptools

or these (Python 2)::
                          python-dev \
                          python-tz \
                          python-numpy \
                          python-matplotlib \
                          python-scipy \
			  python-mysqldb \
                          python-setuptools

And someone also ``pyephem``; perhaps::

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

Use git to clone the moby2 repository.  Our main copy is a public
repo on github.com:

.. code-block:: shell

  git clone ssh://git@github.com/ACTCollaboration/moby2.git moby2


Compile and install moby2
-------------------------

In the moby2 source directory, you can use ``setup.py`` in the usual
way to build and install moby2.  I.e.::

  python setup.py build

and then **one** of::

  python setup.py install --user
  python setup.py install
  python setup.py install --prefix=/path/to/moby2_stuff/

There is a ``Makefile`` as well, with targets:

* ``build`` (default): builds.
* ``install``: installs with --prefix=${PREFIX}
* ``install-user``: installs with --user

If you want to semi-permanently override either PREFIX or PYTHON (to
select a particular version or installation), you can set those in a
file ``Makefile.local``, and they will be included at the start of the
Makefile.  For example, you might populate Makefile.local with::

  PYTHON = python3
  PREFIX = /home/shared/software/moby2/

Test the installation with: ``python -c 'import moby2'`` or similar.

You will probably need to direct the system to a valid ``.moby2``
file.  Define the environment variable DOT_MOBY2 to point to that
file, or just put something valid in ~/.moby2.
