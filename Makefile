# To override PYTHON or PREFIX semi-permanently, set them in Makefile.local
ifeq (Makefile.local,$(wildcard Makefile.local))
  include Makefile.local
endif

PYTHON ?= python
PREFIX ?= /usr/local

default: .build

.build:
	$(PYTHON) setup.py build

install:
	$(PYTHON) setup.py install --prefix=$(PREFIX)

install-user:
	$(PYTHON) setup.py install --user

clean:
	-$(PYTHON) setup.py clean
	-rm -r build
