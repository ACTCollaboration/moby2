ifeq ($(wildcard actpol_config.py),)
# $(error actpol_config.py not found in this directory.)
endif	

PYTHON ?= python
PREFIX = $(shell python -c 'from __future__ import print_function; import actpol_config; print(actpol_config.BASE)')

default: .build

.build:
	$(PYTHON) setup.py build

install:
	$(PYTHON) setup.py install --prefix=$(PREFIX)

clean:
	-$(PYTHON) setup.py clean
	-rm -r build
