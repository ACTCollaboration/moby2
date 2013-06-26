# analysis/ is where independent (or weakly interdependent) analysis
# packages are stored.  To protect modularity and speed loading time,
# do not import any analysis modules by default!
#
# Practically speaking this means:
# * All modules should be stored in a directory beneath analysis (not
#   in analysis/module_name.py
# * analysis/__init__.py should not import anything.
