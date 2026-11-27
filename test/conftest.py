# hdf5plugin 1.4.1 (the python 2 pin in setup.py) sets HDF5_PLUGIN_PATH and
# refuses to import once h5py is loaded. pytest loads this file before any
# test module, so import it here first.
try:
    import hdf5plugin  # noqa: F401
except ImportError:
    pass
