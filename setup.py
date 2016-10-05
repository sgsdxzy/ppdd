from distutils.core import setup
import py2exe
import numpy
import scipy

setup(windows=[{"script":"ppdd_gui.py"}],
      options={"py2exe":{"includes":['sip',
                                     'cunwrap._cunwrap',
                                     'scipy.special._ufuncs_cxx',
                                     'scipy.linalg.cython_blas',
                                     'scipy.linalg.cython_lapack',
                                     'scipy.sparse.csgraph._validation']}})
