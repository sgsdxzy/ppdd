from distutils.core import setup
import py2exe
import numpy
import scipy
import os

setup(windows=[{"script":"ppdd_gui.py", "icon_resources": [(1, "ppdd.ico")]}],
      zipfile=os.path.join('pylibs', 'library.zip'),
      options={"py2exe":{"optimize": 2,
                         "includes":['sip',
                                     'cunwrap._cunwrap',
                                     'scipy.special._ufuncs_cxx',
                                     'scipy.linalg.cython_blas',
                                     'scipy.linalg.cython_lapack',
                                     'scipy.sparse.csgraph._validation']}})
