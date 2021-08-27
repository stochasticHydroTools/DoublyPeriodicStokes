import ctypes
import numpy as np
from config import num_threads
libTransform = ctypes.CDLL('libtransform.so')
class Transformer(object):
  """
  Python wrappers for C library Transform routines.
  
  This class can be thought of as a generator and
  manager of the underlying C++ Transform struct.
  
  See TransformWrapper.cpp.
  
  Attributes:
    in_real (double array) - real part of the input data.
    in_imag (double array or None) - complex part of the input data.
    Nx, Ny, Nz - number of points in x, y and z.
    N - Nx * Ny * Nz.
    dof - degrees of freedom of the data.
    Ntotal - N * dof.
    out_real (double array) - real part of output transform.
    out_imag (double array) - complex part of output transform.
    tType - type of transform (0 for fourier, 1 for fourier-cheb)
    transform (ptr to C++ struct) - a pointer to the generated C++ Transform struct
  """
  def __init__(self, _Nx, _Ny, _Nz, _dof, _tType):
    """ 
    The constructor for the Transformer class.
    
    Parameters:
      in_real (doubles) - real part of input.
      in_imag (doubles) - complex part of input.
      Nx, Ny, Nz (int) - number of points in x,y,z
      dof (int) - degrees of freedom.

    Side Effects: None

    """ 
    # number of points in x,y,z
    self.Nx = _Nx
    self.Ny = _Ny
    self.Nz = _Nz
    # degrees of freedom
    self.dof = _dof
    # get total nums
    self.N = self.Nx * self.Ny * self.Nz
    self.Ntotal = self.N * self.dof
    # outputs
    self.out_real = None
    self.out_imag = None
    # type of transform (0 for fourier, 1 for fourier-cheb)
    self.tType = _tType
    # pointer to initialized c++ Transform struct
    self.transform = libTransform.InitTransforms(self.Nx, self.Ny, self.Nz, self.dof, num_threads, self.tType)

  def Ftransform(self, in_real):
    """
    Python wrapper for the Ftransform(...) C lib routine.

    This executes a forward transform on the input data, assuming that it is real.

    Parameters: None
    Side Effects:
      self.out_real is populated with the real part of the output transform
      self.out_imag is populated with the complex part of the output transform

    """
    libTransform.SetFData(self.transform, in_real.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    libTransform.Ftransform(self.transform)
    self.out_real = np.ctypeslib.as_array(libTransform.getRealOut(self.transform), shape=(self.Ntotal,))
    self.out_imag = np.ctypeslib.as_array(libTransform.getComplexOut(self.transform), shape=(self.Ntotal,))

  def Btransform(self, in_real, in_imag):
    """
    Python wrapper for the Btransform(...) C lib routine.

    This computes the backward plan and executes a backward
    transform on the input data. The output is normalized 
    by self.N

    Parameters: None
    Side Effects:
      self.transform is assigned the pointer to the C++ Transform instance
      self.out_real is populated with the real part of the output transform
      self.out_imag is populated with the complex part of the output transform
    """
    libTransform.SetBData(self.transform, in_real.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                          in_imag.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    libTransform.Btransform(self.transform)
    self.out_real = np.ctypeslib.as_array(libTransform.getRealOut(self.transform), shape=(self.Ntotal,)) 
    self.out_imag = np.ctypeslib.as_array(libTransform.getComplexOut(self.transform), shape=(self.Ntotal,))

  def Clean(self):
    """
    Python wrapper for the CleanTransform(..) C lib routine.
    This cleans the Transform struct returned by B/Ftransform,
    frees any memory internally allocated, and deletes the
    pointer to the transform struct stored in the class.

    Parameters: None
    Side Effects:
      self.transform is deleted and nullified
      self.out_real is deleted and nullified
      self.out_imag is deleted and nullified
    """
    libTransform.CleanTransform(self.transform)
    libTransform.DeleteTransform(self.transform)

"""
The prototypes for relevant functions from the 
C++ Transform library are declared. Any functions added
to the "extern" definition in Transform.h should be
declared here. The attributes out_real, out_imag 
and transform are set to None.
"""


libTransform.InitTransforms.argtypes = [ctypes.c_uint, ctypes.c_uint,\
                                        ctypes.c_uint, ctypes.c_uint,\
                                        ctypes.c_uint, ctypes.c_int]
libTransform.InitTransforms.restype = ctypes.c_void_p

libTransform.Ftransform.argtypes = [ctypes.c_void_p]
libTransform.Ftransform.restype = None

libTransform.Btransform.argtypes = [ctypes.c_void_p]
libTransform.Btransform.restype = None


libTransform.SetFData.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)]
libTransform.SetFData.restype = None

libTransform.SetBData.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double),\
                                  ctypes.POINTER(ctypes.c_double)]
libTransform.SetBData.restype = None

libTransform.getRealOut.argtypes = [ctypes.c_void_p]
libTransform.getRealOut.restype = ctypes.POINTER(ctypes.c_double)

libTransform.getComplexOut.argtypes = [ctypes.c_void_p]
libTransform.getComplexOut.restype = ctypes.POINTER(ctypes.c_double)

libTransform.CleanTransform.argtypes = [ctypes.c_void_p]
libTransform.CleanTransform.restype = None

libTransform.DeleteTransform.argtypes = [ctypes.c_void_p]
libTransform.DeleteTransform.restype = None
