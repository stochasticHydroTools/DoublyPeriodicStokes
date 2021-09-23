import ctypes
import numpy as np
from Ghost import *

"""
Python wrappers for C library Spread/Interp routines

See SpreadInterpWrapper.cpp

The prototypes for relevant functions from the 
C++ SpreadInterp library are declared. Any functions added
to the "extern" definition in SpreadInterpWrapper.cpp should be
declared here.
"""
libSpreadInterp = ctypes.CDLL('libspreadInterp.so')
libSpreadInterpDerivX = ctypes.CDLL('libspreadInterpDerivX.so')
libSpreadInterpDerivY = ctypes.CDLL('libspreadInterpDerivY.so')
libSpreadInterpDerivZ = ctypes.CDLL('libspreadInterpDerivZ.so')

libSpreadInterp.Spread.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libSpreadInterp.Spread.restype = None
libSpreadInterp.Interpolate.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libSpreadInterp.Interpolate.restype = None
libSpreadInterp.Resetdof.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_uint]
libSpreadInterp.Resetdof.restype = None

libSpreadInterp.addDx.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libSpreadInterp.addDx.restype = None
libSpreadInterp.addDy.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libSpreadInterp.addDy.restype = None
libSpreadInterp.addDz.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libSpreadInterp.addDz.restype = None

libSpreadInterpDerivX.Spread.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libSpreadInterpDerivX.Spread.restype = None
libSpreadInterpDerivX.Interpolate.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libSpreadInterpDerivX.Interpolate.restype = None

libSpreadInterpDerivY.Spread.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libSpreadInterpDerivY.Spread.restype = None
libSpreadInterpDerivY.Interpolate.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libSpreadInterpDerivY.Interpolate.restype = None

libSpreadInterpDerivZ.Spread.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libSpreadInterpDerivZ.Spread.restype = None
libSpreadInterpDerivZ.Interpolate.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libSpreadInterpDerivZ.Interpolate.restype = None

def Spread(s, g):
  """
  Spread data from the particles s onto the grid g.
  
  Parameters:
    s - a valid ParticleList
    g - a valid GridGen
  
  Returns: Pointer to spread data  
  Side Effects:
    The C++ Grid data member g.fG is populated with the spread data
  """
  # initialize the extended grid data 
  g.ZeroExtGrid()
  # spread forces on the particles (C lib)
  libSpreadInterp.Spread(s.particles, g.grid)
  # enforce boundary conditions on spread data
  DeGhostify(g.grid, s.particles)
  return g.GetData()

def SpreadDx(s, g):
  """ 
  same as above but for x-derivative of kernel, and
  we return a copy of the data. this is so that 
  the same grid/particle combo can be re-used immediately
  after executing spread. 

  Returns: A separate copy of the spread data
  
  """
  g.ZeroExtGrid()
  libSpreadInterpDerivX.Spread(s.particles, g.grid)
  DeGhostify(g.grid, s.particles)
  return g.GetData()

def SpreadDy(s, g):
  """ same as above but for y-derivative of kernel """
  g.ZeroExtGrid()
  libSpreadInterpDerivY.Spread(s.particles, g.grid)
  DeGhostify(g.grid, s.particles)
  return g.GetData()

def SpreadDz(s, g):
  """ same as above but for z-derivative of kernel """
  g.ZeroExtGrid()
  libSpreadInterpDerivZ.Spread(s.particles, g.grid)
  DeGhostify(g.grid, s.particles)
  return g.GetData()

def Interpolate(s, g, uP):
  """
  Interpolate data from the grid g onto the particles s.
  
  Parameters:
    s - A valid ParticleList
    g - A valid GridGen
    uP - output array with interpolated data 
  Side Effects:
    The C++ ParticleList data member s.fP is populated with the interpolated data
  """
  # reinitialize forces on particles for interpolation
  s.ZeroData()
  # copy data to ghost cells to enforce boundary conditions before interp
  Ghostify(g.grid, s.particles)
  # interpolate velocities on the particles (C lib)
  libSpreadInterp.Interpolate(s.particles, g.grid)
  s.CopyData(uP)

def InterpolateDx(s, g, dxuP):
  """
  Interpolate data from the grid g onto the particles s using x deriv of kernel.
  
  Parameters:
    s - A valid ParticleList
    g - A valid GridGen
    dxuP - output array with interpolated data 
  
  Returns: copy of interpolated data 
  Side Effects:
    The C++ ParticleList data member s.fP is populated with the interpolated data
  """
  s.ZeroData()
  Ghostify(g.grid, s.particles)
  libSpreadInterpDerivX.Interpolate(s.particles, g.grid)
  s.CopyData(dxuP)

def InterpolateDy(s, g, dyuP):
  """
  Interpolate data from the grid g onto the particles s using x deriv of kernel.
  
  Parameters:
    s - A valid ParticleList
    g - A valid GridGen
    dyuP - output array with interpolated data 
  
  Returns: copy of interpolated data 
  Side Effects:
    The C++ ParticleList data member s.fP is populated with the interpolated data
  """
  s.ZeroData()
  Ghostify(g.grid, s.particles)
  libSpreadInterpDerivY.Interpolate(s.particles, g.grid)
  s.CopyData(dyuP)

def InterpolateDz(s, g, dzuP):
  """
  Interpolate data from the grid g onto the particles s using x deriv of kernel.
  
  Parameters:
    s - A valid ParticleList
    g - A valid GridGen
    dzuP - output array with interpolated data 
  
  Returns: copy of interpolated data 
  Side Effects:
    The C++ ParticleList data member s.fP is populated with the interpolated data
  """
  s.ZeroData()
  Ghostify(g.grid, s.particles)
  libSpreadInterpDerivZ.Interpolate(s.particles, g.grid)
  s.CopyData(dzuP)

def ResetDOF(s, g, dof):
  """
  Reset dofs for particles and grid so we can reuse
  search structures for spread/interp with any dof
  """
  libSpreadInterp.Resetdof(g.grid, s.particles, dof)
