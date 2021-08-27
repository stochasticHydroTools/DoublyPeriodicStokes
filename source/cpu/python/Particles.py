import ctypes
import numpy as np
from config import num_threads
libParticles = ctypes.CDLL('libparticles.so')

class ParticlesGen(object):
  """
  Python wrappers for C library ParticlesList routines.
  
  This class can be thought of as a generator and
  manager of the underlying C++ ParticlesList struct.
  
  See ParticlesListWrapper.cpp.
  
  Attributes:
    nP (int) - number of particles, must be specified.
    dof (int) - degrees of freedom, must be specified.
    xP (doubles or None) - particle positions.
    fP (doubles or None) - forces on particles.
    radP (doubles or None) - radii of particles.
    wfP (unsigned shorts or None) - width of kernel for each particle.
    cwfP (doubles or None) - dimensionless radii of particles.
    betafP (doubles or None) - ES kernel beta parameter for each particle.
    particles (ptr to C++ struct) - a pointer to the generated C++ ParticlesList struct
    
    If any of the inputs that can be None are None, the assumption is that
    they will be populated by a call to the C library at a later time
  """
  def __init__(self, _nP, _dof, \
              _xP = None, _fP = None, _radP = None, _wfP = None, _cwfP = None, _betafP = None):
    """ 
    The constructor for the ParticlesGen class.
    
    Parameters:
      nP (int) - number of particles, must be specified.
      dof (int) - degrees of freedom, must be specified.
      xP (doubles or None) - particle positions.
      fP (doubles or None) - forces on particles.
      radP (doubles or None) - radii of particles.
      wfP (unsigned shorts or None) - width of kernel for each particle.
      cwfP (doubles or None) - dimensionless radii of particles.
      betafP (doubles or None) - ES kernel beta parameter for each particle.

    Side Effects: None
    """

    # number of particles
    self.nP = _nP
    # deg of freedom
    self.dof = _dof
    # particle positions
    self.xP = _xP
    # particle forces
    self.fP = _fP
    # beta for ES kernel for each particle (from table)
    self.betafP = _betafP
    # dimensionless radii given ES kernel for each particle (from table)
    self.cwfP = _cwfP
    # width of ES kernel given dimensionless radii (from table)
    self.wfP = _wfP
    # actual radii of the particles
    # - if using the Make_norad() routine, this an array of kernel supports (w * h / 2)
    self.radP = _radP
    # pointer to c++ ParticlesList struct
    self.particles = None
    self.isDipole = False  
    self.useCbeta = False

  def Make(self):
    """
    Python wrapper for the MakeParticles(...) C lib routine.

    This instantiates a ParticlesList object and stores a pointer to it.

    Parameters: None
    Side Effects:
      self.particles is assigned the pointer to the C++ ParticlesList instance
    """
    self.particles = libParticles.MakeParticles(self.xP.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                                          self.fP.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                                          self.radP.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                                          self.betafP.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                                          self.cwfP.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                                          self.wfP.ctypes.data_as(ctypes.POINTER(ctypes.c_ushort)), \
                                          self.nP, self.dof, self.isDipole, self.useCbeta)

    libParticles.SetNumThreads(self.particles, num_threads)


  def Make_norad(self):
    """
    Python wrapper for the MakeParticles(...) C lib routine.

    This instantiates a ParticlesList object and stores a pointer to it.

    NOTE: for this function, cwfp is alphafP - the support of each particle's kernel
          - this should only be used when particle radii are undetermined
          - as is the case if we do not call configure_grid_and_kernels() from Solvers.py
          - primarily for development/effective hydrodynamic radius calculations
    Parameters: None
    Side Effects:
      self.particles is assigned the pointer to the C++ ParticlesList instance
    """
    self.particles = libParticles.MakeParticles_norad(self.xP.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                                          self.fP.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                                          self.betafP.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                                          self.radP.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                                          self.wfP.ctypes.data_as(ctypes.POINTER(ctypes.c_ushort)), \
                                          self.nP, self.dof, self.isDipole)

    libParticles.SetNumThreads(self.particles, num_threads)

  def SetData(self, _fP):
    """
    The python wrapper for setting forces/other data on the particles. This
    modifies the data pointed to by self.particles.

    Parameters: _fP (doubles) - dof x nP fortran ordered array of data
    Side Effects:
      self.particles.fP is created or overwritten with data in _fP 
    """
    libParticles.SetData(self.particles, _fP.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), self.dof)   


  
  def ZeroData(self):
    """
    The python wrapper for zeroing forces/other data on the particles. This
    modifies the data pointed to by self.particles, and should be called
    before interpolation.

    Parameters: None
    Side Effects:
      self.particles.fP is overwritten with 0s, and the program exits if fP is null 
    """
    libParticles.ZeroData(self.particles)   
  
  def GetData(self):
    """
    The python wrapper for getting forces/other data on the particles.

    Parameters: none
    Side Effects: none
    Returns:
      self.particles.fP is returned and encapsulated in numpy array
    """
    return np.ctypeslib.as_array(libParticles.GetData(self.particles), shape=(self.dof * self.nP, )) 
  
  def CopyData(self, data):
    """ same as above, but deep copies the particle data """
    libParticles.CopyData(self.particles, data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

  def Setup(self, g):
    """
    The python wrapper for the Setup(particles,grid) C lib routine
    This function computes internal data structures, like
    the particles-grid locator, effective kernel widths, etc.

    Parameters:
      g - a pointer to a valid GridGen object
    Side Effects:
      The data pointed to by self.particles is modified and extended with
      additional information given the grid.
      The data pointed to by grid is modified and extended with
      additional information given the particles.
    """
    libParticles.Setup(self.particles, g.grid)

  def IsDipole(self):
    """
    Indicate that this instance is of dipoles,
    so derivative of the kernel (and different truncation thresholds)
    will be used for computing unwrapped coordinates/grid locators
    """
    self.isDipole = True
  
  def UseCbeta(self):
    """
    Indicate that this instance uses cwfP as c(beta) = R_h/(wh)
    as opposed to c(w) = R_h/h
    """
    self.useCbeta = True

  def Update(self, xP_new, g):
    """
    Python wrapper for updating particles on ghe grid
  
    Parameters: 
      g - valid GridGen object 
      xP_new - array of new particle positions (must be same size as old)
    Side Effects:
      The data pointed to by self.particles is modified with the new
      particle positions, and the firstn,nextn,number arrays (for particle lookup)
      contained in grid are updated.
    """
    libParticles.Update(self.particles, g.grid, xP_new.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

  def WriteParticles(self, fname):
    """
    Python wrapper for the WriteParticles(particles,fname) C lib routine
    This writes the current state of the ParticlesList to file  
 
    Parameters: 
      fname (string) - desired name of file
    Side Effects: None, besides file creation and write 
    """
    b_fname = fname.encode('utf-8')
    libParticles.WriteParticles(self.particles, b_fname)    
  
  def Clean(self):
    """
    Python wrapper for the CleanParticles(..) C lib routine.
    This cleans the ParticlesList struct returned by Make(),
    frees any memory internally allocated, and deletes the
    pointer to the ParticlesList struct stored in the class.

    Parameters: None
    Side Effects:
      self.particles is deleted (along with underlying data) and nullified
    """
    if self.particles is not None:
      libParticles.CleanParticles(self.particles)
      libParticles.DeleteParticles(self.particles)

"""
  The prototypes for relevant functions from the 
  C++ ParticlesList library are declared. Any functions added
  to the "extern" definition in ParticlesList.h should be
  declared here.
""" 
libParticles.MakeParticles.argtypes = [ctypes.POINTER(ctypes.c_double), \
                                   ctypes.POINTER(ctypes.c_double), \
                                   ctypes.POINTER(ctypes.c_double), \
                                   ctypes.POINTER(ctypes.c_double), \
                                   ctypes.POINTER(ctypes.c_double), \
                                   ctypes.POINTER(ctypes.c_ushort), \
                                   ctypes.c_uint, ctypes.c_uint,\
                                   ctypes.c_bool, ctypes.c_bool]
libParticles.MakeParticles.restype = ctypes.c_void_p 

libParticles.MakeParticles_norad.argtypes = [ctypes.POINTER(ctypes.c_double), \
                                   ctypes.POINTER(ctypes.c_double), \
                                   ctypes.POINTER(ctypes.c_double), \
                                   ctypes.POINTER(ctypes.c_double), \
                                   ctypes.POINTER(ctypes.c_ushort), \
                                   ctypes.c_uint, ctypes.c_uint, ctypes.c_bool]
libParticles.MakeParticles_norad.restype = ctypes.c_void_p 

libParticles.Setup.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libParticles.Setup.restype = None  

libParticles.SetData.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double),\
                                   ctypes.c_uint]
libParticles.SetData.restype = None

libParticles.ZeroData.argtypes = [ctypes.c_void_p]
libParticles.ZeroData.restype = None

libParticles.GetData.argtypes = [ctypes.c_void_p]
libParticles.GetData.restype = ctypes.POINTER(ctypes.c_double)   

libParticles.Update.argtypes = [ctypes.c_void_p, ctypes.c_void_p,\
                                         ctypes.POINTER(ctypes.c_double)]
libParticles.Update.restype = None

libParticles.CleanParticles.argtypes = [ctypes.c_void_p]
libParticles.CleanParticles.restype = None

libParticles.DeleteParticles.argtypes = [ctypes.c_void_p] 
libParticles.DeleteParticles.restype = None

libParticles.WriteParticles.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
libParticles.WriteParticles.restype = None

libParticles.SetNumThreads.argtypes = [ctypes.c_void_p, ctypes.c_int]
libParticles.SetNumThreads.restype = None

libParticles.CopyData.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)]
libParticles.CopyData.restype = None
