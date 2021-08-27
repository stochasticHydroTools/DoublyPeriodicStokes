import ctypes

"""
Python wrappers for C library Ghost point handling routines

See BCWrapper.cpp 

The prototypes for relevant functions from the 
C++ BC library are declared. Any functions added
to the "extern" definition in BCWrapper.cpp should be
declared here.
"""

libBC = ctypes.CDLL('libBC.so')

libBC.Ghostify.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libBC.Ghostify.restype = None
libBC.DeGhostify.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
libBC.DeGhostify.restype = None

def Ghostify(g, s):
  """
  Fold spread data from ghost region of extended grid into the interior
  according to periodicity or boundary condition for each data component

  Parameters:
    g - a pointer to the C++ Grid struct
    s - a pointer to the C++ ParticleList struct
  Side Effects:
    - The member g.fG is modified according to the BC convention and ghost data on g.fG_unwrap
  Returns: none
  """
  libBC.Ghostify(g, s)

def DeGhostify(g,s):
  """
  copy spread data from interior grid to ghost region of extended grid
  according to periodicity or boundary condition for each data component

  Parameters:
    g - a pointer to the C++ Grid struct
    s - a pointer to the C++ ParticleList struct
  Side Effects:
    - The member g.fG_unwrap is modified according to the BC convention and interior data on g.fG
  Returns: none
  """
  libBC.DeGhostify(g, s) 
