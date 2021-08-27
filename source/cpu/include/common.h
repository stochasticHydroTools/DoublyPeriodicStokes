#ifndef _COMMON_H
#define _COMMON_H

/* this file contains enumerations for BC labels and some common functions  */

/* Enumeration for boundary condition types 
   These must be specified at the ends of
   each axis. Note, if periodic is applied
   at the end of one axis, it must be applied
   at the other end as well. */
enum BC {mirror, mirror_inv, none};

// flattened index into 3D array
inline unsigned int at(unsigned int i, unsigned int j,unsigned int k,\
                             const unsigned int Nx, const unsigned int Ny)
{
  return i + Nx * (j + Ny * k);
}


#endif
