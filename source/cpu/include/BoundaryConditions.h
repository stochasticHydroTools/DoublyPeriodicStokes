#ifndef _BOUNDARY_CONDITION_H
#define _BOUNDARY_CONDITION_H
#include<omp.h>
#include"common.h"
#include<iostream>
/* this file contains fold and copy operations for ghost data */

/* implements fold operation to de-ghostify spread data according to BCs
   Specifically, this function
    - copies data in the ghost region of the extended grid Fe to correct region of 
      the interior grid Fe_wrap
    - if periodic[i] = true for axis i = 0,1 or 2, a periodic folding convention
      where data adjacent to a ghost region are copied to the interior region 
      adjacent to the ghost region at the other end of the axis.
    - if periodic[i] = false, then the fold op copies according to convention 
      dispatched by boundary condition on a given data component (dof). data in 
      the ghost region will be copied to the adjacent interior region. 
    - supported conventions are mirror, mirror_inv or none. mirror will copy
      the data with no modification. mirror_inv will copy the negative of the
      ghost data, and none will perform no copy
*/
inline void fold(double* Fe, double* Fe_wrap, const unsigned int Nx_wrap, 
                 const unsigned int Ny_wrap, const unsigned int Nz_wrap, 
                 const unsigned int Nx, const unsigned int Ny, const unsigned int Nz, 
                 const unsigned int dof, bool* periodic, const BC* BCs)
{
  //unsigned int lend = wx, Nx_wrap = Nx - 2 * wx;
  //unsigned int lend = (wx + (wx % 2)) / 2, Nx_wrap = Nx - wx - (wx % 2);
  //unsigned int lend = (wx + (wx % 2)) / 2, Nx_wrap = Nx - (wx - 2 + (wx % 2));
  unsigned int lend = (Nx - Nx_wrap) / 2;
  unsigned int rbeg = Nx - lend;
  //unsigned int bend = wy, Ny_wrap = Ny - 2 * wy;
  //unsigned int bend = (wy + (wy % 2)) / 2, Ny_wrap = Ny - wy - (wy % 2);
  //unsigned int bend = (wy + (wy % 2)) / 2, Ny_wrap = Ny - (wy - 2 + (wy % 2));
  unsigned int bend = (Ny - Ny_wrap) / 2;
  unsigned int tbeg = Ny - bend;
  //unsigned int dend = ext_up, Nz_wrap = Nz - ext_up - ext_down; 
  //unsigned int ubeg = Nz - ext_down;
  unsigned int dend = (Nz - Nz_wrap) / 2;
  unsigned int ubeg = Nz - dend;


  // periodic fold in x 
  if (periodic[0])
  {
    #pragma omp parallel
    {
      // fold eulerian data in y-z plane in ghost region to periodic index
      #pragma omp for collapse(3)
      for (unsigned int k = 0; k < Nz; ++k)
      {
        for (unsigned int j = 0; j < Ny; ++j)
        {
          // first do left
          for (unsigned int i = 0; i < lend; ++i)
          {
            unsigned int ipb = i + Nx_wrap;
            #pragma omp simd 
            for (unsigned int d = 0; d < dof; ++d)
            { 
              Fe[d + dof * at(ipb, j, k, Nx, Ny)] += Fe[d + dof * at(i, j, k, Nx, Ny)]; 
            }
          }
        }
      }
      #pragma omp for collapse(3)
      for (unsigned int k = 0; k < Nz; ++k)
      {
        for (unsigned int j = 0; j < Ny; ++j)
        {
          // now do right
          for (unsigned int i = rbeg; i < Nx; ++i)
          {
            unsigned int ipb = i - Nx_wrap;
            #pragma omp simd 
            for (unsigned int d = 0; d < dof; ++d)
            { 
              Fe[d + dof * at(ipb, j, k, Nx, Ny)] += Fe[d + dof * at(i, j, k, Nx, Ny)]; 
            }
          }  
        }
      }
    }
  }
  // handle BC for each end of x as specified
  else
  {
    // get bc for left end of x 
    const BC* bc_xl = &(BCs[0]);
    // multiplier to enforce bc
    double s;
    for (unsigned int d = 0; d < dof; ++d)
    {
      // we only do something if bc is not none
      if (bc_xl[d] != none)
      { 

        if (bc_xl[d] == mirror) s = 1.0;
        if (bc_xl[d] == mirror_inv) s = -1.0;
        // fold eulerian data in y-z plane in ghost region to adjacent interior region
        // at left end of x axis according to mirror or mirror_inv
        #pragma omp parallel for collapse(3)
        for (unsigned int k = 0; k < Nz; ++k)
        {
          for (unsigned int j = 0; j < Ny; ++j)
          {
            for (unsigned int i = 0; i <= lend; ++i)
            {
              unsigned int ipb = 2 * lend - i;
              Fe[d + dof * at(ipb, j, k, Nx, Ny)] += s * Fe[d + dof * at(i, j, k, Nx, Ny)]; 
            }
          }
        }
      }
    }
    // get bc for right end of x
    const BC* bc_xr = &(BCs[dof]); 
    // no slip wall on x right
    for (unsigned int d = 0; d < dof; ++d)
    {
      // we only do something if bc is not none
      if (bc_xr[d] != none)
      { 
        if (bc_xr[d] == mirror) s = 1.0;
        if (bc_xr[d] == mirror_inv) s = -1.0;
        // fold eulerian data in y-z plane in ghost region to adjacent interior region
        // at right end of x axis according to mirror or mirror_inv
        #pragma omp parallel for collapse(3)
        for (unsigned int k = 0; k < Nz; ++k)
        {
          for (unsigned int j = 0; j < Ny; ++j)
          {
            for (unsigned int i = rbeg - 1; i < Nx; ++i)
            {
              unsigned int ipb = 2 * rbeg - i - 2;
              Fe[d + dof * at(ipb, j, k, Nx, Ny)] += s * Fe[d + dof * at(i, j, k, Nx, Ny)]; 
            }
          }
        }
      }
    }
  }
  // periodic fold in y
  if (periodic[1])
  {
    #pragma omp parallel
    {
      // fold eulerian data in x-z plane in ghost region to periodic index
      #pragma omp for collapse(3)
      for (unsigned int k = 0; k < Nz; ++k)
      {
        // first do bottom
        for (unsigned int j = 0; j < bend; ++j)
        {
          for (unsigned int i = 0; i < Nx; ++i)
          {
            unsigned int jpb = j + Ny_wrap;
            #pragma omp simd 
            for (unsigned int d = 0; d < dof; ++d)
            { 
              Fe[d + dof * at(i, jpb, k, Nx, Ny)] += Fe[d + dof * at(i, j, k, Nx, Ny)]; 
            }
          }
        } 
      } 
      #pragma omp for collapse(3)
      for (unsigned int k = 0; k < Nz; ++k)
      {
        // now do top
        for (unsigned int j = tbeg; j < Ny; ++j)
        {
          for (unsigned int i = 0; i < Nx; ++i)
          {
            unsigned int jpb = j - Ny_wrap;
            #pragma omp simd
            for (unsigned int d = 0; d < dof; ++d)
            { 
              Fe[d + dof * at(i, jpb, k, Nx, Ny)] += Fe[d + dof * at(i, j, k, Nx, Ny)]; 
            }
          }
        }
      }
    }
  }
  // handle BC for each end of y as specified
  else
  {
    // get bc for left end of y 
    const BC* bc_yl = &(BCs[2 * dof]);
    // multiplier to enforce bc
    double s;
    for (unsigned int d = 0; d < dof; ++d)
    {
      // we only do something if bc is not none
      if (bc_yl[d] != none)
      { 
        if (bc_yl[d] == mirror) s = 1.0;
        if (bc_yl[d] == mirror_inv) s = -1.0;
        // fold eulerian data in x-z plane in ghost region to adjacent interior region
        // at bottom end of y axis according to mirror or mirror_inv
        #pragma omp parallel for collapse(3)
        for (unsigned int k = 0; k < Nz; ++k)
        {
          for (unsigned int j = 0; j <= bend; ++j)
          {
            for (unsigned int i = 0; i < Nx; ++i)
            {
              unsigned int jpb = 2 * bend - j;
              Fe[d + dof * at(i, jpb, k, Nx, Ny)] += s * Fe[d + dof * at(i, j, k, Nx, Ny)]; 
            }
          }
        }
      }
    }
    // get bc for right end of y 
    const BC* bc_yr = &(BCs[3 * dof]);
    for (unsigned int d = 0; d < dof; ++d)
    {
      // we only do something if bc is not none
      if (bc_yr[d] != none)
      { 
        if (bc_yr[d] == mirror) s = 1.0;
        if (bc_yr[d] == mirror_inv) s = -1.0;
        // fold eulerian data in x-z plane in ghost region to adjacent interior region
        // at top end of y axis according to mirror or mirror_inv
        #pragma omp parallel for collapse(3)
        for (unsigned int k = 0; k < Nz; ++k)
        {
          for (unsigned int j = tbeg - 1; j < Ny; ++j)
          {
            for (unsigned int i = 0; i < Nx; ++i)
            {
              unsigned int jpb = 2 * tbeg - i - 2;
              Fe[d + dof * at(i, jpb, k, Nx, Ny)] += s * Fe[d + dof * at(i, j, k, Nx, Ny)]; 
            }
          }
        }
      }
    }
  }
  // periodic fold in z
  if (periodic[2])
  {
    #pragma omp parallel
    {
      // fold eulerian data in x-y plane in ghost region to periodic index
      // first do down
      #pragma omp for collapse(3)
      for (unsigned int k = 0; k < dend; ++k)
      {
        for (unsigned int j = 0; j < Ny; ++j)
        {
          for (unsigned int i = 0; i < Nx; ++i)
          {
            unsigned int kpb = k + Nz_wrap;
            #pragma omp simd 
            for (unsigned int d = 0; d < dof; ++d)
            { 
              Fe[d + dof * at(i, j, kpb, Nx, Ny)] += Fe[d + dof * at(i, j, k, Nx, Ny)]; 
            }
          }
        } 
      } 
  
      // now do up
      #pragma omp for collapse(3)
      for (unsigned int k = ubeg; k < Nz; ++k)
      {
        for (unsigned int j = 0; j < Ny; ++j)
        {
          for (unsigned int i = 0; i < Nx; ++i)
          {
            unsigned int kpb = k - Nz_wrap;
            #pragma omp simd 
            for (unsigned int d = 0; d < dof; ++d)
            { 
              Fe[d + dof * at(i, j, kpb, Nx, Ny)] += Fe[d + dof * at(i, j, k, Nx, Ny)]; 
            }
          }
        }
      }
    }
  }
  // handle BC for each end of z as specified
  else
  {
    // get bc for left end of z 
    const BC* bc_zl = &(BCs[4 * dof]);
    // multiplier to enforce bc
    double s;
    for (unsigned int d = 0; d < dof; ++d)
    {
      // we only do something if bc is not none
      if (bc_zl[d] != none)
      { 
        if (bc_zl[d] == mirror) s = 1.0;
        if (bc_zl[d] == mirror_inv) s = -1.0;

        // fold eulerian data in x-y plane in ghost region to adjacent interior region
        // at lower end of z axis according to mirror or mirror_inv
        #pragma omp parallel for collapse(3)
        for (unsigned int k = 0; k <= dend; ++k)
        {
          for (unsigned int j = 0; j < Ny; ++j)
          {
            for (unsigned int i = 0; i < Nx; ++i)
            {
              unsigned int kpb = 2 * dend - k;
              Fe[d + dof * at(i, j, kpb, Nx, Ny)] += s * Fe[d + dof * at(i, j, k, Nx, Ny)]; 
            }
          } 
        }
      }
    }
    // get bc for right end of z
    const BC* bc_zr = &(BCs[5 * dof]);  
    for (unsigned int d = 0; d < dof; ++d)
    {
      // we only do something if bc is not none
      if (bc_zr[d] != none)
      { 
        if (bc_zr[d] == mirror) s = 1.0;
        if (bc_zr[d] == mirror_inv) s = -1.0;

        // fold eulerian data in x-y plane in ghost region to adjacent interior region
        // at upper end of z axis according to mirror or mirror_inv
        #pragma omp parallel for collapse(3)
        for (unsigned int k = ubeg - 1; k < Nz; ++k)
        {
          for (unsigned int j = 0; j < Ny; ++j)
          {
            for (unsigned int i = 0; i < Nx; ++i)
            {
              unsigned int kpb = 2 * ubeg - k - 2;
              Fe[d + dof * at(i, j, kpb, Nx, Ny)] += s * Fe[d + dof * at(i, j, k, Nx, Ny)]; 
            }
          }
        }
      }
    }
  }
  // copy data on extended grid to wrapped grid
  // NOTE: the layout of the wrapped grid is a (2,0,1) permutation of extended
  #pragma omp parallel for collapse(3)
  for (unsigned int j = bend; j < tbeg; ++j)
  {
    for (unsigned int i = lend; i < rbeg; ++i)
    {
      for (unsigned int k = dend; k < ubeg; ++k) 
      {
        unsigned int ii = i - lend, jj = j - bend, kk = k - dend;
        #pragma omp simd 
        for (unsigned int d = 0; d < dof; ++d)
        { 
          Fe_wrap[d + dof * at(kk, ii, jj, Nz_wrap, Nx_wrap)] 
            = Fe[d + dof * at(i, j, k, Nx, Ny)];
        }
      }
    }
  }
}

/* implements copy opertion to enforce periodicity of eulerian data before interpolation
   Specifically, this function
    - copies data on the interior grid Fe_wrap to the ghost region of the
      extended grid Fe
    - if periodic[i] = true for axis i = 0,1 or 2, a periodic folding convention
      where interior data adjacent to a ghost region are copied to the ghost 
      region at the other end of the axis.
    - if periodic[i] = false, then the fold op copies according to convention 
      dispatched by boundary condition on a given data component (dof) 
    - supported conventions are mirror, mirror_inv or none. mirror will copy
      the adjacent data to the ghost region. mirror_inv will copy the negative of the
      adjacent data to the ghost region, and none will perform no copy
*/

inline void copy(double* Fe, const double* Fe_wrap, const unsigned int Nx_wrap, 
                 const unsigned int Ny_wrap, const unsigned int Nz_wrap,
                 const unsigned int Nx, const unsigned int Ny, const unsigned int Nz, 
                 const unsigned int dof, bool* periodic, const BC* BCs)
{
  unsigned int lend = (Nx - Nx_wrap) / 2;
  unsigned int rbeg = Nx - lend;
  //unsigned int bend = wy, Ny_wrap = Ny - 2 * wy;
  //unsigned int bend = (wy + (wy % 2)) / 2, Ny_wrap = Ny - wy - (wy % 2);
  //unsigned int bend = (wy + (wy % 2)) / 2, Ny_wrap = Ny - (wy - 2 + (wy % 2));
  unsigned int bend = (Ny - Ny_wrap) / 2;
  unsigned int tbeg = Ny - bend;
  //unsigned int dend = ext_up, Nz_wrap = Nz - ext_up - ext_down; 
  //unsigned int ubeg = Nz - ext_down;
  unsigned int dend = (Nz - Nz_wrap) / 2;
  unsigned int ubeg = Nz - dend;

  // copy data on wrapped grid to extended periodic grid
  // NOTE: the layout of the wrapped grid is a (2,0,1) permutation of extended
  #pragma omp parallel for collapse(3)
  for (unsigned int k = 0; k < Nz_wrap; ++k)
  {
    for (unsigned int j = 0; j < Ny_wrap; ++j)
    {
      for (unsigned int i = 0; i < Nx_wrap; ++i)
      {
        unsigned int ii = i + lend, jj = j + bend, kk = k + dend;
        for (unsigned int d = 0; d < dof; ++d)
        {
          Fe[d + dof * at(ii, jj, kk, Nx, Ny)] = Fe_wrap[d + dof * at(k, i, j, Nz_wrap, Nx_wrap)];
        }
      }
    }
  }

  // periodic copy in x
  if (periodic[0])
  {
    #pragma omp parallel
    {  
      // copy eulerian data in y-z plane in periodic region to ghost
      #pragma omp for collapse(3)
      for (unsigned int k = dend; k < ubeg; ++k)
      {
        for (unsigned int j = bend; j < tbeg; ++j)
        {
          // first copy right to left
          for (unsigned int i = 0; i < lend; ++i)
          {
            unsigned int ipb = i + Nx_wrap;
            #pragma omp simd 
            for (unsigned int d = 0; d < dof; ++d)
            {
              Fe[d + dof * at(i, j, k, Nx, Ny)] = Fe[d + dof * at(ipb, j, k, Nx, Ny)]; 
            }
          }  
        }
      }
      #pragma omp for collapse(3)
      for (unsigned int k = dend; k < ubeg; ++k)
      {
        for (unsigned int j = bend; j < tbeg; ++j)
        {
          // now copy left to right
          for (unsigned int i = rbeg; i < Nx; ++i)
          {
            unsigned int ipb = i - Nx_wrap;
            #pragma omp simd 
            for (unsigned int d = 0 ; d < dof; ++d)
            {
              Fe[d + dof * at(i, j, k, Nx, Ny)] = Fe[d + dof * at(ipb, j, k, Nx, Ny)]; 
            }
          }
        }
      }
    }
  }
  else
  {
    // get bc for left end of x 
    const BC* bc_xl = &(BCs[0]);
    // multiplier to enforce bc
    double s;
    for (unsigned int d = 0; d < dof; ++d)
    {
      // we only do something if bc is not none
      if (bc_xl[d] != none)
      { 

        if (bc_xl[d] == mirror) s = 1.0;
        if (bc_xl[d] == mirror_inv) s = -1.0;
        // copy eulerian data in y-z plane in interior region to 
        // adjacent ghost region at left end of x axis according to mirror or mirror_inv
        #pragma omp parallel for collapse(3)
        for (unsigned int k = dend; k < ubeg; ++k)
        {
          for (unsigned int j = bend; j < tbeg; ++j)
          {
            for (unsigned int i = 0; i < lend; ++i)
            {
              unsigned int ipb = 2 * lend - i;
              Fe[d + dof * at(i, j, k, Nx, Ny)] = s * Fe[d + dof * at(ipb, j, k, Nx, Ny)]; 
            }
          }
        }
      }
    }
    // get bc for right end of x
    const BC* bc_xr = &(BCs[dof]); 
    // no slip wall on x right
    for (unsigned int d = 0; d < dof; ++d)
    {
      // we only do something if bc is not none
      if (bc_xr[d] != none)
      { 
        if (bc_xr[d] == mirror) s = 1.0;
        if (bc_xr[d] == mirror_inv) s = -1.0;
        // copy eulerian data in y-z plane in interior region to 
        // adjacent ghost region at right end of x axis according to mirror or mirror_inv
        #pragma omp parallel for collapse(3)
        for (unsigned int k = dend; k < ubeg; ++k)
        {
          for (unsigned int j = bend; j < tbeg; ++j)
          {
            for (unsigned int i = rbeg; i < Nx; ++i)
            {
              unsigned int ipb = Nx_wrap - 2 - i + rbeg;
              Fe[d + dof * at(i, j, k, Nx, Ny)] = s * Fe[d + dof * at(ipb, j, k, Nx, Ny)]; 
            }
          }
        }
      }
    }
  }
  // periodic copy in y
  if (periodic[1])
  {
    #pragma omp parallel
    { 
      // copy eulerian data in x-z plane in periodic region to ghost
      #pragma omp for collapse(3)
      for (unsigned int k = dend; k < ubeg; ++k)
      {
        for (unsigned int j = tbeg; j < Ny; ++j)
        {
          // first copy bottom to top
          for (unsigned int i = 0; i < Nx; ++i)
          {
            unsigned int jpb = j - Ny_wrap;
            #pragma omp simd 
            for (unsigned int d = 0; d < dof; ++d)
            {
              Fe[d + dof * at(i, j, k, Nx, Ny)] = Fe[d + dof * at(i, jpb, k, Nx, Ny)]; 
            }
          }
        }
      }
      #pragma omp for collapse(3)
      for (unsigned int k = dend; k < ubeg; ++k)
      {
        for (unsigned int j = 0; j < bend; ++j)
        {
          // now copy top to bottom
          for (unsigned int i = 0; i < Nx; ++i)
          {
            unsigned int jpb = j + Ny_wrap;
            #pragma omp simd 
            for (unsigned int d = 0; d < dof; ++d)
            {
              Fe[d + dof * at(i, j, k, Nx, Ny)] = Fe[d + dof * at(i, jpb, k, Nx, Ny)]; 
            }
          }  
        }
      }
    }
  }
  else
  {
    // get bc for left end of y 
    const BC* bc_yl = &(BCs[2 * dof]);
    // multiplier to enforce bc
    double s;
    for (unsigned int d = 0; d < dof; ++d)
    {
      // we only do something if bc is not none
      if (bc_yl[d] != none)
      { 
        if (bc_yl[d] == mirror) s = 1.0;
        if (bc_yl[d] == mirror_inv) s = -1.0;
        // copy eulerian data in x-z plane in interior region to 
        // adjacent ghost region at bottom end of y axis according to mirror or mirror_inv
        #pragma omp parallel for collapse(3)
        for (unsigned int k = dend; k < ubeg; ++k)
        {
          for (unsigned int j = 0; j < bend; ++j)
          {
            for (unsigned int i = 0; i < Nx; ++i)
            {
              unsigned int jpb = 2 * bend - j;
              Fe[d + dof * at(i, j, k, Nx, Ny)] = s * Fe[d + dof * at(i, jpb, k, Nx, Ny)]; 
            }
          }
        }
      }
    }
    // get bc for right end of y 
    const BC* bc_yr = &(BCs[3 * dof]);
    for (unsigned int d = 0; d < dof; ++d)
    {
      // we only do something if bc is not none
      if (bc_yr[d] != none)
      { 
        if (bc_yr[d] == mirror) s = 1.0;
        if (bc_yr[d] == mirror_inv) s = -1.0;
        // copy of eulerian data in x-z plane in interior region to 
        // adjacent ghost region at top end of y axis according to mirror or mirror_inv
        #pragma omp parallel for collapse(3)
        for (unsigned int k = dend; k < ubeg; ++k)
        {
          for (unsigned int j = tbeg; j < Ny; ++j)
          {
            for (unsigned int i = 0; i < Nx; ++i)
            {
              unsigned int jpb = Ny_wrap - 2 - j + tbeg;
              Fe[d + dof * at(i, j, k, Nx, Ny)] = s * Fe[d + dof * at(i, jpb, k, Nx, Ny)]; 
            }
          }
        }
      }
    }
  }
  // periodic copy in z
  if (periodic[2])
  {
    #pragma omp parallel
    {
      // copy eulerian data in x-y plane in periodic region to ghost
      #pragma omp for collapse(3)
      for (unsigned int k = ubeg; k < Nz; ++k)
      {
        for (unsigned int j = 0; j < Ny; ++j)
        {
          // first copy down to up
          for (unsigned int i = 0; i < Nx; ++i)
          {
            unsigned int kpb = k - Nz_wrap;
            #pragma omp simd
            for (unsigned int d = 0; d < dof; ++d)
            {
              Fe[d + dof * at(i, j, k, Nx, Ny)] = Fe[d + dof * at(i, j, kpb, Nx, Ny)]; 
            }
          }
        }
      }
      #pragma omp for collapse(3)
      for (unsigned int k = 0; k < dend; ++k)
      {
        for (unsigned int j = 0; j < Ny; ++j)
        {
          // now copy up to down
          for (unsigned int i = 0; i < Nx; ++i)
          {
            unsigned int kpb = k + Nz_wrap;
            #pragma omp simd 
            for (unsigned int d = 0; d < dof; ++d)
            {
              Fe[d + dof * at(i, j, k, Nx, Ny)] = Fe[d + dof * at(i, j, kpb, Nx, Ny)]; 
            }
          }
        }
      }
    }
  }
  else
  {
    // get bc for left end of z 
    const BC* bc_zl = &(BCs[4 * dof]);
    // multiplier to enforce bc
    double s;
    for (unsigned int d = 0; d < dof; ++d)
    {
      // we only do something if bc is not none
      if (bc_zl[d] != none)
      { 
        if (bc_zl[d] == mirror) s = 1.0;
        if (bc_zl[d] == mirror_inv) s = -1.0;
        // copy eulerian data in x-y plane in interior region to 
        // adjacent ghost region at down end of z axis according to mirror or mirror_inv
        #pragma omp parallel 
        {
          #pragma omp for collapse(3)
          for (unsigned int k = 0; k < dend; ++k)
          {
            for (unsigned int j = 0; j < Ny; ++j)
            {
              for (unsigned int i = 0; i < Nx; ++i)
              {
                unsigned int kpb = 2 * dend - k;
                Fe[d + dof * at(i, j, k, Nx, Ny)] = s * Fe[d + dof * at(i, j, kpb, Nx, Ny)]; 
              }
            }
          }
          // now handle k=dend (for which we add, not copy)
          #pragma omp for collapse(2)
          for (unsigned int j = 0; j < Ny; ++j)
          {
            for (unsigned int i = 0; i < Nx; ++i)
            {
              Fe[d + dof * at(i, j, dend, Nx, Ny)] += s * Fe[d + dof * at(i, j, dend, Nx, Ny)]; 
            }
          }
        }
        
      }
    }
    // get bc for right end of z
    const BC* bc_zr = &(BCs[5 * dof]);  
    for (unsigned int d = 0; d < dof; ++d)
    {
      // we only do something if bc is not none
      if (bc_zr[d] != none)
      { 
        if (bc_zr[d] == mirror) s = 1.0;
        if (bc_zr[d] == mirror_inv) s = -1.0;
        // copy eulerian data in x-y plane in interior region to 
        // adjacent ghost region at up end of z axis according to mirror or mirror_inv
        #pragma omp parallel
        {
          #pragma omp for collapse(3)
          for (unsigned int k = ubeg - 1; k < Nz; ++k)
          {
            for (unsigned int j = 0; j < Ny; ++j)
            {
              for (unsigned int i = 0; i < Nx; ++i)
              {
                unsigned int kpb = 2 * ubeg - k - 2;
                Fe[d + dof * at(i, j, k, Nx, Ny)] = s * Fe[d + dof * at(i, j, kpb, Nx, Ny)];
              }
            }
          }
          // now handle k=ubeg-1 (for which we add, not copy)
          #pragma omp for collapse(2)
          for (unsigned int j = 0; j < Ny; ++j)
          {
            for (unsigned int i = 0; i < Nx; ++i)
            {
              Fe[d + dof * at(i, j, ubeg-1, Nx, Ny)] += s * Fe[d + dof * at(i, j, ubeg-1, Nx, Ny)];
            }
          } 
        }
      }
    }
  }
}
#endif
