#ifndef SPREADINTERP_H 
#define SPREADINTERP_H
#include"ESKernels.h"
#include"common.h"

/* Main routines and helpers for spreading and interpolation, with
   the appropriate routines dispatched based on
   specified boundary conditions. 
   NOTES: These routines do not reinitialize the state of the grid
          or the particles. They also do not know about boundaries.
          The caller is responsible for:
          - initializing data on the extended 
            grid (eg. using grid.zeroExtGrid()) before calls to spread().
          - handling boundary conditions (eg. using DeGhostify()) after calls
            to spread(). 
          - initializing data on the particles (eg. using particles.zeroForces()) 
            before calls to interpolate()
          - handling boundary conditions (eg. using Ghostify()) before calls to
            interpolate()
   Side effects:
          - the member grid.fG_unwrap is overwritten during both spread and interp
          - the member particle.fP is overwritten with interp data
*/


// forward declarations
struct Grid;
struct ParticleList;

// spread with z uniform or not
void spreadUnifZ(ParticleList& particles, Grid& grid);
void spreadNonUnifZ(ParticleList& particles, Grid& grid);
// interpolate with z uniform or not
void interpUnifZ(ParticleList& particles, Grid& grid);
void interpNonUnifZ(ParticleList& particles, Grid& grid);


// gather data from src at inds into trg
template<typename T> 
inline void gather(unsigned int N, T* trg, T const* src, 
                   const unsigned int* inds, const unsigned int dof)
{
  for (unsigned int i = 0; i < N; ++i) 
  {
    #pragma omp simd
    for (unsigned int j = 0; j < dof; ++j)
    {
      trg[j + dof * i] = src[j + dof * inds[i]];
    }
  }
}

// scatter data from trg into src at inds
template<typename T>
inline void scatter(unsigned int N, T const* trg, T* src, 
                    const unsigned int* inds, const unsigned int dof)
{
  for (unsigned int i = 0; i < N; ++i) 
  {
    #pragma omp simd
    for (unsigned int j = 0; j < dof; ++j)
    {
      src[j + dof * inds[i]] = trg[j + dof * i];
    }
  }
}

// evaluate the delta function weights for the current column for UnifZ = True
inline void delta_eval_col(double* delta, const double* betafPc,
                           const unsigned short* wfPc, const double* normfPc, 
                           const double* xunwrap, const double* yunwrap, 
                           const double* zunwrap, const double* alphafP, const int npts, 
                           const unsigned short wx, const unsigned short wy, 
                           const unsigned short wz, const unsigned short wfxP_max,
                           const unsigned short wfyP_max, const unsigned short wfzP_max)
{
  double x[3];
  for (unsigned int k = 0; k < wz; ++k)
  {
    for (unsigned int j = 0; j < wy; ++j)
    {
      for (unsigned int i = 0; i < wx; ++i)
      {
        #pragma omp simd
        for (unsigned int ipt = 0; ipt < npts; ++ipt)
        {
          unsigned int m = at(i, j, k, wx, wy);
          double norm = normfPc[ipt]; norm *= norm * norm;;
          x[0] = xunwrap[i + ipt * wfxP_max];
          x[1] = yunwrap[j + ipt * wfyP_max];
          x[2] = zunwrap[k + ipt * wfzP_max];
          delta[ipt + m * npts] = esKernel(x, betafPc[ipt] * wfPc[ipt], alphafP[ipt]) / norm;
        }
      }
    }
  }
}

// evaluate the delta function weights for the current column for UnifZ = false
inline void delta_eval_col(double* delta, const double* betafPc,
                           const unsigned short* wfPc, const double* normfPc, 
                           const double* xunwrap, const double* yunwrap, 
                           const double* zunwrap, const double* alphafP, const int npts, 
                           const unsigned short wx, const unsigned short wy, 
                           const unsigned short* wz, const unsigned short wfxP_max,
                           const unsigned short wfyP_max, const unsigned short wfzP_max)
{
  double x[3];
  for (unsigned int ipt = 0; ipt < npts; ++ipt)
  {
    for (unsigned int k = 0; k < wz[ipt]; ++k)
    {
      for (unsigned int j = 0; j < wy; ++j)
      {
        #pragma omp simd
        for (unsigned int i = 0; i < wx; ++i)
        {
          unsigned int m = at(i, j, k, wx, wy);
          double norm = normfPc[ipt]; norm *= norm * norm;;
          x[0] = xunwrap[i + ipt * wfxP_max];
          x[1] = yunwrap[j + ipt * wfyP_max];
          x[2] = zunwrap[k + ipt * wfzP_max];
          delta[ipt + m * npts] = esKernel(x, betafPc[ipt] * wfPc[ipt], alphafP[ipt]) / norm;
        }
      }
    }
  }
}

// spread the delta functions weights for the column for UnifZ = true
inline void spread_col(double* Fec, const double* delta, const double* flc,
                       const unsigned int* zoffset, const int npts,
                       const int w3, const int dof)
{
  for (unsigned int ipt = 0; ipt < npts; ++ipt)
  {
    for (unsigned int i = 0; i < w3; ++i)
    {
      #pragma omp simd
      for (unsigned int j = 0; j < dof; ++j)
      {
        Fec[j + dof * (i + zoffset[ipt])] += delta[ipt + i * npts] * flc[j + dof * ipt];
      }
    }
  }
}

// spread with forces and weights for the column for UnifZ = false
inline void spread_col(double* Fec, const double* delta, const double* flc,
                       const unsigned int* zoffset, const int npts,
                       const int w2, const unsigned short* wz, const int dof)
{
  for (unsigned int ipt = 0; ipt < npts; ++ipt)
  {
    for (unsigned int i = 0; i < w2 * wz[ipt]; ++i)
    {
      #pragma omp simd
      for (unsigned int j = 0; j < dof; ++j)
      {
        Fec[j + dof * (i + zoffset[ipt])] += delta[ipt + i * npts] * flc[j + dof * ipt];
      }
    }
  }
}

// interpolate with the forces and weights for the current column for UNIFORM Z
inline void interp_col(const double* Fec, const double* delta, double* flc, 
                       const unsigned int* zoffset, const int npts, 
                       const int w3, const int dof, const double weight)
{
  for (unsigned ipt = 0; ipt < npts; ++ipt)
  {
    for (unsigned int i = 0; i < w3; ++i)
    {
      #pragma omp simd
      for (unsigned int j = 0; j < dof; ++j)
      { 
        flc[j + dof * ipt] += Fec[j + dof * (i + zoffset[ipt])] * 
                                delta[ipt + i * npts] * weight; 
      }
    }
  }
}


// interpolate with the forces and weights for the current column for NON-UNIFORM Z
inline void interp_col(const double* Fec, const double* delta, double* flc, 
                       const unsigned int* zoffset, const int npts, 
                       const unsigned short wx, const unsigned short wy, 
                       const unsigned short* wz, const unsigned short wfzP_max, 
                       const int dof, const double* weight)
{
  for (unsigned ipt = 0; ipt < npts; ++ipt)
  {
    for (unsigned int k = 0; k < wz[ipt]; ++k)
    {
      for (unsigned int j = 0; j < wy; ++j)
      {
        for (unsigned int i = 0; i < wx; ++i)
        {
          unsigned int m = at(i, j, k, wx, wy);
          #pragma omp simd
          for (unsigned int d = 0; d < dof; ++d)
          {
            flc[d + dof * ipt] += Fec[d + dof * (m + zoffset[ipt])] * 
                                    delta[ipt + m * npts] * weight[k + ipt * wfzP_max];
          } 
        }
      }
    }
  }
}

#endif
