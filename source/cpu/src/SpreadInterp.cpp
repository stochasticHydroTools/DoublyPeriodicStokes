#include"SpreadInterp.h"
#include"Grid.h"
#include"ParticleList.h"
#include"exceptions.h"
#include<omp.h>
#include<algorithm>
#ifndef USE_STACK
#include<fftw3.h>
#endif

void spreadUnifZ(ParticleList& particles, Grid& grid)
{
  // loop over unique alphas
  for (const double& alphaf : particles.unique_threshP)
  {
    const unsigned short wx = std::round(2 * alphaf / grid.hx);
    const unsigned short wy = std::round(2 * alphaf / grid.hy);
    const unsigned short wz = std::round(2 * alphaf / grid.hz);
    const unsigned short w2 = wx * wy;
    const unsigned int kersz = w2 * wz; 
    const unsigned int subsz = w2 * grid.Nzeff;
    const int evenx = -1 * (wx % 2) + 1, eveny = -1 * (wy % 2) + 1;
    // loop over w^2 groups of columns
    for (unsigned int izero = 0; izero < wx; ++izero)
    {
      for (unsigned int jzero = 0; jzero < wy; ++jzero)
      {
        // parallelize over the N^2/w^2 columns in a group
        #pragma omp parallel for collapse(2)
        for (unsigned int ii = izero; ii < grid.Nxeff; ii += wx)
        {
          for (unsigned int jj = jzero; jj < grid.Nyeff; jj += wy)
          {
            // number of pts in this column
            unsigned int npts = grid.number[jj + ii * grid.Nyeff];
            // find first particle in column(ii,jj) with matching alpha 
            int l = grid.firstn[jj + ii * grid.Nyeff];
            while (l >= 0 && particles.threshP[l] != alphaf) 
            {
              l = grid.nextn[l];
              npts -= 1;
            }
            // continue if it's there
            if (l >= 0 && particles.threshP[l] == alphaf)
            {
              // global indices of wx x wy x Nz subarray influenced by column(i,j)
              #ifndef USE_STACK
              unsigned int* indc3D = (unsigned int*) fftw_malloc(subsz * sizeof(unsigned int));
              #else
              unsigned int indc3D[subsz];
              #endif
              for (int k3D = 0; k3D < grid.Nzeff; ++k3D)
              {
                for (int j = 0; j < wy; ++j)
                {
                  int j3D = jj + j - wy / 2 + eveny;
                  for (int i = 0; i < wx; ++i) 
                  {
                    int i3D = ii + i - wx / 2 + evenx;
                    indc3D[at(i,j,k3D,wx,wy)] = at(i3D, j3D, k3D, grid.Nxeff, grid.Nyeff);
                  }
                }
              }
              // gather forces from grid subarray
              #ifndef USE_STACK
              double* fGc = (double*) fftw_malloc(subsz * grid.dof * sizeof(double));
              #else
              double fGc[subsz * grid.dof];
              #endif
              gather(subsz, fGc, grid.fG_unwrap, indc3D, grid.dof);
              // particle indices
              unsigned int npts_match = 1, count  = 1; int ltmp = l;
              // get other particles in col with this alphaf
              for (unsigned int ipt = 1; ipt < npts; ++ipt) 
              {
                ltmp = grid.nextn[ltmp];
                if (particles.threshP[ltmp] == alphaf) {npts_match += 1;}
              }
              #ifndef USE_STACK
              unsigned int* indx = (unsigned int*) fftw_malloc(npts_match * sizeof(unsigned int));
              #else
              unsigned int indx[npts_match];
              #endif
              indx[0] = l;
              for (unsigned int ipt = 1; ipt < npts; ++ipt)
              {
                l = grid.nextn[l];
                if (particles.threshP[l] == alphaf) {indx[count] = l; count += 1;}
              }

              // gather particle pts, betas, forces etc. for this column
              #ifndef USE_STACK
              double *fPc, *betafPc, *normfPc, *xunwrap, *yunwrap, *zunwrap, *alphafPc;
              unsigned int* zoffset; unsigned short* wfPc;
              fPc = (double*) fftw_malloc(npts_match * particles.dof * sizeof(double));
              betafPc = (double*) fftw_malloc(npts_match * sizeof(double));
              alphafPc= (double*) fftw_malloc(npts_match * sizeof(double));
              wfPc = (unsigned short*) fftw_malloc(npts_match * sizeof(unsigned short));
              normfPc = (double*) fftw_malloc(npts_match * sizeof(double)); 
              xunwrap = (double*) fftw_malloc(particles.wfxP_max * npts_match * sizeof(double)); 
              yunwrap = (double*) fftw_malloc(particles.wfyP_max * npts_match * sizeof(double)); 
              zunwrap = (double*) fftw_malloc(particles.wfzP_max * npts_match * sizeof(double));
              zoffset = (unsigned int*) fftw_malloc(npts_match * sizeof(unsigned int));
              #else
              double fPc[npts_match * particles.dof];
              double betafPc[npts_match];
              double alphafPc[npts_match];
              unsigned short wfPc[npts_match];
              double normfPc[npts_match];
              double xunwrap[particles.wfxP_max * npts_match];
              double yunwrap[particles.wfyP_max * npts_match];
              double zunwrap[particles.wfzP_max * npts_match];
              unsigned int zoffset[npts_match];
              #endif
              gather(npts_match, betafPc, particles.betafP, indx, 1);
              gather(npts_match, alphafPc, particles.alphafP, indx, 1);
              gather(npts_match, fPc, particles.fP, indx, particles.dof);
              gather(npts_match, normfPc, particles.normfP, indx, 1);
              gather(npts_match, wfPc, particles.wfP, indx, 1);
              gather(npts_match, xunwrap, particles.xunwrap, indx, particles.wfxP_max);
              gather(npts_match, yunwrap, particles.yunwrap, indx, particles.wfyP_max);
              gather(npts_match, zunwrap, particles.zunwrap, indx, particles.wfzP_max);
              gather(npts_match, zoffset, particles.zoffset, indx, 1);

              // get the kernel w x w x w kernel weights for each particle in col
              #ifndef USE_STACK 
              double* delta = (double*) fftw_malloc(kersz * npts_match * sizeof(double));
              #else
              double delta[kersz * npts_match];
              #endif
              delta_eval_col(delta, betafPc, wfPc, normfPc, xunwrap, yunwrap, 
                             zunwrap, alphafPc, npts_match, wx, wy, wz, particles.wfxP_max,
                             particles.wfyP_max, particles.wfzP_max);

              // spread the particle forces with the kernel weights
              spread_col(fGc, delta, fPc, zoffset, npts_match, kersz, grid.dof);

              // scatter back to global eulerian grid
              scatter(subsz, fGc, grid.fG_unwrap, indc3D, grid.dof);
              #ifndef USE_STACK
              fftw_free(fPc); fPc = 0; fftw_free(betafPc); 
              betafPc = 0; fftw_free(wfPc); wfPc = 0; fftw_free(normfPc); normfPc = 0; 
              fftw_free(xunwrap); xunwrap = 0; fftw_free(yunwrap); yunwrap = 0;
              fftw_free(zunwrap); zunwrap = 0; fftw_free(zoffset); zoffset = 0; 
              fftw_free(delta); delta = 0; fftw_free(indc3D); indc3D = 0; fftw_free(fGc); 
              fGc = 0; fftw_free(indx); indx = 0; fftw_free(alphafPc); alphafPc = 0;
              #endif
            } // finished with column
          } 
        } // finished with group of columns
      }
    } // finished with all groups
  } // finished with this alphaf
}

void interpUnifZ(ParticleList& particles, Grid& grid)
{
  // loop over unique alphas
  for (const double& alphaf : particles.unique_threshP)
  {
    const unsigned short wx = std::round(2 * alphaf / grid.hx);
    const unsigned short wy = std::round(2 * alphaf / grid.hy);
    const unsigned short wz = std::round(2 * alphaf / grid.hz);
    const unsigned short w2 = wx * wy;
    const unsigned int kersz = w2 * wz; 
    const unsigned int subsz = w2 * grid.Nzeff;
    const int evenx = -1 * (wx % 2) + 1, eveny = -1 * (wy % 2) + 1;
    const double weight = grid.hx * grid.hy * grid.hz;
    // loop over w^2 groups of columns
    for (unsigned int izero = 0; izero < wx; ++izero)
    {
      for (unsigned int jzero = 0; jzero < wy; ++jzero)
      {
        // parallelize over the N^2/w^2 columns in a group
        #pragma omp parallel for collapse(2)
        for (unsigned int ii = izero; ii < grid.Nxeff; ii += wx)
        {
          for (unsigned int jj = jzero; jj < grid.Nyeff; jj += wy)
          { 
            // number of pts in this column
            unsigned int npts = grid.number[jj + ii * grid.Nyeff];
            // find first particle in column(ii,jj) with matching alpha 
            int l = grid.firstn[jj + ii * grid.Nyeff];
            while (l >= 0 && particles.threshP[l] != alphaf) 
            {
              l = grid.nextn[l];
              npts -= 1;
            }
            // continue if it's there
            if (l >= 0 && particles.threshP[l] == alphaf)
            {
              // global indices of wx x wy x Nz subarray influenced by column(i,j)
              #ifndef USE_STACK
              unsigned int* indc3D = (unsigned int*) fftw_malloc(subsz * sizeof(unsigned int));
              #else
              unsigned int indc3D[subsz];
              #endif
              for (int k3D = 0; k3D < grid.Nzeff; ++k3D)
              {
                for (int j = 0; j < wy; ++j)
                {
                  int j3D = jj + j - wy / 2 + eveny;
                  for (int i = 0; i < wx; ++i) 
                  {
                    int i3D = ii + i - wx / 2 + evenx;
                    indc3D[at(i,j,k3D,wx,wy)] = at(i3D, j3D, k3D, grid.Nxeff, grid.Nyeff);
                  }
                }
              }
              // gather forces from grid subarray
              #ifndef USE_STACK
              double* fGc = (double*) fftw_malloc(subsz * grid.dof * sizeof(double));
              #else
              double fGc[subsz * grid.dof];
              #endif
              gather(subsz, fGc, grid.fG_unwrap, indc3D, grid.dof);
              // particle indices
              unsigned int npts_match = 1, count  = 1; int ltmp = l;
              // get other particles in col with this alphaf
              for (unsigned int ipt = 1; ipt < npts; ++ipt) 
              {
                ltmp = grid.nextn[ltmp];
                if (particles.threshP[ltmp] == alphaf) {npts_match += 1;}
              }
              #ifndef USE_STACK
              unsigned int* indx = (unsigned int*) fftw_malloc(npts_match * sizeof(unsigned int));
              #else
              unsigned int indx[npts_match];
              #endif
              indx[0] = l;
              for (unsigned int ipt = 1; ipt < npts; ++ipt)
              {
                l = grid.nextn[l];
                if (particles.threshP[l] == alphaf) {indx[count] = l; count += 1;}
              }

              // gather particle pts, betas, forces etc. for this column
              #ifndef USE_STACK
              double *fPc, *betafPc, *normfPc, *xunwrap, *yunwrap, *zunwrap, *alphafPc;
              unsigned int* zoffset; unsigned short* wfPc;
              fPc = (double*) fftw_malloc(npts_match * particles.dof * sizeof(double));
              betafPc = (double*) fftw_malloc(npts_match * sizeof(double));
              alphafPc= (double*) fftw_malloc(npts_match * sizeof(double));
              wfPc = (unsigned short*) fftw_malloc(npts_match * sizeof(unsigned short));
              normfPc = (double*) fftw_malloc(npts_match * sizeof(double)); 
              xunwrap = (double*) fftw_malloc(particles.wfxP_max * npts_match * sizeof(double)); 
              yunwrap = (double*) fftw_malloc(particles.wfyP_max * npts_match * sizeof(double)); 
              zunwrap = (double*) fftw_malloc(particles.wfzP_max * npts_match * sizeof(double));
              zoffset = (unsigned int*) fftw_malloc(npts_match * sizeof(unsigned int));
              #else
              double fPc[npts_match * particles.dof];
              double betafPc[npts_match];
              double alphafPc[npts_match];
              unsigned short wfPc[npts_match];
              double normfPc[npts_match];
              double xunwrap[particles.wfxP_max * npts_match];
              double yunwrap[particles.wfyP_max * npts_match];
              double zunwrap[particles.wfzP_max * npts_match];
              unsigned int zoffset[npts_match];
              #endif
              gather(npts_match, betafPc, particles.betafP, indx, 1);
              gather(npts_match, alphafPc, particles.alphafP, indx, 1);
              gather(npts_match, fPc, particles.fP, indx, particles.dof);
              gather(npts_match, normfPc, particles.normfP, indx, 1);
              gather(npts_match, wfPc, particles.wfP, indx, 1);
              gather(npts_match, xunwrap, particles.xunwrap, indx, particles.wfxP_max);
              gather(npts_match, yunwrap, particles.yunwrap, indx, particles.wfyP_max);
              gather(npts_match, zunwrap, particles.zunwrap, indx, particles.wfzP_max);
              gather(npts_match, zoffset, particles.zoffset, indx, 1);

              // get the kernel w x w x w kernel weights for each particle in col 
              #ifndef USE_STACK
              double* delta = (double*) fftw_malloc(kersz * npts_match * sizeof(double));
              #else
              double delta[kersz * npts_match];
              #endif
              delta_eval_col(delta, betafPc, wfPc, normfPc, xunwrap, yunwrap, 
                             zunwrap, alphafPc, npts_match, wx, wy, wz, particles.wfxP_max,
                             particles.wfyP_max, particles.wfzP_max);

              // interpolate on the particles with the kernel weights
              interp_col(fGc, delta, fPc, zoffset, npts_match, kersz, grid.dof, weight);

              // scatter back to global lagrangian grid
              scatter(npts_match, fPc, particles.fP, indx, particles.dof);
              #ifndef USE_STACK
              fftw_free(fPc); fPc = 0; fftw_free(betafPc); 
              betafPc = 0; fftw_free(wfPc); wfPc = 0; fftw_free(normfPc); normfPc = 0; 
              fftw_free(xunwrap); xunwrap = 0; fftw_free(yunwrap); yunwrap = 0;
              fftw_free(zunwrap); zunwrap = 0; fftw_free(zoffset); zoffset = 0; 
              fftw_free(delta); delta = 0; fftw_free(indc3D); indc3D = 0; fftw_free(fGc); 
              fGc = 0; fftw_free(indx); indx = 0; fftw_free(alphafPc); alphafPc = 0;
              #endif
            } // finished with column
          } 
        } // finished with group of columns
      }
    } // finished with all groups
  } // finished with this alphaf
}

void spreadNonUnifZ(ParticleList& particles, Grid& grid)
{
  // loop over unique alphas
  for (const double& alphaf : particles.unique_threshP)
  {
    const unsigned short wx = std::round(2 * alphaf / grid.hx);
    const unsigned short wy = std::round(2 * alphaf / grid.hy);
    const unsigned short w2 = wx * wy;
    const unsigned int subsz = w2 * grid.Nzeff;
    const int evenx = -1 * (wx % 2) + 1, eveny = -1 * (wy % 2) + 1;
    // loop over w^2 groups of columns
    for (unsigned int izero = 0; izero < wx; ++izero)
    {
      for (unsigned int jzero = 0; jzero < wy; ++jzero)
      {
        // parallelize over the N^2/w^2 columns in a group
        #pragma omp parallel for collapse(2)
        for (unsigned int ii = izero; ii < grid.Nxeff; ii += wx)
        {
          for (unsigned int jj = jzero; jj < grid.Nyeff; jj += wy)
          {
            // number of pts in this column
            unsigned int npts = grid.number[jj + ii * grid.Nyeff];
            // find first particle in column(ii,jj) with matching alpha 
            int l = grid.firstn[jj + ii * grid.Nyeff];
            while (l >= 0 && particles.threshP[l] != alphaf) 
            {
              l = grid.nextn[l];
              npts -= 1;
            }
            // continue if it's there
            if (l >= 0 && particles.threshP[l] == alphaf)
            {
              // global indices of wx x wy x Nz subarray influenced by column(i,j)
              #ifndef USE_STACK
              unsigned int* indc3D = (unsigned int*) fftw_malloc(subsz * sizeof(unsigned int));
              #else
              unsigned int indc3D[subsz];
              #endif
              for (int k3D = 0; k3D < grid.Nzeff; ++k3D)
              {
                for (int j = 0; j < wy; ++j)
                {
                  int j3D = jj + j - wy / 2 + eveny;
                  for (int i = 0; i < wx; ++i) 
                  {
                    int i3D = ii + i - wx / 2 + evenx;
                    indc3D[at(i,j,k3D,wx,wy)] = at(i3D, j3D, k3D, grid.Nxeff, grid.Nyeff);
                  }
                }
              }
              // gather forces from grid subarray
              #ifndef USE_STACK
              double* fGc = (double*) fftw_malloc(subsz * grid.dof * sizeof(double));
              #else
              double fGc[subsz * grid.dof];
              #endif
              gather(subsz, fGc, grid.fG_unwrap, indc3D, grid.dof);
              // particle indices
              unsigned int npts_match = 1, count  = 1; int ltmp = l;
              // get other particles in col with this alphaf
              for (unsigned int ipt = 1; ipt < npts; ++ipt) 
              {
                ltmp = grid.nextn[ltmp];
                if (particles.threshP[ltmp] == alphaf) {npts_match += 1;}
              }
              #ifndef USE_STACK
              unsigned int* indx = (unsigned int*) fftw_malloc(npts_match * sizeof(unsigned int));
              #else
              unsigned int indx[npts_match];
              #endif
              indx[0] = l;
              for (unsigned int ipt = 1; ipt < npts; ++ipt)
              {
                l = grid.nextn[l];
                if (particles.threshP[l] == alphaf) {indx[count] = l; count += 1;}
              }

              // gather particle pts, betas, forces etc. for this column
              #ifndef USE_STACK
              double *fPc, *betafPc, *normfPc, *xunwrap, *yunwrap, *zunwrap, *alphafPc;
              unsigned int* zoffset; unsigned short *wfPc, *wz;
              fPc = (double*) fftw_malloc(npts_match * particles.dof * sizeof(double));
              betafPc = (double*) fftw_malloc(npts_match * sizeof(double));
              alphafPc= (double*) fftw_malloc(npts_match * sizeof(double));
              wfPc = (unsigned short*) fftw_malloc(npts_match * sizeof(unsigned short));
              wz = (unsigned short*) fftw_malloc(npts_match * sizeof(unsigned short));
              normfPc = (double*) fftw_malloc(npts_match * sizeof(double)); 
              xunwrap = (double*) fftw_malloc(particles.wfxP_max * npts_match * sizeof(double)); 
              yunwrap = (double*) fftw_malloc(particles.wfyP_max * npts_match * sizeof(double)); 
              zunwrap = (double*) fftw_malloc(particles.wfzP_max * npts_match * sizeof(double));
              zoffset = (unsigned int*) fftw_malloc(npts_match * sizeof(unsigned int));
              #else
              double fPc[npts_match * particles.dof]; 
              double betafPc[npts_match];
              double alphafPc[npts_match];
              unsigned short wfPc[npts_match]; 
              unsigned short wz[npts_match]; 
              double normfPc[npts_match]; 
              double xunwrap[particles.wfxP_max * npts_match]; 
              double yunwrap[particles.wfyP_max * npts_match]; 
              double zunwrap[particles.wfzP_max * npts_match]; 
              unsigned int zoffset[npts_match]; 
              #endif
              gather(npts_match, betafPc, particles.betafP, indx, 1);
              gather(npts_match, alphafPc, particles.alphafP, indx, 1);
              gather(npts_match, fPc, particles.fP, indx, particles.dof);
              gather(npts_match, normfPc, particles.normfP, indx, 1);
              gather(npts_match, wfPc, particles.wfP, indx, 1);
              gather(npts_match, wz, particles.wfzP, indx, 1);
              gather(npts_match, xunwrap, particles.xunwrap, indx, particles.wfxP_max);
              gather(npts_match, yunwrap, particles.yunwrap, indx, particles.wfyP_max);
              gather(npts_match, zunwrap, particles.zunwrap, indx, particles.wfzP_max);
              gather(npts_match, zoffset, particles.zoffset, indx, 1);
              
              const unsigned int kersz = w2 * (*std::max_element(wz, wz + npts_match));

              // get the kernel w x w x w kernel weights for each particle in col 
              #ifndef USE_STACK
              double* delta = (double*) fftw_malloc(kersz * npts_match * sizeof(double));
              #else
              double delta[kersz * npts_match];
              #endif
              delta_eval_col(delta, betafPc, wfPc, normfPc, xunwrap, yunwrap, 
                             zunwrap, alphafPc, npts_match, wx, wy, wz, particles.wfxP_max,
                             particles.wfyP_max, particles.wfzP_max);

              // spread the particle forces with the kernel weights
              spread_col(fGc, delta, fPc, zoffset, npts_match, w2, wz, grid.dof);

              // scatter back to global eulerian grid
              scatter(subsz, fGc, grid.fG_unwrap, indc3D, grid.dof);
              #ifndef USE_STACK
              fftw_free(fPc); fPc = 0; fftw_free(betafPc); fftw_free(wz); wz = 0; 
              betafPc = 0; fftw_free(wfPc); wfPc = 0; fftw_free(normfPc); normfPc = 0; 
              fftw_free(xunwrap); xunwrap = 0; fftw_free(yunwrap); yunwrap = 0;
              fftw_free(zunwrap); zunwrap = 0; fftw_free(zoffset); zoffset = 0; 
              fftw_free(delta); delta = 0; fftw_free(indc3D); indc3D = 0; fftw_free(fGc); 
              fGc = 0; fftw_free(indx); indx = 0; fftw_free(alphafPc); alphafPc = 0;
              #endif
            } // finished with column
          } 
        } // finished with group of columns
      }
    } // finished with all groups
  } // finished with this alphaf
}

void interpNonUnifZ(ParticleList& particles, Grid& grid)
{
  // loop over unique alphas
  for (const double& alphaf : particles.unique_threshP)
  {
    const unsigned short wx = std::round(2 * alphaf / grid.hx);
    const unsigned short wy = std::round(2 * alphaf / grid.hy);
    const unsigned short w2 = wx * wy;
    const unsigned int subsz = w2 * grid.Nzeff;
    const int evenx = -1 * (wx % 2) + 1, eveny = -1 * (wy % 2) + 1;
    // loop over w^2 groups of columns
    for (unsigned int izero = 0; izero < wx; ++izero)
    {
      for (unsigned int jzero = 0; jzero < wy; ++jzero)
      {
        // parallelize over the N^2/w^2 columns in a group
        #pragma omp parallel for collapse(2)
        for (unsigned int ii = izero; ii < grid.Nxeff; ii += wx)
        {
          for (unsigned int jj = jzero; jj < grid.Nyeff; jj += wy)
          {
            // number of pts in this column
            unsigned int npts = grid.number[jj + ii * grid.Nyeff];
            // find first particle in column(ii,jj) with matching alpha 
            int l = grid.firstn[jj + ii * grid.Nyeff];
            while (l >= 0 && particles.threshP[l] != alphaf) 
            {
              l = grid.nextn[l];
              npts -= 1;
            }
            // continue if it's there
            if (l >= 0 && particles.threshP[l] == alphaf)
            {
              // global indices of wx x wy x Nz subarray influenced by column(i,j)
              #ifndef USE_STACK
              unsigned int* indc3D = (unsigned int*) fftw_malloc(subsz * sizeof(unsigned int));
              #else
              unsigned int indc3D[subsz];
              #endif
              for (int k3D = 0; k3D < grid.Nzeff; ++k3D)
              {
                for (int j = 0; j < wy; ++j)
                {
                  int j3D = jj + j - wy / 2 + eveny;
                  for (int i = 0; i < wx; ++i) 
                  {
                    int i3D = ii + i - wx / 2 + evenx;
                    indc3D[at(i,j,k3D,wx,wy)] = at(i3D, j3D, k3D, grid.Nxeff, grid.Nyeff);
                  }
                }
              }
              // gather forces from grid subarray
              #ifndef USE_STACK 
              double* fGc = (double*) fftw_malloc(subsz * grid.dof * sizeof(double));
              #else
              double fGc[subsz * grid.dof];
              #endif
              gather(subsz, fGc, grid.fG_unwrap, indc3D, grid.dof);
              // particle indices
              unsigned int npts_match = 1, count  = 1; int ltmp = l;
              // get other particles in col with this alphaf
              for (unsigned int ipt = 1; ipt < npts; ++ipt) 
              {
                ltmp = grid.nextn[ltmp];
                if (particles.threshP[ltmp] == alphaf) {npts_match += 1;}
              }
              #ifndef USE_STACK
              unsigned int* indx = (unsigned int*) fftw_malloc(npts_match * sizeof(unsigned int));
              #else
              unsigned int indx[npts_match];
              #endif
              indx[0] = l;
              for (unsigned int ipt = 1; ipt < npts; ++ipt)
              {
                l = grid.nextn[l];
                if (particles.threshP[l] == alphaf) {indx[count] = l; count += 1;}
              }

              // gather particle pts, betas, forces etc. for this column
              #ifndef USE_STACK 
              double *fPc, *betafPc, *normfPc, *xunwrap, *yunwrap, *zunwrap, *pt_wts, *alphafPc;
              unsigned int* zoffset; unsigned short *wfPc, *wz;
              fPc = (double*) fftw_malloc(npts_match * particles.dof * sizeof(double));
              betafPc = (double*) fftw_malloc(npts_match * sizeof(double));
              alphafPc = (double*) fftw_malloc(npts_match * sizeof(double));
              wfPc = (unsigned short*) fftw_malloc(npts_match * sizeof(unsigned short));
              wz = (unsigned short*) fftw_malloc(npts_match * sizeof(unsigned short));
              normfPc = (double*) fftw_malloc(npts_match * sizeof(double)); 
              xunwrap = (double*) fftw_malloc(particles.wfxP_max * npts_match * sizeof(double)); 
              yunwrap = (double*) fftw_malloc(particles.wfyP_max * npts_match * sizeof(double)); 
              zunwrap = (double*) fftw_malloc(particles.wfzP_max * npts_match * sizeof(double));
              pt_wts = (double*) fftw_malloc(particles.wfzP_max * npts_match * sizeof(double));
              zoffset = (unsigned int*) fftw_malloc(npts_match * sizeof(unsigned int));
              #else
              double fPc[npts_match * particles.dof]; 
              double betafPc[npts_match];
              double alphafPc[npts_match];
              unsigned short wfPc[npts_match]; 
              unsigned short wz[npts_match]; 
              double normfPc[npts_match]; 
              double xunwrap[particles.wfxP_max * npts_match]; 
              double yunwrap[particles.wfyP_max * npts_match]; 
              double zunwrap[particles.wfzP_max * npts_match]; 
              double pt_wts[particles.wfzP_max * npts_match];
              unsigned int zoffset[npts_match]; 
              #endif 
              gather(npts_match, betafPc, particles.betafP, indx, 1);
              gather(npts_match, alphafPc, particles.alphafP, indx, 1);
              gather(npts_match, fPc, particles.fP, indx, particles.dof);
              gather(npts_match, normfPc, particles.normfP, indx, 1);
              gather(npts_match, wfPc, particles.wfP, indx, 1);
              gather(npts_match, wz, particles.wfzP, indx, 1);
              gather(npts_match, xunwrap, particles.xunwrap, indx, particles.wfxP_max);
              gather(npts_match, yunwrap, particles.yunwrap, indx, particles.wfyP_max);
              gather(npts_match, zunwrap, particles.zunwrap, indx, particles.wfzP_max);
              gather(npts_match, pt_wts, particles.pt_wts, indx, particles.wfzP_max);
              gather(npts_match, zoffset, particles.zoffset, indx, 1);

              const unsigned int kersz = w2 * (*std::max_element(wz, wz + npts_match));

              // get the kernel w x w x w kernel weights for each particle in col 
              // double delta[kersz * npts_match];
              #ifndef USE_STACK
              double* delta = (double*) fftw_malloc(kersz * npts_match * sizeof(double));
              #else
              double delta[kersz * npts_match];
              #endif
              delta_eval_col(delta, betafPc, wfPc, normfPc, xunwrap, yunwrap, 
                             zunwrap, alphafPc, npts_match, wx, wy, wz, particles.wfxP_max,
                             particles.wfyP_max, particles.wfzP_max);

              // spread the particle forces with the kernel weights
              interp_col(fGc, delta, fPc, zoffset, npts_match, wx, wy, wz, particles.wfzP_max, grid.dof, pt_wts);

              // scatter back to global lagrangian grid
              scatter(npts_match, fPc, particles.fP, indx, particles.dof);
              #ifndef USE_STACK
              fftw_free(fPc); fPc = 0; fftw_free(betafPc); fftw_free(wz); wz = 0; 
              betafPc = 0; fftw_free(wfPc); wfPc = 0; fftw_free(normfPc); normfPc = 0; 
              fftw_free(xunwrap); xunwrap = 0; fftw_free(yunwrap); yunwrap = 0;
              fftw_free(zunwrap); zunwrap = 0; fftw_free(zoffset); zoffset = 0; 
              fftw_free(delta); delta = 0; fftw_free(indc3D); indc3D = 0; fftw_free(fGc); 
              fGc = 0; fftw_free(indx); indx = 0; fftw_free(alphafPc); alphafPc = 0;
              #endif
            } // finished with column
          } 
        } // finished with group of columns
      }
    } // finished with all groups
  } // finished with this alphaf
}
