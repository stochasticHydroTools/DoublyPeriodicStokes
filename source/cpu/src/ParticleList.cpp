#include<unordered_set>
#include<algorithm>
#include<fstream>
#include<iomanip>
#include<random>
#include<omp.h>
#include<math.h>
#include<complex.h>
#include<fftw3.h>
#include"ParticleList.h"
#include"ESKernels.h"
#include"Grid.h"
#include"Quadrature.h"
#include"exceptions.h"

// null initialization
ParticleList::ParticleList() : xP(0), fP(0), betafP(0), alphafP(0), threshP(0),
                             radP(0), normfP(0), wfP(0), wfxP(0), wfyP(0),
                             wfzP(0), nP(0), normalized(false), rad_unknown(true), dof(0), 
                             unique_monopoles(ESParticleSet(20,esparticle_hash)),
                             xunwrap(0), yunwrap(0), zunwrap(0), zoffset(0), pt_wts(0),
                             n_threads(1), isDipole(false), isCbeta(false), indl(0), indr(0)
{}

/* construct with external data by copy */
ParticleList::ParticleList(const double* _xP, const double* _fP, const double* _radP, 
                         const double* _betafP, const double* _cwfP, const unsigned short* _wfP, 
                         const unsigned int _nP, const unsigned int _dof, const bool _isDipole, 
                         const bool isCbeta) :
  nP(_nP), dof(_dof), normalized(false), unique_monopoles(ESParticleSet(20,esparticle_hash)), 
  xunwrap(0), yunwrap(0), zunwrap(0), zoffset(0), pt_wts(0), n_threads(1), rad_unknown(false), 
  isDipole(_isDipole), isCbeta(isCbeta)
{
  xP = (double*) fftw_malloc(nP * 3 * sizeof(double));
  fP = (double*) fftw_malloc(nP * dof * sizeof(double));
  betafP = (double*) fftw_malloc(nP * sizeof(double));
  radP = (double*) fftw_malloc(nP * sizeof(double));
  cwfP = (double*) fftw_malloc(nP * sizeof(double));
  alphafP = (double*) fftw_malloc(nP * sizeof(double));
  threshP = (double*) fftw_malloc(nP * sizeof(double));
  normfP = (double*) fftw_malloc(nP * sizeof(double));
  wfP = (unsigned short*) fftw_malloc(nP * sizeof(unsigned short));
  wfxP = (unsigned short*) fftw_malloc(nP * sizeof(unsigned short));
  wfyP = (unsigned short*) fftw_malloc(nP * sizeof(unsigned short));
  wfzP = (unsigned short*) fftw_malloc(nP * sizeof(unsigned short));
  for (unsigned int i = 0; i < nP; ++i)
  {
    xP[3 * i] = _xP[3 * i];
    xP[1 + 3 * i] = _xP[1 + 3 * i];
    xP[2 + 3 * i] = _xP[2 + 3 * i];
    radP[i] = _radP[i];
    betafP[i] = _betafP[i];
    cwfP[i] = _cwfP[i];
    wfP[i] = _wfP[i];
    for (unsigned int j = 0; j < dof; ++j)
    {
      fP[j + dof * i] = _fP[j + dof * i];
    }
  }
  this->setup();
}


/* construct partially with external data by copy */
ParticleList::ParticleList(const double* _xP, const double* _fP, const double* _betafP, const double* _alphafP,
                           const unsigned short* _wfP, const unsigned int _nP, const unsigned int _dof, 
                           const bool _isDipole) :
  nP(_nP), dof(_dof), normalized(false), unique_monopoles(ESParticleSet(0,esparticle_hash)), 
  xunwrap(0), yunwrap(0), zunwrap(0), zoffset(0), pt_wts(0), n_threads(1), rad_unknown(true), 
  isDipole(_isDipole), isCbeta(false)
{
  xP = (double*) fftw_malloc(nP * 3 * sizeof(double));
  fP = (double*) fftw_malloc(nP * dof * sizeof(double));
  betafP = (double*) fftw_malloc(nP * sizeof(double));
  radP = (double*) fftw_malloc(nP * sizeof(double));
  cwfP = (double*) fftw_malloc(nP * sizeof(double));
  alphafP = (double*) fftw_malloc(nP * sizeof(double));
  threshP = (double*) fftw_malloc(nP * sizeof(double));
  normfP = (double*) fftw_malloc(nP * sizeof(double));
  wfP = (unsigned short*) fftw_malloc(nP * sizeof(unsigned short));
  wfxP = (unsigned short*) fftw_malloc(nP * sizeof(unsigned short));
  wfyP = (unsigned short*) fftw_malloc(nP * sizeof(unsigned short));
  wfzP = (unsigned short*) fftw_malloc(nP * sizeof(unsigned short)); 

  for (unsigned int i = 0; i < nP; ++i)
  {
    xP[3 * i] = _xP[3 * i];
    xP[1 + 3 * i] = _xP[1 + 3 * i];
    xP[2 + 3 * i] = _xP[2 + 3 * i];
    betafP[i] = _betafP[i];
    alphafP[i] = _alphafP[i];
    wfP[i] = _wfP[i];
    for (unsigned int j = 0; j < dof; ++j)
    {
      fP[j + dof * i] = _fP[j + dof * i];
    }
  }
  // cheb grid size, pts, weights and f for kernel vals
  unsigned int N = 1000;
  double* cpts = (double*) fftw_malloc(N * sizeof(double));
  double* cwts = (double*) fftw_malloc(N * sizeof(double)); 
  double* f = (double*) fftw_malloc(N * sizeof(double));
  double alpha, beta, betaw, norm; unsigned short w;
  // iterate over unique tuples of (w, beta, c(w), Rh)
  for (unsigned int i = 0; i < nP; ++i)
  {
    w = wfP[i]; beta = betafP[i];
    betaw = beta * w;
    alpha = alphafP[i];
    clencurt(cpts, cwts, -1 * alpha, alpha, N);
    #pragma omp simd
    for (unsigned int j = 0; j < N; ++j)
    {
      f[j] = eskernel(cpts[j], betaw, alpha);
    }
    norm = 0;
    // get normalization for this particle type (unique tuple)
    #pragma omp simd
    for (unsigned int j = 0; j < N; ++j) {norm += f[j] * cwts[j];} 
    // assign the normalization to particles 
    normfP[i] = norm;
    // if we use the es kernel, phi(z), we need z <= alpha 
    if (not this->isDipole){threshP[i] = alpha;}
    // if we use the derivative phi'(z), we need z <= thresh
    else
    {
      /*
        Deriv of es kernel phi'(z) = -beta * z * phi(z) / (alpha^2 * sqrt(1 - (z/alpha)^2))

        We must have z <= thresh so that we avoid the singularity at z = alpha. 
        That is, 1 \approx thresh/alpha < 1, and we calculate thresh for each 
        particle below
      */
      threshP[i] =  creal(csqrt(0.6e1 * cpow(((8 * betaw * betaw) + 
                          0.12e2 * csqrt((-12 * betaw * betaw + 81)) -  
                          0.108e3) * betaw, 0.1e1 / 0.3e1) / betaw +  
                          0.24e2 * betaw * cpow(((8 * betaw * betaw) + 0.12e2 *  
                          csqrt((-12 * betaw * betaw + 81)) - 0.108e3) * betaw, 
                          -0.1e1 / 0.3e1) + 0.12e2) * alpha / 0.6e1);
    }
    this->unique_threshP.emplace(threshP[i]);
  }
  fftw_free(cpts);
  fftw_free(cwts);
  fftw_free(f);
  this->normalized = true;
}

void ParticleList::setData(const double* _fP, unsigned int _dof)
{
  if (!dof) this->dof = _dof;
  else if (this->dof != _dof) exitErr("DOF does not match current.");
  if(!this->fP)
  {
    this->fP = (double*) fftw_malloc(nP * dof * sizeof(double));
  }
  memcpy(fP, _fP, nP * dof * sizeof(double));
}

void ParticleList::zeroData()
{
  if (this->fP)
  {
    memset(fP, 0, dof * nP * sizeof(double));
  }
  else exitErr("Forces have not been allocated.");
}

void ParticleList::setup()
{
  if (this->validState())
  {
    this->findUniqueKernels();
    this->normalizeKernels();
  }
  else
  {
    exitErr("ParticleList is invalid.");
  }
}

void ParticleList::setup(Grid& grid)
{
  if (grid.validState())
  {
    if (dof != grid.dof)
    {
      exitErr("DOF of ParticleList must match DOF of grid.");
    }
    this->setup();
    this->locateOnGrid(grid);
  }
  else
  {
    exitErr("Particles could not be setup on grid because grid is invalid.");
  }
}

void ParticleList::findUniqueKernels()
{
  if (this->unique_monopoles.size() == 0 and not this->rad_unknown)
  {
    for (unsigned int i = 0; i < nP; ++i)
    {
      this->unique_monopoles.emplace(wfP[i], betafP[i], cwfP[i], radP[i]); 
    }
  }
}

/* normalize ES kernels using clenshaw-curtis quadrature*/
void ParticleList::normalizeKernels()
{
  // proceed if we haven't already normalized
  if (not this->normalized)
  {
    // cheb grid size, pts, weights and f for kernel vals
    unsigned int N = 1000;
    double* cpts = (double*) fftw_malloc(N * sizeof(double));
    double* cwts = (double*) fftw_malloc(N * sizeof(double)); 
    double* f = (double*) fftw_malloc(N * sizeof(double));
    double alpha, beta, betaw, norm, cw, rad; unsigned short w;
    // iterate over unique tuples of (w, beta, c(w), Rh)
    for (const auto& tuple : unique_monopoles)
    {
      w = std::get<0>(tuple); beta = std::get<1>(tuple);
      betaw = beta * w;
      cw = std::get<2>(tuple);
      rad = std::get<3>(tuple); 
      if (this->isCbeta) {alpha = rad / (2 * cw);}
      else {alpha = w * rad / (2 * cw);}

      clencurt(cpts, cwts, -1 * alpha, alpha, N);
      #pragma omp simd
      for (unsigned int j = 0; j < N; ++j)
      {
        f[j] = eskernel(cpts[j], betaw, alpha);
      }
      norm = 0;
      // get normalization for this particle type (unique tuple)
      #pragma omp simd
      for (unsigned int j = 0; j < N; ++j) {norm += f[j] * cwts[j];} 
      // assign the normalization to particles with this type
      for (unsigned int i = 0; i < this->nP; ++i)
      {
        if (wfP[i] == w && betafP[i] == beta && cwfP[i] == cw && radP[i] == rad) 
        {
          normfP[i] = norm; alphafP[i] = alpha;
          // if we use the es kernel, phi(z), we need z <= alpha 
          if (not this->isDipole){threshP[i] = alpha;}
          // if we use the derivative phi'(z), we need z <= thresh
          else
          {
            /*
              Deriv of es kernel phi'(z) = -beta * z * phi(z) / (alpha^2 * sqrt(1 - (z/alpha)^2))

              We must have z <= thresh so that we avoid the singularity at z = alpha. 
              That is, 1 \approx thresh/alpha < 1, and we calculate thresh for each 
              particle below
            */
            threshP[i] =  creal(csqrt(0.6e1 * cpow(((8 * betaw * betaw) + 
                                0.12e2 * csqrt((-12 * betaw * betaw + 81)) - 
                                0.108e3) * betaw, 0.1e1 / 0.3e1) / betaw + 
                                0.24e2 * betaw * cpow(((8 * betaw * betaw) + 0.12e2 * 
                                csqrt((-12 * betaw * betaw + 81)) - 0.108e3) * betaw, 
                                -0.1e1 / 0.3e1) + 0.12e2) * alpha / 0.6e1);
          }
          this->unique_threshP.emplace(threshP[i]);
        }
      }
    }
    fftw_free(cpts);
    fftw_free(cwts);
    fftw_free(f);
    this->normalized = true;
  }
}

void ParticleList::locateOnGrid(Grid& grid)
{
  if (grid.unifZ) {this->locateOnGridUnifZ(grid);}
  else {this->locateOnGridNonUnifZ(grid);}
  grid.has_locator = true; 
}

void ParticleList::locateOnGridUnifZ(Grid& grid)
{
  // get widths on effective uniform grid
  #pragma omp parallel for num_threads(n_threads)
  for (unsigned int i = 0; i < nP; ++i)
  {
    wfxP[i] = std::round(2 * threshP[i] / grid.hx);
    wfyP[i] = std::round(2 * threshP[i] / grid.hy);
    wfzP[i] = std::round(2 * threshP[i] / grid.hz);
    // shift to [0,L]x..
    xP[3 * i] -= grid.minX ;
    xP[1 + 3 * i] -= grid.minY;
    xP[1 + 3 * i] -= grid.minZ;
    if (xP[3 * i] < 0 || xP[3 * i] > grid.Lx || 
        xP[1 + 3 * i] < 0 || xP[1 + 3 * i] > grid.Ly ||
        xP[2 + 3 * i] < 0 || xP[2 + 3 * i] > grid.Lz)
    {
      exitErr("Particles are outside of the unit cell.");
    }
  }
  wfxP_max = *std::max_element(wfxP, wfxP + nP); 
  wfyP_max = *std::max_element(wfyP, wfyP + nP); 
  wfzP_max = *std::max_element(wfzP, wfzP + nP); 

  grid.Nxeff = grid.Nx + wfxP_max + (wfxP_max % 2);
  grid.Nyeff = grid.Ny + wfyP_max + (wfyP_max % 2);
  grid.Nzeff = grid.Nz + wfzP_max + (wfzP_max % 2);  

  unsigned int N2 = grid.Nxeff * grid.Nyeff, N3 = N2 * grid.Nzeff;
  if (not grid.fG_unwrap) grid.fG_unwrap = (double*) fftw_malloc(N3 * grid.dof * sizeof(double)); 
  if (not grid.firstn) grid.firstn = (int*) fftw_malloc(N2 * sizeof(int));
  if (not grid.number) grid.number = (unsigned int*) fftw_malloc(N2 * sizeof(unsigned int));  
  if (not grid.nextn) grid.nextn = (int*) fftw_malloc(nP * sizeof(int));
  if (not xunwrap) xunwrap = (double *) fftw_malloc(wfxP_max * nP * sizeof(double));
  if (not yunwrap) yunwrap = (double *) fftw_malloc(wfyP_max * nP * sizeof(double));
  if (not zunwrap) zunwrap = (double *) fftw_malloc(wfzP_max * nP * sizeof(double));
  if (not zoffset) zoffset = (unsigned int *) fftw_malloc(nP * sizeof(unsigned int));
  std::fill(xunwrap, xunwrap + wfxP_max * nP, 0);
  std::fill(yunwrap, yunwrap + wfyP_max * nP, 0);
  std::fill(zunwrap, zunwrap + wfzP_max * nP, 0);
  
  unsigned int* xclose = (unsigned int*) fftw_malloc(nP * sizeof(unsigned int));
  unsigned int* yclose = (unsigned int*) fftw_malloc(nP * sizeof(unsigned int));
  

 
  #pragma omp parallel for num_threads(n_threads)
  for (unsigned int i = 0; i < nP; ++i)
  {
    unsigned short wx = wfxP[i];
    unsigned short wy = wfyP[i];
    unsigned short wz = wfzP[i];
    const int evenx = -1 * (wx % 2) + 1, eveny = -1 * (wy % 2) + 1;
    const int evenz = -1 * (wz % 2) + 1;
    xclose[i] = (int) (xP[3 * i] / grid.hx);
    yclose[i] = (int) (xP[1 + 3 * i] / grid.hy);
    unsigned int zclose = (int) (xP[2 + 3 * i] / grid.hz);
    xclose[i] += ((wx % 2) && (xP[3 * i] / grid.hx - xclose[i] > 1.0 / 2.0) ? 1 : 0);
    yclose[i] += ((wy % 2) && (xP[1 + 3 * i] / grid.hy - yclose[i] > 1.0 / 2.0) ? 1 : 0);
    zclose    += ((wz % 2) && (xP[2 + 3 * i] / grid.hz - zclose    > 1.0 / 2.0) ? 1 : 0);
    for (unsigned int j = 0; j < wx; ++j)
    {
      xunwrap[j + i * wfxP_max] = ((double) xclose[i] + j - wx / 2 + evenx) * grid.hx - xP[3 * i];
      if (xunwrap[j + i * wfxP_max] > threshP[i] || xunwrap[j + i * wfxP_max] < -threshP[i])
      {
        xunwrap[j + i * wfxP_max] = threshP[i];
      }
    } 
    for (unsigned int j = 0; j < wy; ++j)
    {
      yunwrap[j + i * wfyP_max] = ((double) yclose[i] + j - wy / 2 + eveny) * grid.hy - xP[1 + 3 * i];
      if (threshP[i] < yunwrap[j + i * wfyP_max] || yunwrap[j + i * wfyP_max] < -threshP[i])
      {
        yunwrap[j + i * wfyP_max] = threshP[i];
      }
    } 
    for (unsigned int j = 0; j < wz; ++j)
    {
      zunwrap[j + i * wfzP_max] = ((double) zclose + j - wz / 2 + evenz) * grid.hz - xP[2 + 3 * i];
      if (threshP[i] < zunwrap[j + i * wfzP_max] || zunwrap[j + i * wfzP_max] < -threshP[i] ) 
      {
        zunwrap[j + i * wfzP_max] = threshP[i];
      }
    }
    zoffset[i] = wx * wy * (zclose - wz / 2 + evenz + (grid.Nzeff - grid.Nz) / 2);    
    grid.nextn[i] = -1;
  }
  std::fill(grid.firstn, grid.firstn + N2, -1);
  std::fill(grid.number, grid.number + N2, 0);
  int ind, indn;
  for (unsigned int i = 0; i < nP; ++i) 
  {
    ind = (yclose[i] + (grid.Nyeff-grid.Ny) / 2) + (xclose[i] + (grid.Nxeff-grid.Nx) / 2) * grid.Nyeff;
    if (grid.firstn[ind] < 0) {grid.firstn[ind] = i;}
    else
    {
      indn = grid.firstn[ind];
      while (grid.nextn[indn] >= 0)
      {
        indn = grid.nextn[indn];
      }
      grid.nextn[indn] = i;
    }
    grid.number[ind] += 1;
  }
  if (xclose) {fftw_free(xclose); xclose = 0;}
  if (yclose) {fftw_free(yclose); yclose = 0;}
}

void ParticleList::locateOnGridNonUnifZ(Grid& grid)
{
  // get widths on effective uniform grid
  #pragma omp parallel for num_threads(n_threads)
  for (unsigned int i = 0; i < nP; ++i)
  {
    wfxP[i] = std::round(2 * threshP[i] / grid.hx);
    wfyP[i] = std::round(2 * threshP[i] / grid.hy);
    // shift to [0,L]x..
    xP[3 * i] -= grid.minX; 
    xP[1 + 3 * i] -= grid.minY;
    xP[2 + 3 * i] -= grid.minZ;
    if (xP[3 * i] < 0 || xP[3 * i] > grid.Lx || 
        xP[1 + 3 * i] < 0 || xP[1 + 3 * i] > grid.Ly ||
        xP[2 + 3 * i] < 0 || xP[2 + 3 * i] > grid.Lz)
    {
      exitErr("Particles are outside of the unit cell.");
    }
  }
  wfxP_max = *std::max_element(wfxP, wfxP + nP);
  wfyP_max = *std::max_element(wfyP, wfyP + nP);
  grid.Nxeff = grid.Nx + wfxP_max + (wfxP_max % 2);
  grid.Nyeff = grid.Ny + wfyP_max + (wfyP_max % 2);
  double threshP_max = *std::max_element(threshP, threshP + nP);
  
  // define extended z grid
  ext_down = 0; ext_up = 0;
  if (this->indl) this->indl = (unsigned short*) fftw_malloc(nP * sizeof(unsigned short));
  if (this->indr) this->indr = (unsigned short*) fftw_malloc(nP * sizeof(unsigned short));
  unsigned int i = 1;
  while (grid.zG[0] - grid.zG[i] <= threshP_max) {ext_up += 1; i += 1;} 
  i = grid.Nz - 2;
  while (grid.zG[i] - grid.zG[grid.Nz - 1] <= threshP_max) {ext_down += 1; i -= 1;}
  grid.Nzeff = grid.Nz + ext_up + ext_down;
  if (not grid.zG_ext)  grid.zG_ext = (double*) fftw_malloc(grid.Nzeff * sizeof(double));
  if (not grid.zG_ext_wts) grid.zG_ext_wts = (double*) fftw_malloc(grid.Nzeff * sizeof(double));
  // copy z grid
  for (unsigned int i = ext_up; i < grid.Nzeff - ext_down; ++i)
  {
    grid.zG_ext[i] = grid.zG[i - ext_up];
    grid.zG_ext_wts[i] = grid.zG_wts[i - ext_up];
  } 
  // define extended grid for z > b
  unsigned int j = 0;
  for (unsigned int i = ext_up; i > 0; --i) 
  {
    grid.zG_ext[j] = 2.0 * grid.zG[0] - grid.zG[i]; 
    grid.zG_ext_wts[j] = grid.zG_wts[i];
    j += 1;
  }
  // define extended grid for z < a=0
  j = grid.Nzeff - ext_down;
  for (unsigned int i = grid.Nz - 2; i > grid.Nz - 2 - ext_down; --i)
  { 
    grid.zG_ext[j] = -1.0 * grid.zG[i]; 
    grid.zG_ext_wts[j] = grid.zG_wts[i]; 
    j += 1; 
  } 
  // find wz for each particle
  for (unsigned int i = 0; i < nP; ++i)
  {
    // find index of z grid pt w/i alpha below 
    auto high = std::lower_bound(&grid.zG_ext[0], &grid.zG_ext[0] + grid.Nzeff, \
                                 xP[2 + 3 * i] - threshP[i], std::greater<double>());
    auto low = std::lower_bound(&grid.zG_ext[0], &grid.zG_ext[0] + grid.Nzeff, \
                                xP[2 + 3 * i] + threshP[i], std::greater<double>());
    indl[i] = low - &grid.zG_ext[0];  
    indr[i] = high - &grid.zG_ext[0];
    if (indr[i] == grid.Nzeff) {indr[i] -= 1;}
    else if (xP[2 + 3 * i] - threshP[i] > grid.zG_ext[indr[i]]) {indr[i] -= 1;}
    wfzP[i] = indr[i] - indl[i] + 1; 
  } 
  wfzP_max = *std::max_element(wfzP, wfzP + nP); 
  unsigned int N2 = grid.Nxeff * grid.Nyeff, N3 = N2 * grid.Nzeff;
  if (not grid.fG_unwrap) grid.fG_unwrap = (double*) fftw_malloc(N3 * grid.dof * sizeof(double)); 
  if (not grid.firstn) grid.firstn = (int*) fftw_malloc(N2 * sizeof(int));
  if (not grid.number) grid.number = (unsigned int*) fftw_malloc(N2 * sizeof(unsigned int));  
  if (not grid.nextn) grid.nextn = (int*) fftw_malloc(nP * sizeof(int));

  if (not xunwrap) xunwrap = (double*) fftw_malloc(wfxP_max * nP * sizeof(double));
  if (not yunwrap) yunwrap = (double*) fftw_malloc(wfyP_max * nP * sizeof(double));
  if (not zunwrap) zunwrap = (double*) fftw_malloc(wfzP_max * nP * sizeof(double));
  if (not pt_wts) pt_wts = (double*) fftw_malloc(wfzP_max * nP * sizeof(double));  
  if (not zoffset) zoffset = (unsigned int *) fftw_malloc(nP * sizeof(unsigned int));
  std::fill(xunwrap, xunwrap + wfxP_max * nP, 0);
  std::fill(yunwrap, yunwrap + wfyP_max * nP, 0);
  std::fill(zunwrap, zunwrap + wfzP_max * nP, 0);

  std::fill(pt_wts, pt_wts + wfzP_max * nP, 0);
  unsigned int* xclose = (unsigned int*) fftw_malloc(nP * sizeof(unsigned int));
  unsigned int* yclose = (unsigned int*) fftw_malloc(nP * sizeof(unsigned int));

  #pragma omp parallel for num_threads(n_threads)
  for (unsigned int i = 0; i < nP; ++i)
  {
    unsigned short wx = wfxP[i];
    unsigned short wy = wfyP[i];
    unsigned short wz = wfzP[i];
    const int evenx = -1 * (wx % 2) + 1, eveny = -1 * (wy % 2) + 1;
    xclose[i] = (int) (xP[3 * i] / grid.hx);
    yclose[i] = (int) (xP[1 + 3 * i] / grid.hy);
    xclose[i] += ((wx % 2) && (xP[3 * i] / grid.hx - xclose[i] > 1.0 / 2.0) ? 1 : 0);
    yclose[i] += ((wy % 2) && (xP[1 + 3 * i] / grid.hy - yclose[i] > 1.0 / 2.0) ? 1 : 0);
    for (unsigned int j = 0; j < wx; ++j)
    {
      xunwrap[j + i * wfxP_max] = ((double) xclose[i] + j - wx / 2 + evenx) * grid.hx - xP[3 * i];
      if (xunwrap[j + i * wfxP_max] > threshP[i] || xunwrap[j + i * wfxP_max] < -threshP[i])
      {
        xunwrap[j + i * wfxP_max] = threshP[i];
      }
    } 
    for (unsigned int j = 0; j < wy; ++j)
    {
      yunwrap[j + i * wfyP_max] = ((double) yclose[i] + j - wy / 2 + eveny) * grid.hy - xP[1 + 3 * i];
      if (threshP[i] < yunwrap[j + i * wfyP_max] || yunwrap[j + i * wfyP_max] < -threshP[i])
      {
        yunwrap[j + i * wfyP_max] = threshP[i] ;
      }
    }
    unsigned int k = 0;
    for (unsigned int j = indl[i]; j <= indr[i]; ++j)
    {
      zunwrap[k + i * wfzP_max] = grid.zG_ext[j] - xP[2 + 3 * i]; 
      if (threshP[i] < zunwrap[k + i * wfzP_max] || zunwrap[k + i * wfzP_max] < -threshP[i] ) 
      {
        zunwrap[k + i * wfzP_max] = threshP[i];
      }
      pt_wts[k + i * wfzP_max] = grid.hx * grid.hy * grid.zG_ext_wts[j];
      k += 1;
    }
    zoffset[i] = wx * wy * indl[i];  
    grid.nextn[i] = -1;
  }
  std::fill(grid.firstn, grid.firstn + N2, -1);
  std::fill(grid.number, grid.number + N2, 0);

  int ind, indn;
  for (unsigned int i = 0; i < nP; ++i) 
  {
    ind = (yclose[i] + (grid.Nyeff-grid.Ny) / 2) + (xclose[i] + (grid.Nxeff-grid.Nx) / 2) * grid.Nyeff;
    if (grid.firstn[ind] < 0) {grid.firstn[ind] = i;}
    else
    {
      indn = grid.firstn[ind];
      while (grid.nextn[indn] >= 0)
      {
        indn = grid.nextn[indn];
      }
      grid.nextn[indn] = i;
    }
    grid.number[ind] += 1;
  }
  if (xclose) {fftw_free(xclose); xclose = 0;}
  if (yclose) {fftw_free(yclose); yclose = 0;}
}

void ParticleList::update(double* xP_new, Grid& grid)
{
  #pragma omp parallel num_threads(n_threads)
  for (unsigned int i = 0; i < nP; ++i)
  {
    // shift to [0,L]x..
    xP_new[3 * i] -= grid.minX; 
    xP_new[1 + 3 * i] -= grid.minY;
    xP_new[2 + 3 * i] -= grid.minZ;
    int img;
    // periodically wrap
    if (xP_new[3 * i] < 0)  
    {
      img = 1 - ((int) (xP_new[3 * i] / grid.Lx));
      xP_new[3 * i] += img * grid.Lx;
    }
    else if (xP_new[3 * i] > grid.Lx) 
    {
      img = (int) (xP_new[3 * i] / grid.Lx);
      xP_new[3 * i] -= img * grid.Lx;
    }
    if (xP_new[1 + 3 * i] < 0) 
    {
      img = 1 - ((int) (xP_new[1 + 3 * i] / grid.Ly)); 
      xP_new[1 + 3 * i] += img * grid.Ly;
    }
    else if  (xP_new[1 + 3 * i] > grid.Ly) 
    {
      img = (int) (xP_new[1 + 3 * i] / grid.Ly);
      xP_new[1 + 3 * i] -= img * grid.Ly;
    }
    if (grid.isperiodic[2])
    {
      if (xP_new[2 + 3 * i] < 0) 
      {
        img = 1 - ((int) (xP_new[2 + 3 * i] / grid.Lz));
        xP_new[2 + 3 * i] += img * grid.Lz;
      }
      else if  (xP_new[2 + 3 * i] > grid.Lz) 
      {
        img = (int) (xP_new[2 + 3 * i] / grid.Lz);
        xP_new[2 + 3 * i] -= img * grid.Lz;
      }
    }
    // abort if new pos is out of z bounds for DP
    else if (xP_new[2 + 3 * i] < 0 || xP_new[2 + 3 * i] > grid.Lz) 
    {
      #pragma omp critical
      {
        exitErr("Position is outside of z extent");
      }
    }
    xP[3 * i] = xP_new[3 * i];
    xP[1 + 3 * i] = xP_new[1 + 3 * i];
    xP[2 + 3 * i] = xP_new[2 + 3 * i];
    // update z kernel coords
    if (grid.unifZ) // z uniform
    {
      unsigned short wx = wfxP[i];
      unsigned short wy = wfyP[i];
      unsigned short wz = wfzP[i];
      const int evenz = -1 * (wz % 2) + 1;
      unsigned int zclose = (int) (xP[2 + 3 * i] / grid.hz);
      zclose += ((wz % 2) && (xP[2 + 3 * i] / grid.hz - zclose    > 1.0 / 2.0) ? 1 : 0);
      for (unsigned int j = 0; j < wz; ++j)
      {
        zunwrap[j + i * wfzP_max] = ((double) zclose + j - wz / 2 + evenz) * grid.hz - xP[2 + 3 * i];
        if (threshP[i] < zunwrap[j + i * wfzP_max] || zunwrap[j + i * wfzP_max] < -threshP[i] ) 
        {
          zunwrap[j + i * wfzP_max] = threshP[i];
        }
      }
      zoffset[i] = wx * wy * (zclose - wz / 2 + evenz + (grid.Nzeff - grid.Nz) / 2);    
    }
    else // recompute wz max for non-uniform z
    {
      // find index of z grid pt w/i alpha below 
      auto high = std::lower_bound(&grid.zG_ext[0], &grid.zG_ext[0] + grid.Nzeff, \
                                   xP[2 + 3 * i] - threshP[i], std::greater<double>());
      auto low = std::lower_bound(&grid.zG_ext[0], &grid.zG_ext[0] + grid.Nzeff, \
                                  xP[2 + 3 * i] + threshP[i], std::greater<double>());
      indl[i] = low - &grid.zG_ext[0];  
      indr[i] = high - &grid.zG_ext[0];
      if (indr[i] == grid.Nzeff) {indr[i] -= 1;}
      else if (xP[2 + 3 * i] - threshP[i] > grid.zG_ext[indr[i]]) {indr[i] -= 1;}
      wfzP[i] = indr[i] - indl[i] + 1; 
    }
  }
  if (not grid.unifZ) // update z kernel coords for z non-uniform
  {
    unsigned short wfzP_max_old = wfzP_max;
    wfzP_max = *std::max_element(wfzP, wfzP + nP);
    if (wfzP_max != wfzP_max_old)
    {
      fftw_free(zunwrap); fftw_free(pt_wts);
      zunwrap = (double*) fftw_malloc(wfzP_max * nP * sizeof(double));
      pt_wts = (double*) fftw_malloc(wfzP_max * nP * sizeof(double));  
    } 
    #pragma omp parallel for num_threads(n_threads)
    for (unsigned int i = 0; i < nP; ++i)
    {
      unsigned short wx = wfxP[i];
      unsigned short wy = wfyP[i];
      unsigned int k = 0;
      for (unsigned int j = indl[i]; j <= indr[i]; ++j)
      {
        zunwrap[k + i * wfzP_max] = grid.zG_ext[j] - xP[2 + 3 * i]; 
        if (threshP[i] < zunwrap[k + i * wfzP_max] || zunwrap[k + i * wfzP_max] < -threshP[i] ) 
        {
          zunwrap[k + i * wfzP_max] = threshP[i];
        }
        pt_wts[k + i * wfzP_max] = grid.hx * grid.hy * grid.zG_ext_wts[j];
        k += 1;
      }
      zoffset[i] = wx * wy * indl[i];  
    }
  }
  // loop over unique alphas
  for (const double& alphaf : unique_threshP)
  {
    const unsigned short wx = std::round(2 * alphaf / grid.hx);
    const unsigned short wy = std::round(2 * alphaf / grid.hy);
    const int evenx = -1 * (wx % 2) + 1, eveny = -1 * (wy % 2) + 1;
    // loop over w^2 groups of columns
    for (unsigned int izero = 0; izero < wx; ++izero)
    {
      for (unsigned int jzero = 0; jzero < wy; ++jzero)
      {
        // Uncomment pragma if enforcing that particles
        // can move at most to a neighboring column per step

        //#pragma omp parallel for collapse(2)
        for (unsigned int ii = izero; ii < grid.Nxeff; ii += wx)
        {
          for (unsigned int jj = jzero; jj < grid.Nyeff; jj += wy)
          { 
            // column index
            int col = jj + ii * grid.Nyeff;
            // trailing ptr, index of first particle in col
            int nprev = -1, n = grid.firstn[col];
            // if there is a particle
            while (n > -1)
            {
              // if it has matching alpha
              if (threshP[n] == alphaf) 
              {
                // check if particle n has moved out of the column
                int xclose = (int) (xP_new[3 * n] / grid.hx);
                int yclose = (int) (xP_new[1 + 3 * n] / grid.hy);
                xclose += ((wx % 2) && (xP_new[3 * n] / grid.hx - xclose > 1.0 / 2.0) ? 1 : 0);
                yclose += ((wy % 2) && (xP_new[1 + 3 * n] / grid.hy - yclose > 1.0 / 2.0) ? 1 : 0);
                int col_new = (yclose + (grid.Nyeff-grid.Ny) / 2) + (xclose + (grid.Nxeff-grid.Nx) / 2) * grid.Nyeff;
                // if the particle has moved
                if (col_new != col)
                {
                  // index of next pt in column
                  int nnext = grid.nextn[n];
                  // update add particle n to col_new of search struct
                  grid.nextn[n] = grid.firstn[col_new];
                  grid.firstn[col_new] = n;
                  grid.number[col_new] += 1;
                  // delete particle n from col of search struct
                  if (nprev == -1){grid.firstn[col] = nnext;}
                  else {grid.nextn[nprev] = nnext;}
                  grid.number[col] -= 1;
                  n = nnext;
                  // update kernel coords
                  for (unsigned int j = 0; j < wx; ++j)
                  {
                    xunwrap[j + n * wfxP_max] = ((double) xclose + j - wx / 2 + evenx) * grid.hx - xP[3 * n];
                    if (xunwrap[j + n * wfxP_max] > threshP[n] || xunwrap[j + n * wfxP_max] < -threshP[n])
                    {
                      xunwrap[j + n * wfxP_max] = threshP[n];
                    }
                  } 
                  for (unsigned int j = 0; j < wy; ++j)
                  {
                    yunwrap[j + n * wfyP_max] = ((double) yclose + j - wy / 2 + eveny) * grid.hy - xP[1 + 3 * n];
                    if (threshP[n] < yunwrap[j + n * wfyP_max] || yunwrap[j + n * wfyP_max] < -threshP[n])
                    {
                      yunwrap[j + n * wfyP_max] = threshP[n] ;
                    }
                  }
                }
                // if the particle hasn't moved
                else
                {
                  // update trailing pointer and set n to next particle in col
                  nprev = n; n = grid.nextn[n];
                }
              }
              // if the particle does not have matching alpha
              else
              {
                // update trailing pointer and set n to next particle in col
                nprev = n; n = grid.nextn[n]; 
              }
            }
          }
        }
      }
    }
  }
}

/* write current state of ParticleList to ostream */
void ParticleList::writeParticles(std::ostream& outputStream) const
{
  if (this->validState() && outputStream.good()) 
  { 
    for (unsigned int i = 0; i < nP; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        outputStream << std::setprecision(16) << xP[j + i * 3] << " ";
      }
      for (unsigned int j = 0; j < this->dof; ++j)
      {
        outputStream << fP[j + i * dof] << " ";
      }
      outputStream << wfP[i] << " " << betafP[i] << " ";
      if (this->normalized) {outputStream << std::setprecision(16) << normfP[i] << " ";}
      outputStream << std::endl;
    }
  }
  else
  {
    exitErr("Unable to write particles to output stream.");
  }
}

/* write current state of ParticleList to file */
void ParticleList::writeParticles(const char* fname) const
{
  std::ofstream file; file.open(fname);
  writeParticles(file); file.close();
}

bool ParticleList::validState() const
{
  try
  {
    if (not (xP && fP && wfP && betafP && normfP && radP && 
             cwfP && alphafP && wfxP && wfyP && wfzP)) 
    {
      throw Exception("Pointer(s) is null", __func__, __FILE__,__LINE__ );
    } 
    if (not dof)
    {
      throw Exception("Degrees of freedom for data on particles must be specified.",
                      __func__, __FILE__, __LINE__);
    }
  } 
  catch (Exception& e)
  {
    e.getErr();
    return false;
  }
  return true;
}

void ParticleList::cleanup()
{
  if (this->validState())
  {
    if (xP) {fftw_free(xP); xP = 0;}
    if (fP) {fftw_free(fP); fP = 0;}
    if (betafP) {fftw_free(betafP); betafP = 0;}
    if (radP) {fftw_free(radP); radP = 0;}
    if (cwfP) {fftw_free(cwfP); cwfP = 0;}
    if (alphafP) {fftw_free(alphafP); alphafP = 0;}  
    if (normfP) {fftw_free(normfP); normfP = 0;}  
    if (wfP) {fftw_free(wfP); wfP = 0;} 
    if (wfxP) {fftw_free(wfxP); wfxP = 0;} 
    if (wfyP) {fftw_free(wfyP); wfyP = 0;} 
    if (wfzP) {fftw_free(wfzP); wfzP = 0;} 
    if (xunwrap) {fftw_free(xunwrap); xunwrap = 0;}
    if (yunwrap) {fftw_free(yunwrap); yunwrap = 0;}
    if (zunwrap) {fftw_free(zunwrap); zunwrap = 0;}
    if (zoffset) {fftw_free(zoffset); zoffset = 0;}
    if (pt_wts) {fftw_free(pt_wts); pt_wts = 0;}
    if (threshP) {fftw_free(threshP); threshP = 0;}
    if (indl) {fftw_free(indl); indl = 0;} 
    if (indr) {fftw_free(indr); indr = 0;} 
  }
}
