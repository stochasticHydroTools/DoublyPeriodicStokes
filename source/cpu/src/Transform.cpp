#include "Transform.h"
#include<omp.h>
#include<sstream>
#include<sys/stat.h>
#include<sys/types.h>

Transform::Transform() : in_real(0),in_imag(0),out_real(0),out_imag(0),
                         out_real_tr(0), out_imag_tr(0), Nx(0),Ny(0),
                         Nz(0),dof(0),rank(3),pF(0),pB(0),type(unset) {}

void Transform::init(const unsigned int Nx, const unsigned int Ny,
                     const unsigned int _Nz, const unsigned int dof,
                     const unsigned int n_threads, transformType type)
{
  this->Nx = Nx; this->Ny = Ny; this->dof = dof; this->type = type;

  if (type == fourier) {this->Nz = _Nz;}
  if (type == fourierCheb) {this->Nz = 2 * _Nz - 2; }

  // num threads for fftw
  this->n_threads = n_threads; fftw_plan_with_nthreads(n_threads);
  // configure dimensions of problem
  configDims();
  // allocate input arrs
  in_real = (double*) fftw_malloc((1 + 2 * (this->Nz * Ny * Nx * dof)) * sizeof(double));
  // complex part is next chunk in memory
  // this is necessary to use fftw_wisdom
  in_imag = in_real + this->Nz * Ny * Nx * dof + ((this->Nz * Ny * Nx * dof) % 2);
  if (!in_real) {exitErr("alloc failed in Transform");}
  // alias out to in for in-place transform
  out_real = in_real; out_imag = in_imag;
  // create forward and backward plans by 
  // reading from disk or create from scratch and save
  #ifdef ENABLE_WISDOM
  struct stat st = {0};
  if (stat("fftw_wisdom", &st) == -1)
  {
    mkdir("fftw_wisdom", 0777);
  }
  std::stringstream ss;
  ss << "fftw_wisdom/fftw_forward_wisdom_Nx" << this->Nx 
     << "_Ny" << this->Ny << "_Nz" << this->Nz << "_dof" << this->dof << "_nthr" << n_threads;
    #ifdef USE_FFTW_PATIENT
    ss << "_" << "patient";
    #else
    ss << "_" << "measure";
    #endif
  std::string fname1 = ss.str(); ss.str("");
  ss << "fftw_wisdom/fftw_backward_wisdom_Nx" << this->Nx 
     << "_Ny" << this->Ny << "_Nz" << this->Nz << "_dof" << this->dof << "_nthr" << n_threads;
    #ifdef USE_FFTW_PATIENT
    ss << "_" << "patient";
    #else
    ss << "_" << "measure";
    #endif
  std::string fname2 = ss.str(); ss.str("");
  int has_wisdom1 = fftw_import_wisdom_from_filename(fname1.c_str());
  if (has_wisdom1)
  {
    pF = fftw_plan_guru_split_dft(rank, dims, 1, howmany_dims, in_real, in_imag, out_real, out_imag, FFTW_WISDOM_ONLY);
  }
  else
  {
    #ifdef USE_FFTW_PATIENT
    pF = fftw_plan_guru_split_dft(rank, dims, 1, howmany_dims, in_real, in_imag, out_real, out_imag, FFTW_PATIENT);
    #else
    pF = fftw_plan_guru_split_dft(rank, dims, 1, howmany_dims, in_real, in_imag, out_real, out_imag, FFTW_MEASURE);
    #endif
  }
  int has_wisdom2 = fftw_import_wisdom_from_filename(fname2.c_str());
  if (has_wisdom2)
  {
    pB = fftw_plan_guru_split_dft(rank, dims, 1, howmany_dims, out_imag, out_real, in_imag, in_real, FFTW_WISDOM_ONLY);
  }
  else
  {
    #ifdef USE_FFTW_PATIENT
    pB = fftw_plan_guru_split_dft(rank, dims, 1, howmany_dims, out_imag, out_real, in_imag, in_real, FFTW_PATIENT);
    #else
    pB = fftw_plan_guru_split_dft(rank, dims, 1, howmany_dims, out_imag, out_real, in_imag, in_real, FFTW_MEASURE);
    #endif
  }
  if (!pF | !pB) {exitErr("FFTW planning failed");}
  if ((!has_wisdom1 && !fftw_export_wisdom_to_filename(fname1.c_str())) || 
      (!has_wisdom2 && !fftw_export_wisdom_to_filename(fname2.c_str())))
    {exitErr("FFTW could not save wisdom to disk");}
  #else
  pF = fftw_plan_guru_split_dft(rank, dims, 1, howmany_dims, in_real, in_imag, out_real, out_imag, FFTW_ESTIMATE);
  pB = fftw_plan_guru_split_dft(rank, dims, 1, howmany_dims, out_imag, out_real, in_imag, in_real, FFTW_ESTIMATE);
  if (!pF | !pB) {exitErr("FFTW planning failed");}
  #endif
  // truncated buffers for cheb transform
  if (type == fourierCheb)
  {
    out_real_tr = (double*) fftw_malloc(_Nz * Ny * Nx * dof * sizeof(double));
    out_imag_tr = (double*) fftw_malloc(_Nz * Ny * Nx * dof * sizeof(double));
  }
}

void Transform::setFData(const double* in_real)
{
  if (type == fourier)
  {
    unsigned int N = Nz * Ny * Nx * dof;
    #pragma omp parallel for
    for (unsigned int i = 0; i < N; ++i)
    {
      this->in_real[i] = in_real[i];
      this->in_imag[i] = 0;
    }
  }
  else if (type == fourierCheb)
  {
    unsigned int Nu = this->Nz;
    unsigned int _Nz = (unsigned int) (Nu + 2) / 2;
    unsigned int Ntot = Nu * Ny * Nx, Nyxd = dof * Ny * Nx;
    #pragma omp parallel for num_threads(n_threads)
    for (unsigned int ixy = 0; ixy < Ny * Nx; ++ixy)
    {
      double *this_inreal, *this_inimag;
      const double* inreal;
      this_inreal = &(this->in_real[at(0, 0, ixy, dof, Nu)]);
      this_inimag = &(this->in_imag[at(0, 0, ixy, dof, Nu)]);
      inreal = &(in_real[at(0, 0, ixy, dof, _Nz)]);
      #pragma omp simd
      for (unsigned int k = 0; k < dof * _Nz; ++k)
      {
        this_inreal[k] = inreal[k]; 
        this_inimag[k] = 0;
      }
      for (unsigned int k = _Nz; k < Nu; ++k)
      {
        this_inreal = &(this->in_real[at(0, k, ixy, dof, Nu)]);
        this_inimag = &(this->in_imag[at(0, k, ixy, dof, Nu)]);
        inreal = &(in_real[at(0, _Nz - 2 - (k - _Nz), ixy, dof, _Nz)]);
        #pragma omp simd
        for (unsigned int j = 0; j < dof; ++j)
        {
          this_inreal[j] = inreal[j]; 
          this_inimag[j] = 0;
        }
      }
    }
  }
}

void Transform::setBData(const double* in_real, const double* in_imag)
{
  if (type == fourier)
  {
    unsigned int N = this->Nz * this->Ny * this->Nx * this->dof;
    #pragma omp parallel for num_threads(n_threads)
    for (unsigned int i = 0; i < N; ++i)
    {
      this->in_real[i] = in_real[i];
      this->in_imag[i] = in_imag[i];
    }
  }
  else if (type == fourierCheb)
  {
    unsigned int Nu = this->Nz;
    unsigned int _Nz = (unsigned int) (Nu + 2) / 2;
    unsigned int Ntot = Nu * Ny * Nx;
    #pragma omp parallel for num_threads(n_threads)
    for (unsigned int ixy = 0; ixy < Nx * Ny; ++ixy)
    {
      double *this_inreal, *this_inimag;
      const double *inreal, *inimag;
      this_inreal = &(this->in_real[at(0, 0, ixy, dof, Nu)]);
      this_inimag = &(this->in_imag[at(0, 0, ixy, dof, Nu)]);
      inreal = &(in_real[at(0, 0, ixy, dof, _Nz)]);
      inimag = &(in_imag[at(0, 0, ixy, dof, _Nz)]);
      #pragma omp simd
      for (unsigned int j = 0; j < dof; ++j)
      {
        this_inreal[j] = Nu * inreal[j];
        this_inimag[j] = Nu * inimag[j];
      }
      #pragma omp simd
      for (unsigned int k = dof; k < dof * (_Nz - 1); ++k)
      { 
        this_inreal[k] = Nu * inreal[k] / 2; 
        this_inimag[k] = Nu * inimag[k] / 2; 
      }
      #pragma omp simd
      for (unsigned int k = dof * (_Nz - 1); k < dof * _Nz; ++k)
      {
        this_inreal[k] = Nu * inreal[k];
        this_inimag[k] = Nu * inimag[k];
      }
      for (unsigned int k = _Nz; k < Nu; ++k)
      {
        this_inreal = &(this->in_real[at(0, k, ixy, dof, Nu)]);
        this_inimag = &(this->in_imag[at(0, k, ixy, dof, Nu)]);
        inreal = &(in_real[at(0, _Nz - 2 - (k - _Nz), ixy, dof, _Nz)]);
        inimag = &(in_imag[at(0, _Nz - 2 - (k - _Nz), ixy, dof, _Nz)]);
        #pragma omp simd
        for (unsigned int j = 0; j < dof; ++j)
        {
          this_inreal[j] = Nu * inreal[j] / 2; 
          this_inimag[j] = Nu * inimag[j] / 2;
        }
      }
    }
  }
}

void Transform::fTransform()
{
  fftw_execute_split_dft(pF, in_real, in_imag, out_real, out_imag);
  if (type == fourier) return;
  else if (type == fourierCheb)
  {
    unsigned int Nu = this->Nz;
    unsigned int _Nz = (unsigned int) (Nu + 2) / 2;
    #pragma omp parallel for num_threads(n_threads)
    for (unsigned int ixy = 0; ixy < Nx * Ny; ++ixy)
    {
      const double *this_outreal, *this_outimag;
      const double *this_outreal1, *this_outimag1;
      double *outreal, *outimag;

      this_outreal = &(this->out_real[at(0, 0, ixy, dof, Nu)]);
      this_outimag = &(this->out_imag[at(0, 0, ixy, dof, Nu)]);
      outreal = &(out_real_tr[at(0, 0, ixy, dof, _Nz)]);
      outimag = &(out_imag_tr[at(0, 0, ixy, dof, _Nz)]);
      #pragma omp simd
      for (unsigned int j = 0; j < dof; ++j)
      {
        outreal[j] = this_outreal[j] / Nu;
        outimag[j] = this_outimag[j] / Nu;
      }
      
      for (unsigned int k = 1; k < _Nz - 1; ++k) 
      { 
        outreal = &(out_real_tr[at(0, k, ixy, dof, _Nz)]); 
        outimag = &(out_imag_tr[at(0, k, ixy, dof, _Nz)]); 
        this_outreal = &(this->out_real[at(0, k, ixy, dof, Nu)]);
        this_outimag = &(this->out_imag[at(0, k, ixy, dof, Nu)]);
        this_outreal1 = &(this->out_real[at(0, Nu - k, ixy, dof, Nu)]);
        this_outimag1 = &(this->out_imag[at(0, Nu - k, ixy, dof, Nu)]);
        #pragma omp simd
        for (unsigned int j = 0; j < dof; ++j)
        {
          outreal[j] = (this_outreal[j] + this_outreal1[j]) / Nu;
          outimag[j] = (this_outimag[j] + this_outimag1[j]) / Nu;
        }
      }
      this_outreal = &(this->out_real[at(0, 0, ixy, dof, Nu)]);
      this_outimag = &(this->out_imag[at(0, 0, ixy, dof, Nu)]);
      outreal = &(out_real_tr[at(0, 0, ixy, dof, _Nz)]);
      outimag = &(out_imag_tr[at(0, 0, ixy, dof, _Nz)]);
      #pragma omp simd
      for (unsigned int k = dof * (_Nz - 1); k < dof * _Nz; ++k)
      {
        outreal[k] = this_outreal[k] / Nu;
        outimag[k] = this_outimag[k] / Nu;
      }
    }
  }
}

void Transform::bTransform()
{
  fftw_execute_split_dft(pB, out_imag, out_real, in_imag, in_real);
  if (type == fourier)
  {
    unsigned int N = Nx * Ny * Nz;
    #pragma omp parallel for num_threads(n_threads)
    for (unsigned int i = 0; i < N * dof; ++i)
    {
      out_real[i] /= N; out_imag[i] /= N;
    }
  }
  else if (type == fourierCheb)
  { 
    unsigned int Nu = this->Nz;
    unsigned int _Nz = (unsigned int) (Nu + 2) / 2, dof = this->dof;
    unsigned int Ntot = Nu * Ny * Nx;
    #pragma omp parallel for num_threads(n_threads)
    for (unsigned int ixy = 0; ixy < Ny * Nx; ++ixy)
    {
      const double *this_outreal, *this_outimag;
      double *outreal, *outimag;
      this_outreal = &(this->out_real[at(0, 0, ixy, dof, Nu)]);
      this_outimag = &(this->out_imag[at(0, 0, ixy, dof, Nu)]);
      outreal = &(out_real_tr[at(0, 0, ixy, dof, _Nz)]);
      outimag = &(out_imag_tr[at(0, 0, ixy, dof, _Nz)]);
      #pragma omp simd
      for (unsigned int k = 0; k < dof * _Nz; ++k)
      {
        outreal[k] = this_outreal[k] / Ntot;
        outimag[k] = this_outimag[k] / Ntot;
      }
    }
  }
}

void Transform::configDims()
{
  // set up iodims - we store as (k,j,i,l), l = 0:dof 
  dims = (fftw_iodim*) fftw_malloc(rank * sizeof(fftw_iodim));    
  if (!dims) {exitErr("alloc failed in configDims for Transform");}
  // we want to do 1 fft for the entire dof x 3D array     
  howmany_dims = (fftw_iodim*) fftw_malloc(1 * sizeof(fftw_iodim));
  if (!howmany_dims) {exitErr("alloc failed in configDims for Transform");}
  // size of k
  dims[0].n = Ny;
  // stride for k
  dims[0].is = dof * Nz * Nx;
  dims[0].os = dof * Nz * Nx;
  // size of j
  dims[1].n = Nx;
  // stride for j
  dims[1].is = dof * Nz;
  dims[1].os = dof * Nz;
  // size of i
  dims[2].n = Nz;
  // stride for i
  dims[2].is = dof;
  dims[2].os = dof;
  
  // dof component vec field
  howmany_dims[0].n = dof;
  // stride of 1 b/w each component (interleaved)
  howmany_dims[0].is = 1;
  howmany_dims[0].os = 1;
}

void Transform::cleanup()
{
  // destroy plans
  if (pF) {fftw_destroy_plan(pF);}
  if (pB) {fftw_destroy_plan(pB);}
  // free memory
  if (in_real) {fftw_free(in_real); in_real = 0;}
  //if (in_imag) {fftw_free(in_imag); in_imag = 0;}
  if (out_real_tr) {fftw_free(out_real_tr); out_real_tr = 0;}
  if (out_imag_tr) {fftw_free(out_imag_tr); out_imag_tr = 0;}
  if (dims) {fftw_free(dims); dims = 0;}
  if (howmany_dims) {fftw_free(howmany_dims); howmany_dims = 0;}
  fftw_forget_wisdom();
  fftw_cleanup_threads();
}
