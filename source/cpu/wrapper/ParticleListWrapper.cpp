#include "ParticleList.h"

/* C wrapper for calling ParticleList methods from Python. Any functions
   defined here should also have their prototypes 
   and wrappers defined in ParticleList.py */
extern "C"
{
  /* make a ParticleList from external data */
  ParticleList* MakeParticles(const double* xP, const double* fP, const double* radP, 
                              const double* betafP, const double* cwfP, const unsigned short* wfP, 
                              const unsigned int nP, const unsigned int dof, const bool isDipole, const bool isCbeta)
  {
    ParticleList* particles = new ParticleList(xP, fP, radP, betafP, cwfP, wfP, nP, dof, isDipole, isCbeta);
    return particles;
  }

  /* make a ParticleList from external data but undetermined radii*/
  ParticleList* MakeParticles_norad(const double* xP, const double* fP, const double* betafP, 
                                    const double* alphafP, const unsigned short* wfP, 
                                    const unsigned int nP, const unsigned int dof, const bool isDipole)
  {
    ParticleList* particles = new ParticleList(xP, fP, betafP, alphafP, wfP, nP, dof, isDipole);
    return particles;
  }
  
  /* setup the particles on the grid (builds grid locators)*/
  void Setup(ParticleList* particles, Grid* grid)
  {
    particles->setup(*grid);  
  }

  /* set or get data on the particles */
  void SetData(ParticleList* particles, const double* _fP, unsigned int dof)
  {
    particles->setData(_fP, dof);
  }
  
  /* zero the data on particles */
  void ZeroData(ParticleList* particles) {particles->zeroData();}
  /* shallow copy data */
  double* GetData(ParticleList* particles) {return particles->fP;}
  /* deep copy data into provided compatible array */
  void CopyData(ParticleList* p, double* fP) {std::copy(p->fP, p->fP + p->dof * p->nP, fP);}
 
  double* getPoints(ParticleList* particles) {return particles->xP;}
  double* getRadii(ParticleList* particles) {return particles->radP;}
  double* getBetaf(ParticleList* particles) {return particles->betafP;}
  unsigned short* getWf(ParticleList* particles) {return particles->wfP;}
  unsigned short* getWfx(ParticleList* particles) {return particles->wfxP;}
  unsigned short* getWfy(ParticleList* particles) {return particles->wfyP;}
  unsigned short* getWfz(ParticleList* particles) {return particles->wfzP;}
  double* getNormf(ParticleList* particles) 
  {
    if (not particles->normalized) particles->normalizeKernels();
    return particles->normfP;
  }
  void Update(ParticleList* s, Grid* g, double* x_new) {s->update(x_new, *g);}
  void CleanParticles(ParticleList* s) {s->cleanup();}
  void DeleteParticles(ParticleList* s) {if(s) {delete s; s = 0;}}
  void WriteParticles(ParticleList* s, const char* fname) {s->writeParticles(fname);}
  void SetNumThreads(ParticleList* s, int num_threads) {s->n_threads =  num_threads;}
}
