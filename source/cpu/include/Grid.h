#ifndef GRID_H
#define GRID_H
#include<ostream>
#include"common.h"
/* Grid is an SoA describing the domain and its data

 * fG                     - forces on the grid
 * fG_unwrap              - forces on extended grid (used internally for BCs)
 * xG, yG, zG             - grids for each axis (sorted in inc or dec order) (see below)
 * Lx, Ly, Lz, hx, hy, hz - length and grid spacing in each dimension 
 *                        - if hx > 0, xG should be Null (same for y,z)
 *                        - if hx = 0, xG must be allocated (same for y,z)
 * Nxeff, Nyeff, Nzeff    - num points in each dimension for EXTENDED grid
 * isperiodic             - bool array indicating whether periodicity is on or off for each axis
 * has_bc                 - bool array indicating whether BCs for each dof are specified
 * firstn, nextn          - enables the lookup of particles in terms of columns of the grid 
                            for column ind, grid.firstn[ind] = i1 is the index of the first particle in the column
                            grid.nextn[i1] = i2 is the index of the next particle in the column, and so on.
 * BCs                    - enum for boundary conditions for each dof at the ends of each axis (dof x 6)
 * n_threads              - number of threads to use for any omp loop (default to 1, set externally with c wrapper)
*/ 


struct Grid
{
  double *fG, *fG_unwrap, *xG, *yG, *zG, *zG_wts; 
  double minX, minY, minZ;
  int *firstn, *nextn;
  unsigned int* number;
  unsigned int Nx, Ny, Nz, dof;
  double Lx, Ly, Lz;
  double hx, hy, hz;
  unsigned int Nxeff, Nyeff, Nzeff;
  bool has_locator, has_bc, unifZ;
  // bool array specifying if grid is periodic in direction i
  // and another bool to make sure this array is populated
  bool isperiodic[3], has_periodicity;
  // enum for boundary conditions for each dof at the ends of each axis (dof x 6)
  BC* BCs;
  // number of threads for routines using omp
  int n_threads;
    
  /* empty/null ctor */
  Grid();
  /* set up Grid based on what caller has provided */
  void setup();
  /* zero the extended grid */
  void zeroExtGrid();
  /* zero the extended grid ghost cells only*/
  void zeroExtGrid_ghost();
  /* zero interior grid */
  void zeroIntGrid();
  /* set wrapped grid data */
  void setData(double* data); 
  /* clean memory */
  void cleanup();
  /* check validity of current state */
  bool validState() const;
  
  /* write Grid data F to ostream */
  void writeGrid(std::ostream& outputStream) const;
  void writeGrid(const char* fname) const;
  void writeCoords(std::ostream& outputStream) const;
  void writeCoords(const char* fname) const;
};


#endif
