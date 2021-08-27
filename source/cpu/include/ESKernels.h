#ifndef _ESKERNEL_H_
#define _ESKERNEL_H_
#include<math.h>
#include<iostream>
#if !defined(TORQUE_X) && !defined(TORQUE_Y) && !defined(TORQUE_Z)
// ES kernel definition, where x_i <= alpha 

#pragma omp declare simd
inline double esKernel(const double x[3], const double beta, const double alpha)
{
  double alpha2 = alpha * alpha;
  return exp(beta * (sqrt(1 - x[0] * x[0] / alpha2) + \
                     sqrt(1 - x[1] * x[1] / alpha2) + \
                     sqrt(1 - x[2] * x[2] / alpha2) - 3));
}
#else

#if defined(TORQUE_X)
#define COORD 0
#elif defined(TORQUE_Y)
#define COORD 1
#elif defined(TORQUE_Z)
#define COORD 2
#endif

/* 
  deriv of es kernel phi'(z) = -beta * z * phi(z) / (alpha^2 * sqrt(1 - (z/alpha)^2))

  We must have z <= thresh = creal(csqrt(0.6e1 * cpow(((8 * beta * beta) + 
                             0.12e2 * csqrt((-12 * beta * beta + 81)) - 
                             0.108e3) * beta, 0.1e1 / 0.3e1) / beta + 
                             0.24e2 * beta * cpow(((8 * beta * beta) + 0.12e2 * 
                             csqrt((-12 * beta * beta + 81)) - 0.108e3) * beta, 
                             -0.1e1 / 0.3e1) + 0.12e2) * alpha / 0.6e1)
  so that we avoid the singularity at z = alpha. That is, 1 \approx thresh/alpha < 1
*/
#pragma omp declare simd
inline double esKernel(const double x[3], const double beta, const double alpha)
{
  double alpha2 = alpha *  alpha;
  return -beta * x[COORD] * \
          exp(beta * (sqrt(1 - x[0] * x[0] / alpha2) + \
                      sqrt(1 - x[1] * x[1] / alpha2) + \
                      sqrt(1 - x[2] * x[2] / alpha2) - 3)) / \
          (alpha2 * sqrt(1 - x[COORD] * x[COORD] / (alpha2)));
}
#endif

/* singleton version of kernel for computing normalizations for particles using 
   a given beta,alpha. */
inline double eskernel(const double x, const double beta, const double alpha)
{
  return exp(beta * (sqrt(1 - x * x / (alpha * alpha)) - 1));
}

#endif
