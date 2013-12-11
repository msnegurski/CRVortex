#ifndef __FFT_FUNCTIONS_H_
#define __FFT_FUNCTIONS_H_

#include <stdbool.h>

// x is pow(2, k), k=1,2, ...
#define  NUMBER_IS_2_POW_K(x)   ((!((x)&((x)-1)))&&((x)>1))

#define  FT_DIRECT  -1 // Direct transform.
#define  FT_INVERSE  1 // Inverse transform.

bool FFT(double *Rdat, double *Idat, int N, int LogN, int Ft_Flag);
bool FFTX(double **Rdat, double **Idat, int N, int LogN, int Ft_Flag, int y);
bool  FFTY(double **Rdat, double **Idat, int N, int LogN, int Ft_Flag, int x);
bool FFT2D(double **Rdat, double **Idat, int N, int LogN, int Ft_Flag);

#endif
