#include <math.h>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#define  NUMBER_IS_2_POW_K(x)   ((!((x)&((x)-1)))&&((x)>1))  // x is pow(2, k), k=1,2, ...
#define  FT_DIRECT        -1    // Direct transform.
#define  FT_INVERSE        1    // Inverse transform.




bool  FFT(double *Rdat, double *Idat, int N, int LogN, int Ft_Flag)
{
  if((Rdat == NULL) || (Idat == NULL))                  return false;
  if((N > 16384) || (N < 1))                            return false;
  if(!NUMBER_IS_2_POW_K(N))                             return false;
  if((LogN < 2) || (LogN > 14))                         return false;
  if((Ft_Flag != FT_DIRECT) && (Ft_Flag != FT_INVERSE)) return false;
 
  register int  i, j, n, k, io, ie, in, nn;
  double         ru, iu, rtp, itp, rtq, itq, rw, iw, sr;
 
  static const double Rcoef[14] =
  {  -1.0000000000000000F,  0.0000000000000000F,  0.7071067811865475F,
      0.9238795325112867F,  0.9807852804032304F,  0.9951847266721969F,
      0.9987954562051724F,  0.9996988186962042F,  0.9999247018391445F,
      0.9999811752826011F,  0.9999952938095761F,  0.9999988234517018F,
      0.9999997058628822F,  0.9999999264657178F
  };
  static const double Icoef[14] =
  {   0.0000000000000000F, -1.0000000000000000F, -0.7071067811865474F,
     -0.3826834323650897F, -0.1950903220161282F, -0.0980171403295606F,
     -0.0490676743274180F, -0.0245412285229122F, -0.0122715382857199F,
     -0.0061358846491544F, -0.0030679567629659F, -0.0015339801862847F,
     -0.0007669903187427F, -0.0003834951875714F
  };
 
  nn = N >> 1;
  ie = N;
  for(n=1; n<=LogN; n++)
  {
    rw = Rcoef[LogN - n];
    iw = Icoef[LogN - n];
    if(Ft_Flag == FT_INVERSE) iw = -iw;
    in = ie >> 1;
    ru = 1.0F;
    iu = 0.0F;
    for(j=0; j<in; j++)
    {
      for(i=j; i<N; i+=ie)
      {
        io       = i + in;
        rtp      = Rdat[i]  + Rdat[io];
        itp      = Idat[i]  + Idat[io];
        rtq      = Rdat[i]  - Rdat[io];
        itq      = Idat[i]  - Idat[io];
        Rdat[io] = rtq * ru - itq * iu;
        Idat[io] = itq * ru + rtq * iu;
        Rdat[i]  = rtp;
        Idat[i]  = itp;
      }
 
      sr = ru;
      ru = ru * rw - iu * iw;
      iu = iu * rw + sr * iw;
    }
 
    ie >>= 1;
  }
 
  for(j=i=1; i<N; i++)
  {
    if(i < j)
    {
      io       = i - 1;
      in       = j - 1;
      rtp      = Rdat[in];
      itp      = Idat[in];
      Rdat[in] = Rdat[io];
      Idat[in] = Idat[io];
      Rdat[io] = rtp;
      Idat[io] = itp;
    }
 
    k = nn;
 
    while(k < j)
    {
      j   = j - k;
      k >>= 1;
    }
 
    j = j + k;
  }
 
  
  rw=1.0/(sqrt((double)N));
 
  for(i=0; i<N; i++)
  {
    Rdat[i]*=rw;
    Idat[i]*=rw;
  }
 
  return true;
}
 

bool FFTX(double **Rdat, double **Idat, int N, int LogN, int Ft_Flag, int y)
{
  if((Rdat == NULL) || (Idat == NULL))                  return false;
  if((N > 16384) || (N < 1))                            return false;
  if(!NUMBER_IS_2_POW_K(N))                             return false;
  if((LogN < 2) || (LogN > 14))                         return false;
  if((Ft_Flag != FT_DIRECT) && (Ft_Flag != FT_INVERSE)) return false;
 
  register int  i, j, n, k, io, ie, in, nn;
  double         ru, iu, rtp, itp, rtq, itq, rw, iw, sr;
 
  static const double Rcoef[14] =
  {  -1.0000000000000000F,  0.0000000000000000F,  0.7071067811865475F,
      0.9238795325112867F,  0.9807852804032304F,  0.9951847266721969F,
      0.9987954562051724F,  0.9996988186962042F,  0.9999247018391445F,
      0.9999811752826011F,  0.9999952938095761F,  0.9999988234517018F,
      0.9999997058628822F,  0.9999999264657178F
  };
  static const double Icoef[14] =
  {   0.0000000000000000F, -1.0000000000000000F, -0.7071067811865474F,
     -0.3826834323650897F, -0.1950903220161282F, -0.0980171403295606F,
     -0.0490676743274180F, -0.0245412285229122F, -0.0122715382857199F,
     -0.0061358846491544F, -0.0030679567629659F, -0.0015339801862847F,
     -0.0007669903187427F, -0.0003834951875714F
  };


  nn = N >> 1;
  ie = N;
  for(n=1; n<=LogN; n++)
  {
    rw = Rcoef[LogN - n];
    iw = Icoef[LogN - n];
    if(Ft_Flag == FT_INVERSE) iw = -iw;
    in = ie >> 1;
    ru = 1.0F;
    iu = 0.0F;
    for(j=0; j<in; j++)
    {
      for(i=j; i<N; i+=ie)
      {
        io       = i + in;
        rtp      = Rdat[i][y]  + Rdat[io][y];
        itp      = Idat[i][y]  + Idat[io][y];
        rtq      = Rdat[i][y]  - Rdat[io][y];
        itq      = Idat[i][y]  - Idat[io][y];
        Rdat[io][y] = rtq * ru - itq * iu;
        Idat[io][y] = itq * ru + rtq * iu;
        Rdat[i][y]  = rtp;
        Idat[i][y]  = itp;
      }
 
      sr = ru;
      ru = ru * rw - iu * iw;
      iu = iu * rw + sr * iw;
    }
 
    ie >>= 1;
  }
 
  for(j=i=1; i<N; i++)
  {
    if(i < j)
    {
      io       = i - 1;
      in       = j - 1;
      rtp      = Rdat[in][y];
      itp      = Idat[in][y];
      Rdat[in][y] = Rdat[io][y];
      Idat[in][y] = Idat[io][y];
      Rdat[io][y] = rtp;
      Idat[io][y] = itp;
    }
 
    k = nn;
 
    while(k < j)
    {
      j   = j - k;
      k >>= 1;
    }
 
    j = j + k;
  }
 
  
  rw=1.0/(sqrt((double)N));
 
  for(i=0; i<N; i++)
  {
    Rdat[i][y]*=rw;
    Idat[i][y]*=rw;
  }
 
  return true;
}

bool  FFTY(double **Rdat, double **Idat, int N, int LogN, int Ft_Flag, int x)
{
  if((Rdat == NULL) || (Idat == NULL))                  return false;
  if((N > 16384) || (N < 1))                            return false;
  if(!NUMBER_IS_2_POW_K(N))                             return false;
  if((LogN < 2) || (LogN > 14))                         return false;
  if((Ft_Flag != FT_DIRECT) && (Ft_Flag != FT_INVERSE)) return false;
 
  register int  i, j, n, k, io, ie, in, nn;
  double         ru, iu, rtp, itp, rtq, itq, rw, iw, sr;
 
  static const double Rcoef[14] =
  {  -1.0000000000000000F,  0.0000000000000000F,  0.7071067811865475F,
      0.9238795325112867F,  0.9807852804032304F,  0.9951847266721969F,
      0.9987954562051724F,  0.9996988186962042F,  0.9999247018391445F,
      0.9999811752826011F,  0.9999952938095761F,  0.9999988234517018F,
      0.9999997058628822F,  0.9999999264657178F
  };
  static const double Icoef[14] =
  {   0.0000000000000000F, -1.0000000000000000F, -0.7071067811865474F,
     -0.3826834323650897F, -0.1950903220161282F, -0.0980171403295606F,
     -0.0490676743274180F, -0.0245412285229122F, -0.0122715382857199F,
     -0.0061358846491544F, -0.0030679567629659F, -0.0015339801862847F,
     -0.0007669903187427F, -0.0003834951875714F
  };
 
  nn = N >> 1;
  ie = N;
  for(n=1; n<=LogN; n++)
  {
    rw = Rcoef[LogN - n];
    iw = Icoef[LogN - n];
    if(Ft_Flag == FT_INVERSE) iw = -iw;
    in = ie >> 1;
    ru = 1.0F;
    iu = 0.0F;
    for(j=0; j<in; j++)
    {
      for(i=j; i<N; i+=ie)
      {
        io       = i + in;
        rtp      = Rdat[x][i]  + Rdat[x][io];
        itp      = Idat[x][i]  + Idat[x][io];
        rtq      = Rdat[x][i]  - Rdat[x][io];
        itq      = Idat[x][i]  - Idat[x][io];
        Rdat[x][io] = rtq * ru - itq * iu;
        Idat[x][io] = itq * ru + rtq * iu;
        Rdat[x][i]  = rtp;
        Idat[x][i]  = itp;
      }
 
      sr = ru;
      ru = ru * rw - iu * iw;
      iu = iu * rw + sr * iw;
    }
 
    ie >>= 1;
  }
 
  for(j=i=1; i<N; i++)
  {
    if(i < j)
    {
      io       = i - 1;
      in       = j - 1;
      rtp      = Rdat[x][in];
      itp      = Idat[x][in];
      Rdat[x][in] = Rdat[x][io];
      Idat[x][in] = Idat[x][io];
      Rdat[x][io] = rtp;
      Idat[x][io] = itp;
    }
 
    k = nn;
 
    while(k < j)
    {
      j   = j - k;
      k >>= 1;
    }
 
    j = j + k;
  }
 
  
  rw=1.0/(sqrt((double)N));
 
  for(i=0; i<N; i++)
  {
    Rdat[x][i]*=rw;
    Idat[x][i]*=rw;
  }
 
  return true;
}

bool FFT2D(double **Rdat, double **Idat, int N, int LogN, int Ft_Flag)
{
	if((Rdat == NULL) || (Idat == NULL))                  return false;
	if((N > 16384) || (N < 1))                            return false;
	if(!NUMBER_IS_2_POW_K(N))                             return false;
	if((LogN < 2) || (LogN > 14))                         return false;
	if((Ft_Flag != FT_DIRECT) && (Ft_Flag != FT_INVERSE)) return false;


	int i;

	for (i=0; i<N; i++)
		FFTX(Rdat, Idat, N, LogN, Ft_Flag, i);
	for (i=0; i<N; i++)
		FFTY(Rdat, Idat, N, LogN, Ft_Flag, i);


	return true;
}
