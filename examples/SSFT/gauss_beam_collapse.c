#include "fft_functions.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double gradx_sq(double **Rdat, double **Idat, int i, int j, double dx)
{
    return ((Rdat[i][j] - Rdat[i][j-1]) * (Rdat[i][j] - Rdat[i][j-1]) +
            (Idat[i][j] - Idat[i][j-1]) * (Idat[i][j] - Idat[i][j-1])) / (dx * dx);
}

double grady_sq(double **Rdat, double **Idat, int i, int j, double dx)
{
    return ((Rdat[i][j] - Rdat[i-1][j]) * (Rdat[i][j] - Rdat[i-1][j]) +
            (Idat[i][j] - Idat[i-1][j]) * (Idat[i][j] - Idat[i-1][j])) / (dx * dx);
}

double psi_sq(double **Rdat, double **Idat, int i, int j)
{
    return (Rdat[i][j] * Rdat[i][j] + Idat[i][j] * Idat[i][j]);
}

// Возвращает полную мощность луча
double power(double **Rdat, double **Idat, int N, double dx)
{
    double p = 0.0;
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            p += Rdat[i][j] * Rdat[i][j] * dx * dx;
            p += Idat[i][j] * Idat[i][j] * dx * dx;
        }
    }
    return p;
}

// Возвращает полную энергию луча
double energy(double **Rdat, double **Idat, int N, double dx)
{
    double e = 0.0;
    int i, j;
    for (i = 1; i < N; i++)
    {
        for (j = 1; j < N; j++)
        {
            e += gradx_sq(Rdat, Idat, i, j, dx) +
                 grady_sq(Rdat, Idat, i, j, dx) -
             0.5 * psi_sq(Rdat, Idat, i, j) * psi_sq(Rdat, Idat, i, j);
            //e+ = (gradx_sq(Rdat,Idat,i,j,dx)+grady_sq(Rdat,Idat,i,j,dx)) * dx * dx;
            //e+ = psi_sq(Rdat,Idat,i,j) * psi_sq(Rdat,Idat,i,j) * dx * dx;
        }
    }

    e *= dx * dx;
    return e;
}

void initiate(double **Rdat, double **Idat, int N,double dx,double ampl)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            // Начальные условия - гауссоида
            Rdat[i][j] = ampl * exp(-((i - N/2) * (i - N/2) +
                        (j - N/2) * (j - N/2)) * dx * dx/2);
            Idat[i][j] = 0.0;

            /*
            RePsi[i][j] = sin((i+j) * dx);
            ImPsi[i][j] = 0;
            */
        }
    }
}

// Шаг эволюции системы во времени, зависящий от нелинейности (кубической)
void cubic_evol(double **Rdat, double **Idat, int N, double dt)
{
    int i, j;
    double buf1, buf2, buf3;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            buf1 = Rdat[i][j];
            buf2 = Idat[i][j];
            buf3 = -psi_sq(Rdat, Idat, i, j) * dt;
            Rdat[i][j] = buf1 * cos(buf3) - buf2 * sin(buf3);
            Idat[i][j] = buf1 * sin(buf3) + buf2 * cos(buf3);
        }
    }
}

// Шаг эволюции системы во времени, зависящий от линейного члена
void linear_evol(double **Rdat, double **Idat, int N, int LogN, double dt, double dx)
{
    int i, j;
    double kx, ky;
    double p = 2 * 3.141592653589 / N;

    double t = p * p / dx / dx;
    double buf1, buf2, buf3;

    FFT2D(Rdat, Idat, N, LogN, -1);

    for (i = 0; i < N; i++)
    {
        for (j = 0; j  <N; j++)
        {
            if (i < N / 2)
            {
                kx = i;
            } else
            {
                kx = (N - i);
            }

            if (j < N / 2)
            {
                ky = j;
            } else
            {
                ky = (N - j);
            }

            buf1 = -(kx * kx+ky * ky) * t * dt;
            buf2 = Rdat[i][j];
            buf3 = Idat[i][j];

            Rdat[i][j] = buf2 * cos(buf1) - buf3 * sin(buf1);
            Idat[i][j] = buf2 * sin(buf1) + buf3 * cos(buf1);
        }
    }
    FFT2D(Rdat, Idat, N, LogN, 1);
}

void total_evol(double **Rdat, double **Idat, int N, int LogN, double dt, double dx,int num_step)
{
    int i;
    linear_evol(Rdat, Idat, N, LogN, dt / 2, dx);


    for (i = 0; i < num_step - 1; i++)
    {
        cubic_evol(Rdat, Idat, N, dt);
        linear_evol(Rdat, Idat, N, LogN, dt, dx);
    }

    cubic_evol(Rdat, Idat, N, dt);
    linear_evol(Rdat, Idat, N, LogN, dt / 2, dx);
}

double evol_check(double **Rdat, double **Idat, int N, int LogN, double dt, double dx,double ampl,int steps)
{
    double in_max, fin_max;

    in_max = sqrt(psi_sq(Rdat, Idat, N / 2, N / 2));

    total_evol(Rdat, Idat, N, LogN, dt, dx, steps);

    fin_max = sqrt(psi_sq(Rdat, Idat, N / 2, N / 2));

    return (in_max / fin_max);
}



double average_width(double **Rdat, double **Idat, int N, double dx)
{
    double s = 0.0;
    double p = 0.0;
    int i, j;
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < N; i++)
        {
            s += ((i - N / 2) * (i - N / 2) + (j - N / 2) *
                  (j - N / 2)) * dx * dx * psi_sq(Rdat, Idat, i, j) * dx;
            p += psi_sq(Rdat, Idat, i, j) * dx;
        }
    }

    s = s / p;
    return s;
}

void find_stat(double **Rdat, double **Idat, int N, int LogN, double dt, double dx, double ampl_min, double ampl_max, int step,int evol_step)
{
    double ampl;
    for (ampl = ampl_min; ampl < ampl_max; ampl += ((ampl_max - ampl_min) / step))
    {
        initiate(Rdat, Idat, N, dx, ampl);
        printf("%f %f %f\n", ampl, power(Rdat, Idat, N, dx),
                evol_check(Rdat, Idat, N, LogN, dt, dx, ampl, evol_step));
    }
}

int main()
{

    int i;
    int deg = 8;
    int cnum = 1;
    double **RePsi, **ImPsi;
    double size = 10.0;
    double ampl = 2.0;
    FILE *final, *start;
    double dt = 0.001;
    int n = 10;

    final = fopen("final.txt","w");
    start = fopen("start.txt","w");

    for (i = 0; i < deg; i++)
    {
        cnum *= 2;
    }

    double dx = size / cnum;

    RePsi = (double **) malloc (cnum * sizeof(double *));
    ImPsi = (double **) malloc (cnum * sizeof(double *));

    for (i = 0; i < cnum; i++)
    {
        RePsi[i] = (double *) malloc (cnum * sizeof(double));
        ImPsi[i] = (double *) malloc (cnum * sizeof(double));
    }

    //-------------------------------------------------------------------------
    initiate(RePsi,ImPsi,cnum,dx,ampl);

    for (i = 0; i < cnum; i++)
    {
        fprintf(start, "%d %f\n", i, sqrt(psi_sq(RePsi, ImPsi, cnum / 2, i)));
    }

    i = 0;

    printf("%f %f %f %f\n", i * n * dt, average_width(RePsi, ImPsi, cnum, dx),
            energy(RePsi, ImPsi, cnum, dx), power(RePsi, ImPsi, cnum, dx));

    for (i = 0;i<100;i++)
    {
        total_evol(RePsi, ImPsi, cnum, deg, dt, dx, n);
        printf("%f %.8f %f %f\n", (i + 1) * n * dt, average_width(RePsi, ImPsi, cnum, dx), energy(RePsi, ImPsi, cnum, dx), power(RePsi, ImPsi, cnum, dx));
    }


    for (i = 0; i < cnum; i++)
    {
        fprintf(final,"%d %f\n", i, sqrt(psi_sq(RePsi, ImPsi, cnum / 2, i)));
    }

    //-------------------------------------------------------------------------
    fclose(start);
    fclose(final);

    for (i = 0; i<cnum; i++)
    {
        free(RePsi[i]);
        free(ImPsi[i]);
    }
    free(RePsi);
    free(ImPsi);

    return 0;
}

