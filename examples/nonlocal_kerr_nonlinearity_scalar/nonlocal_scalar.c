#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    const int n = 1000;

    FILE * fp_1;
    FILE * fp_t;
    fp_1 = fopen("sequence1.txt", "w");
    fp_t = fopen("sequencet.txt", "w");
    int i;

    double delta = 20.0 / n;
    double lambda = 0.5;
    double alfa = 1.0;

    double m = 1.0;

    double bot_1[n], mid_1[n], top_1[n], diag_1[n];
    double bot_3[n], mid_3[n], top_3[n], diag_3[n];

    double norm_1, norm_p_1;
    double norm_3, norm_p_3;

    double col_1[n], psi_1[n], psi_p_1[n];
    double col_3[n], theta[n], theta_p[n];

    double S1 = 0.0, S3 = 0.0;

    int fl = 0;

    for (i = 0; i < n; i++)
    {
        psi_1[i] = pow((i * delta), m) * pow(2.71, -i * delta * (i * delta));
        theta[i] = pow(2.71, -i * delta * i * delta);
    }

    for (i = 1; i < n; i++)
    {
        bot_1[i] =  -1.0 / (delta * delta) + 1.0 / (2.0 * i * delta * delta);
        mid_1[i] = lambda + 2 / (delta * delta) + m * m / (i * i * delta * delta);
        top_1[i] =  -1.0 / (delta * delta) - 1.0 / (2.0 * i * delta * delta);

        bot_3[i] =  -1.0 / (delta * delta) + 1.0 / (2.0 * i * delta * delta);
        top_3[i] =  -1.0 / (delta * delta) - 1.0 / (2.0 * i * delta * delta);
        mid_3[i] = alfa + 2 / (delta * delta);
    }

    bot_1[0] = 0.0;
    top_1[0] = -2 / (delta * delta);

    bot_3[0] = 0.0;
    top_3[0] = -2 / (delta * delta);
    mid_3[0] = alfa + 2 / (delta * delta);

    while (((S1 - 1) * (S1 - 1)) > 0.0000001)
    {
        norm_p_1 = 0.0;
        norm_1 = 0.0;
        norm_p_3 = 0.0;
        norm_3 = 0.0;
        fl++;

        for (i = 0; i < n; i++)
        {
            diag_1[i] = mid_1[i];

            col_1[i] = psi_1[i] * theta[i];

            psi_p_1[i] = psi_1[i];
        }

        for (i = 2; i < n; i++)
        {
            diag_1[i] = diag_1[i] - top_1[i-1] * bot_1[i] / diag_1[i-1];
            col_1[i] = col_1[i] - col_1[i-1] * bot_1[i] / diag_1[i-1];
        }

        psi_1[n-1] = col_1[n-1] / diag_1[n-1];

        for (i = n - 2; i >= 1; i--)
        {
            col_1[i] = col_1[i] - col_1[i+1] * top_1[i] / diag_1[i+1];
            psi_1[i] = col_1[i] / diag_1[i];
        }

        psi_1[0] = 0.0;

        for (i = 1; i < n; i++)
        {
            norm_p_1 += (psi_p_1[i] * psi_p_1[i] * i * delta * delta);
            norm_1 += (psi_1[i] * psi_1[i] * i * delta * delta);
        }

        S1 = pow((norm_p_1 / norm_1), 0.75);

        for (i = 0; i < n; i++)
        {
            psi_1[i] *= S1;
        }

        for (i = 0; i < n; i++)
        {
            diag_3[i] = mid_3[i];
            col_3[i] = psi_1[i] * psi_1[i];
            theta_p[i] = theta[i];
        }

        for (i = 1; i < n; i++)
        {
            diag_3[i] = diag_3[i] - top_3[i-1] * bot_3[i] / diag_3[i-1];
            col_3[i] = col_3[i] - col_3[i-1] * bot_3[i] / diag_3[i-1];
        }

        theta[n-1] = col_3[n-1] / diag_3[n-1];

        for (i = n-2; i >= 0; i--)
        {
            col_3[i] = col_3[i] - col_3[i+1] * top_3[i] / diag_3[i+1];
            theta[i] = col_3[i] / diag_3[i];
        }

        for (i = 1; i < n; i++)
        {
            norm_p_3 += (theta_p[i] * theta_p[i] * i * delta * delta);
            norm_3 += (theta[i] * theta[i] * i * delta * delta);
        }

        S3 = norm_p_3 / norm_3;

        /* 
        for (i = 0; i < n; i++)
        {
            theta[i] *  = S1;
        }
        */ 
        printf("%lf %lf\n", S1,  S3);
    }

    printf("%d\n", fl);

    for (i = 0; i < n; i++)
    {
        fprintf(fp_1, "%lf %lf\n", (i * delta), psi_1[i]);
        fprintf(fp_t, "%lf %lf\n", (i * delta), theta[i]);
    }

    fclose(fp_1);
    return 0;
}

