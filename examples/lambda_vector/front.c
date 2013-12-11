#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    const int n = 1000;

    FILE* fp_1;
    FILE* fp_2;

    fp_1 = fopen("front_1.txt", "w");
    fp_2 = fopen("front_2.txt", "w");

    int i;

    double delta = 10.0 / n;
    double lambda_1 = 0.0;
    double lambda_2 = 1.0;
    double alfa = 1.0;

    double m_1 = 1.0;
    double m_2 = 1.0;

    double bot_1[n], mid_1[n], top_1[n], diag_1[n];
    double bot_2[n], mid_2[n], top_2[n], diag_2[n];
    double bot_3[n], mid_3[n], top_3[n], diag_3[n];

    double norm_1, norm_p_1;
    double norm_2, norm_p_2;


    double col_1[n], psi_1[n], psi_p_1[n];
    double col_2[n], psi_2[n], psi_p_2[n];
    double col_3[n], theta[n];

    double S2 = 0.0;
    int fl = 0;

    while (lambda_2 >= 0.5)
    {
		S2 = 0.0;

		for (i = 0; i < n; i++)
        {
			psi_1[i] = pow((i * delta), abs(m_1)) * pow(2.71, -i * delta * (i * delta));
			psi_2[i] = pow((i * delta), abs(m_2)) *	pow(2.71, -i * delta * (i * delta));
			theta[i] = pow(2.71, -i * delta * i * delta);
        }

		for (i = 1; i < n; i++)
        {
			bot_1[i] = -1.0 / (delta * delta) + 1.0 / (2.0 * i * delta * delta);
			mid_1[i] = lambda_1 + 2 / (delta * delta) + m_1 * m_1 / (i * i * delta * delta);
			top_1[i] = -1.0 / (delta * delta) -1.0 / (2.0 * i * delta * delta);

			bot_2[i] = -1.0 / (delta * delta) + 1.0 / (2.0 * i * delta * delta);
			mid_2[i] = lambda_2 + 2 / (delta * delta) + m_2 * m_2 / (i * i * delta * delta);
			top_2[i] = -1.0 / (delta * delta) - 1.0 / (2.0 * i * delta * delta);

			bot_3[i] = -1.0 / (delta * delta) + 1.0 / (2.0 * i * delta * delta);
			top_3[i] = -1.0 / (delta * delta) - 1.0 / (2.0 * i * delta * delta);
			mid_3[i] = alfa + 2 / (delta * delta);
        }

		bot_1[0] = 0.0;
		top_1[0] = -2 / (delta * delta);

		bot_2[0] = 0.0;
		top_2[0] = -2 / (delta * delta);

		bot_3[0] = 0.0;
		top_3[0] = -2 / (delta * delta);
		mid_3[0] = alfa + 2 / (delta * delta);

		while (((S2 - 1) * (S2 - 1)) > 0.00001)
		{
			norm_p_1 = 0.0;
			norm_1   = 0.0;
			norm_p_2 = 0.0;
			norm_2   = 0.0;

			fl++;

			for (i = 0; i < n; i++)
			{
				diag_1[i]  = mid_1[i];
				diag_2[i]  = mid_2[i];

				col_1[i]   = psi_1[i] * theta[i];
				col_2[i]   = psi_2[i] * theta[i];

				psi_p_1[i] = psi_1[i];
				psi_p_2[i] = psi_2[i];
			}

			for (i = 2; i < n; i++)
			{
				diag_1[i] = diag_1[i] - top_1[i-1] * bot_1[i] / diag_1[i-1];
				diag_2[i] = diag_2[i] - top_2[i-1] * bot_2[i] / diag_2[i-1];

				col_1[i] = col_1[i] - col_1[i-1] * bot_1[i] / diag_1[i-1];
				col_2[i] = col_2[i] - col_2[i-1] * bot_2[i] / diag_2[i-1];
			}

			psi_1[n-1] = col_1[n-1] / diag_1[n-1];
			psi_2[n-1] = col_2[n-1] / diag_2[n-1];

			for (i = n - 2; i >= 1; i--)
			{
				col_1[i] = col_1[i] - col_1[i+1] * top_1[i] / diag_1[i+1];
				col_2[i] = col_2[i] - col_2[i+1] * top_2[i] / diag_2[i+1];

				psi_1[i] = col_1[i] / diag_1[i];
				psi_2[i] = col_2[i] / diag_2[i];
			}

			psi_1[0] = 0.0;

			for (i = 1; i < n; i++)
			{
				norm_p_1 += (psi_p_1[i] * psi_p_1[i] * i * delta * delta);
				norm_p_2 += (psi_p_2[i] * psi_p_2[i] * i * delta * delta);

				norm_1 += (psi_1[i] * psi_1[i] * i * delta * delta);
				norm_2 += (psi_2[i] * psi_2[i] * i * delta * delta);
			}

			S2 = pow(((norm_p_2 + norm_p_1) / (norm_2 + norm_1)), 0.75);

			for (i = 0; i < n; i++)
			{
				psi_1[i] *= S2;
				psi_2[i] *= S2;
			}

			for (i = 0; i < n; i++)
			{
				diag_3[i]  = mid_3[i];
				col_3[i]   = psi_1[i] * psi_1[i] + psi_2[i] * psi_2[i];
			}

			for (i = 1; i < n; i++)
			{
				diag_3[i] = diag_3[i] - top_3[i-1] * bot_3[i] / diag_3[i-1];
				col_3[i]  = col_3[i] - col_3[i-1] * bot_3[i] / diag_3[i-1];
			}

			theta[n-1] = col_3[n-1] / diag_3[n-1];

			for (i = n - 2; i >= 0; i--)
			{
				col_3[i] = col_3[i] - col_3[i+1] * top_3[i] / diag_3[i+1];
				theta[i] = col_3[i] / diag_3[i];
			}
		}

		printf("%lf %lf\n", lambda_1, norm_1);
		printf("%lf %lf\n", lambda_2, norm_2);

		fprintf(fp_1, "%lf %lf\n", lambda_1, norm_1);
		fprintf(fp_2, "%lf %lf\n", lambda_2, norm_2);
		lambda_2 -= 0.001;
		lambda_1 += 0.001;
    }

    fclose(fp_1);
    fclose(fp_2);

    return 0;
}
