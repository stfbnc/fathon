//    cLoops.c - C loops of fathon package
//    Copyright (C) 2019-2021  Stefano Bianchi
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "cFuncs.h"
#include "cLoops.h"
#include "omp.h"

//main loop for DFA (computes fluctuations starting from the beginning of the array y)
double flucDFAForwCompute(double *y, double *t, int curr_win_size, int N, int pol_ord)
{
    int N_s = N / curr_win_size;
    double f = 0.0;
#ifdef _WIN64
    int v = 0;
#endif

    #pragma omp parallel for reduction(+ : f)
#ifdef _WIN64
    for(v = 0; v < N_s; v++)
#else
    for(int v = 0; v < N_s; v++)
#endif
    {
        int start_lim = v * curr_win_size;
        double *fit_coeffs = malloc((pol_ord + 1) * sizeof(double));
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y + start_lim, fit_coeffs);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var = y[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var -= fit_coeffs[k] * pow(t[start_lim + j], k);
            }
            f += pow(var, 2.0);
        }

        free(fit_coeffs);
    }

    f = sqrt(f / (N_s * curr_win_size));

    return f;
}

//main loop for DFA (computes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the array y)
double flucDFAForwBackwCompute(double *y, double *t, int curr_win_size, int N, int pol_ord)
{
    int N_s = N / curr_win_size;
    double f = 0.0;
#ifdef _WIN64
    int v = 0;
#endif

    #pragma omp parallel for reduction(+ : f)
#ifdef _WIN64
    for(v = 0; v < N_s; v++)
#else
    for(int v = 0; v < N_s; v++)
#endif
    {
        int start_lim = v * curr_win_size;
        double *fit_coeffs = malloc((pol_ord + 1) * sizeof(double));
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y + start_lim, fit_coeffs);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs[k] * pow(t[start_lim + j], k);
            }
            f += pow(var_1, 2.0);
        }

        start_lim = v * curr_win_size + (N - N_s * curr_win_size);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y + start_lim, fit_coeffs);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_2 = y[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_2 -= fit_coeffs[k] * pow(t[start_lim + j], k);
            }
            f += pow(var_2, 2.0);
        }

        free(fit_coeffs);
    }

    f = sqrt(f / (2.0 * N_s * curr_win_size));

    return f;
}

//main loop for MFDFA (computes fluctuations starting from the beginning of the array y)
double flucMFDFAForwCompute(double *y, double *t, int curr_win_size, double q, int N, int pol_ord)
{
    int N_s = N / curr_win_size;
    double f = 0.0;
#ifdef _WIN64
    int v = 0;
#endif

    #pragma omp parallel for reduction(+ : f)
#ifdef _WIN64
    for(v = 0; v < N_s; v++)
#else
    for(int v = 0; v < N_s; v++)
#endif
    {
        double rms = 0.0;
        int start_lim = v * curr_win_size;
        double *fit_coeffs = malloc((pol_ord + 1) * sizeof(double));
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y + start_lim, fit_coeffs);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var = y[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var -= fit_coeffs[k] * pow(t[start_lim + j], k);
            }
            rms += pow(var, 2.0);
        }

        if(q == 0.0)
        {
            f += log(rms / (double)curr_win_size);
        }
        else
        {
            f += pow(rms / (double)curr_win_size, 0.5 * q);
        }

        free(fit_coeffs);
    }

    if(q == 0.0)
    {
        f = exp(f / (double)(2 * N_s));
    }
    else
    {
        f = pow(f / (double)N_s, 1 / (double)q);
    }

    return f;
}

//main loop for MFDFA (computes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the array y)
double flucMFDFAForwBackwCompute(double *y, double *t, int curr_win_size, double q, int N, int pol_ord)
{
    int N_s = N / curr_win_size;
    double f = 0.0;
#ifdef _WIN64
    int v = 0;
#endif

    #pragma omp parallel for reduction(+ : f)
#ifdef _WIN64
    for(v = 0; v < N_s; v++)
#else
    for(int v = 0; v < N_s; v++)
#endif
    {
        double rms1 = 0.0;
        double rms2 = 0.0;
        int start_lim = v * curr_win_size;
        double *fit_coeffs = malloc((pol_ord + 1) * sizeof(double));
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y + start_lim, fit_coeffs);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs[k] * pow(t[start_lim + j], k);
            }
            rms1 += pow(var_1, 2.0);
        }

        start_lim = v * curr_win_size + (N - N_s * curr_win_size);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y + start_lim, fit_coeffs);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_2 = y[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_2 -= fit_coeffs[k] * pow(t[start_lim + j], k);
            }
            rms2 += pow(var_2, 2.0);
        }

        if(q == 0.0)
        {
            f += (log(rms1 / (double)curr_win_size) + log(rms2 / (double)curr_win_size));
        }
        else
        {
            f += (pow(rms1 / (double)curr_win_size, 0.5 * q) + pow(rms2 / (double)curr_win_size, 0.5 * q));
        }

        free(fit_coeffs);
    }

    if(q == 0.0)
    {
        f = exp(f / (double)(4 * N_s));
    }
    else
    {
        f = pow(f / (double)(2 * N_s), 1 / (double)q);
    }

    return f;
}

//main loop for DCCA (computes fluctuations using absolute values)
double flucDCCAAbsCompute(double *y1, double *y2, double *t, int curr_win_size, int N, int pol_ord)
{
    int N_s = N - curr_win_size;
    double f = 0.0;
#ifdef _WIN64
    int v = 0;
#endif

    #pragma omp parallel for reduction(+ : f)
#ifdef _WIN64
    for(v = 0; v < N_s; v++)
#else
    for(int v = 0; v < N_s; v++)
#endif
    {
        double *fit_coeffs1 = malloc((pol_ord + 1) * sizeof(double));
        double *fit_coeffs2 = malloc((pol_ord + 1) * sizeof(double));
        polynomialFit(curr_win_size + 1, pol_ord + 1, t + v, y1 + v, fit_coeffs1);
        polynomialFit(curr_win_size + 1, pol_ord + 1, t + v, y2 + v, fit_coeffs2);

        for(int j = 0; j <= curr_win_size; j++)
        {
            double var_1 = y1[v + j];
            double var_2 = y2[v + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs1[k] * pow(t[v + j], k);
                var_2 -= fit_coeffs2[k] * pow(t[v + j], k);
            }
            f += fabs(var_1 * var_2);
        }

        free(fit_coeffs1);
        free(fit_coeffs2);
    }

    f = sqrt(f / (N_s * (curr_win_size - 1)));

    return f;
}

//main loop for DCCA (computes fluctuations without using absolute values)
double flucDCCANoAbsCompute(double *y1, double *y2, double *t, int curr_win_size, int N, int pol_ord)
{
    int N_s = N - curr_win_size;
    double f = 0.0;
#ifdef _WIN64
    int v = 0;
#endif

    #pragma omp parallel for reduction(+ : f)
#ifdef _WIN64
    for(v = 0; v < N_s; v++)
#else
    for(int v = 0; v < N_s; v++)
#endif
    {
        double *fit_coeffs1 = malloc((pol_ord + 1) * sizeof(double));
        double *fit_coeffs2 = malloc((pol_ord + 1) * sizeof(double));
        polynomialFit(curr_win_size + 1, pol_ord + 1, t + v, y1 + v, fit_coeffs1);
        polynomialFit(curr_win_size + 1, pol_ord + 1, t + v, y2 + v, fit_coeffs2);

        for(int j = 0; j <= curr_win_size; j++)
        {
            double var_1 = y1[v + j];
            double var_2 = y2[v + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs1[k] * pow(t[v + j], k);
                var_2 -= fit_coeffs2[k] * pow(t[v + j], k);
            }
            f += var_1 * var_2;
        }

        free(fit_coeffs1);
        free(fit_coeffs2);
    }

    f = f / (N_s * (curr_win_size - 1));

    return f;
}

//main loop for HT (computes fluctuations)
double HTCompute(double *y, double *t, int scale, int N, int pol_ord, int v)
{
    double f = 0.0;
    double *fit_coeffs = malloc((pol_ord + 1) * sizeof(double));
    polynomialFit(scale, pol_ord + 1, t + v, y + v, fit_coeffs);

    for(int j = 0; j < scale; j++)
    {
        double var = y[v + j];
        for(int k = 0; k < (pol_ord + 1); k++)
        {
            var -= fit_coeffs[k] * pow(t[v + j], k);
        }
        f += pow(var, 2.0);
    }

    f = sqrt(f / (double)scale);

    free(fit_coeffs);

    return f;
}

//main loop for DCCA without overlap (computes fluctuations starting from the beginning
// of the array y and using absolute values)
double flucDCCAForwAbsComputeNoOverlap(double *y1, double *y2, double *t, int curr_win_size, int N, int pol_ord)
{
    int N_s = N / curr_win_size;
    double f = 0.0;
#ifdef _WIN64
    int v = 0;
#endif

    #pragma omp parallel for reduction(+ : f)
#ifdef _WIN64
    for(v = 0; v < N_s; v++)
#else
    for(int v = 0; v < N_s; v++)
#endif
    {
        int start_lim = v * curr_win_size;
        double *fit_coeffs_1 = malloc((pol_ord + 1) * sizeof(double));
        double *fit_coeffs_2 = malloc((pol_ord + 1) * sizeof(double));
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y1 + start_lim, fit_coeffs_1);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y2 + start_lim, fit_coeffs_2);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y1[start_lim + j];
            double var_2 = y2[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs_1[k] * pow(t[start_lim + j], k);
                var_2 -= fit_coeffs_2[k] * pow(t[start_lim + j], k);
            }
            f += fabs(var_1 * var_2);
        }

        free(fit_coeffs_1);
        free(fit_coeffs_2);
    }

    f = sqrt(f / (N_s * curr_win_size));

    return f;
}

//main loop for DCCA without overlap (ccomputes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the arrays y1 and y2, and using absolute values)
double flucDCCAForwBackwAbsComputeNoOverlap(double *y1, double *y2, double *t, int curr_win_size, int N, int pol_ord)
{
    int N_s = N / curr_win_size;
    double f = 0.0;
#ifdef _WIN64
    int v = 0;
#endif

    #pragma omp parallel for reduction(+ : f)
#ifdef _WIN64
    for(v = 0; v < N_s; v++)
#else
    for(int v = 0; v < N_s; v++)
#endif
    {
        int start_lim = v * curr_win_size;
        double *fit_coeffs_1 = malloc((pol_ord + 1) * sizeof(double));
        double *fit_coeffs_2 = malloc((pol_ord + 1) * sizeof(double));
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y1 + start_lim, fit_coeffs_1);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y2 + start_lim, fit_coeffs_2);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y1[start_lim + j];
            double var_2 = y2[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs_1[k] * pow(t[start_lim + j], k);
                var_2 -= fit_coeffs_2[k] * pow(t[start_lim + j], k);
            }
            f += fabs(var_1 * var_2);
        }

        start_lim = v * curr_win_size + (N - N_s * curr_win_size);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y1 + start_lim, fit_coeffs_1);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y2 + start_lim, fit_coeffs_2);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y1[start_lim + j];
            double var_2 = y2[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs_1[k] * pow(t[start_lim + j], k);
                var_2 -= fit_coeffs_2[k] * pow(t[start_lim + j], k);
            }
            f += fabs(var_1 * var_2);
        }

        free(fit_coeffs_1);
        free(fit_coeffs_2);
    }

    f = sqrt(f / (2.0 * N_s * curr_win_size));

    return f;
}

//main loop for DCCA without overlap (computes fluctuations starting from the beginning
// of the array y)
double flucDCCAForwNoAbsComputeNoOverlap(double *y1, double *y2, double *t, int curr_win_size, int N, int pol_ord)
{
    int N_s = N / curr_win_size;
    double f = 0.0;
#ifdef _WIN64
    int v = 0;
#endif

    #pragma omp parallel for reduction(+ : f)
#ifdef _WIN64
    for(v = 0; v < N_s; v++)
#else
    for(int v = 0; v < N_s; v++)
#endif
    {
        int start_lim = v * curr_win_size;
        double *fit_coeffs_1 = malloc((pol_ord + 1) * sizeof(double));
        double *fit_coeffs_2 = malloc((pol_ord + 1) * sizeof(double));
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y1 + start_lim, fit_coeffs_1);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y2 + start_lim, fit_coeffs_2);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y1[start_lim + j];
            double var_2 = y2[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs_1[k] * pow(t[start_lim + j], k);
                var_2 -= fit_coeffs_2[k] * pow(t[start_lim + j], k);
            }
            f += var_1 * var_2;
        }

        free(fit_coeffs_1);
        free(fit_coeffs_2);
    }

    f = f / (N_s * curr_win_size);

    return f;
}

//main loop for DCCA without overlap (computes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the arrays y1 and y2)
double flucDCCAForwBackwNoAbsComputeNoOverlap(double *y1, double *y2, double *t, int curr_win_size, int N, int pol_ord)
{
    int N_s = N / curr_win_size;
    double f = 0.0;
#ifdef _WIN64
    int v = 0;
#endif

    #pragma omp parallel for reduction(+ : f)
#ifdef _WIN64
    for(v = 0; v < N_s; v++)
#else
    for(int v = 0; v < N_s; v++)
#endif
    {
        int start_lim = v * curr_win_size;
        double *fit_coeffs_1 = malloc((pol_ord + 1) * sizeof(double));
        double *fit_coeffs_2 = malloc((pol_ord + 1) * sizeof(double));
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y1 + start_lim, fit_coeffs_1);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y2 + start_lim, fit_coeffs_2);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y1[start_lim + j];
            double var_2 = y2[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs_1[k] * pow(t[start_lim + j], k);
                var_2 -= fit_coeffs_2[k] * pow(t[start_lim + j], k);
            }
            f += var_1 * var_2;
        }

        start_lim = v * curr_win_size + (N - N_s * curr_win_size);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y1 + start_lim, fit_coeffs_1);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y2 + start_lim, fit_coeffs_2);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y1[start_lim + j];
            double var_2 = y2[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs_1[k] * pow(t[start_lim + j], k);
                var_2 -= fit_coeffs_2[k] * pow(t[start_lim + j], k);
            }
            f += var_1 * var_2;
        }

        free(fit_coeffs_1);
        free(fit_coeffs_2);
    }

    f = f / (2.0 * N_s * curr_win_size);

    return f;
}

//main loop for MFDCCA (computes fluctuations starting from the beginning of the array y)
double flucMFDCCAForwCompute(double *y1, double *y2, double *t, int curr_win_size, double q, int N, int pol_ord)
{
    int N_s = N / curr_win_size;
    double f = 0.0;
#ifdef _WIN64
    int v = 0;
#endif

    #pragma omp parallel for reduction(+ : f)
#ifdef _WIN64
    for(v = 0; v < N_s; v++)
#else
    for(int v = 0; v < N_s; v++)
#endif
    {
        double rms = 0.0;
        int start_lim = v * curr_win_size;
        double *fit_coeffs_1 = malloc((pol_ord + 1) * sizeof(double));
        double *fit_coeffs_2 = malloc((pol_ord + 1) * sizeof(double));
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y1 + start_lim, fit_coeffs_1);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y2 + start_lim, fit_coeffs_2);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y1[start_lim + j];
            double var_2 = y2[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs_1[k] * pow(t[start_lim + j], k);
                var_2 -= fit_coeffs_2[k] * pow(t[start_lim + j], k);
            }
            rms += fabs(var_1 * var_2);
        }

        if(q == 0.0)
        {
            f += log(rms / (double)curr_win_size);
        }
        else
        {
            f += pow(rms / (double)curr_win_size, 0.5 * q);
        }

        free(fit_coeffs_1);
        free(fit_coeffs_2);
    }

    if(q == 0.0)
    {
        f = exp(f / (double)(2 * N_s));
    }
    else
    {
        f = pow(f / (double)N_s, 1 / (double)q);
    }

    return f;
}

//main loop for MFDCCA (computes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the array y)
double flucMFDCCAForwBackwCompute(double *y1, double *y2, double *t, int curr_win_size, double q, int N, int pol_ord)
{
    int N_s = N / curr_win_size;
    double f = 0.0;
#ifdef _WIN64
    int v = 0;
#endif

    #pragma omp parallel for reduction(+ : f)
#ifdef _WIN64
    for(v = 0; v < N_s; v++)
#else
    for(int v = 0; v < N_s; v++)
#endif
    {
        double rms1 = 0.0;
        double rms2 = 0.0;
        int start_lim = v * curr_win_size;
        double *fit_coeffs_1 = malloc((pol_ord + 1) * sizeof(double));
        double *fit_coeffs_2 = malloc((pol_ord + 1) * sizeof(double));
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y1 + start_lim, fit_coeffs_1);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y2 + start_lim, fit_coeffs_2);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y1[start_lim + j];
            double var_2 = y2[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs_1[k] * pow(t[start_lim + j], k);
                var_2 -= fit_coeffs_2[k] * pow(t[start_lim + j], k);
            }
            rms1 += fabs(var_1 * var_2);
        }

        start_lim = v * curr_win_size + (N - N_s * curr_win_size);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y1 + start_lim, fit_coeffs_1);
        polynomialFit(curr_win_size, pol_ord + 1, t + start_lim, y2 + start_lim, fit_coeffs_2);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y1[start_lim + j];
            double var_2 = y2[start_lim + j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs_1[k] * pow(t[start_lim + j], k);
                var_2 -= fit_coeffs_2[k] * pow(t[start_lim + j], k);
            }
            rms2 += fabs(var_1 * var_2);
        }

        if(q == 0.0)
        {
            f += (log(rms1 / (double)curr_win_size) + log(rms2 / (double)curr_win_size));
        }
        else
        {
            f += (pow(rms1 / (double)curr_win_size, 0.5 * q) + pow(rms2 / (double)curr_win_size, 0.5 * q));
        }

        free(fit_coeffs_1);
        free(fit_coeffs_2);
    }

    if(q == 0.0)
    {
        f = exp(f / (double)(4 * N_s));
    }
    else
    {
        f = pow(f / (double)(2 * N_s), 1 / (double)q);
    }

    return f;
}
