//    cLoops.c - C loops of fathon package
//    Copyright (C) 2019-2020  Stefano Bianchi
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
double flucDFAForwCompute(double *y, int curr_win_size, int N, int pol_ord)
{
    double *t = malloc(N * sizeof(double));
    for(int i = 0; i < N; i++)
    {
        t[i] = (double)(i + 1);
    }

    int N_s = N / curr_win_size;
    double f = 0.0;

    #pragma omp parallel for reduction(+ : f)
    for(int v = 0; v < N_s; v++)
    {
        int start_lim = v * curr_win_size;

        double t_fit[curr_win_size], y_fit[curr_win_size];
        for(int i = 0; i < curr_win_size; i++)
        {
            t_fit[i] = t[start_lim + i];
            y_fit[i] = y[start_lim + i];
        }

        double fit_coeffs[pol_ord + 1];
        polynomialFit(curr_win_size, pol_ord+1, t_fit, y_fit, fit_coeffs);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var = y_fit[j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var -= fit_coeffs[k] * pow(t_fit[j], k);
            }
            f += pow(var, 2.0);
        }
    }

    f = sqrt(f / (N_s * curr_win_size));

    free(t);

    return f;
}

//main loop for DFA (computes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the array y)
double flucDFAForwBackwCompute(double *y, int curr_win_size, int N, int pol_ord)
{
    double *t = malloc(N * sizeof(double));
    for(int i = 0; i < N; i++)
    {
        t[i] = (double)(i + 1);
    }

    int N_s = N / curr_win_size;
    double f = 0.0;

    #pragma omp parallel for reduction(+ : f)
    for(int v = 0; v < N_s; v++)
    {
        int start_lim = v * curr_win_size;

        double t_fit[curr_win_size], y_fit[curr_win_size];
        for(int i = 0; i < curr_win_size; i++)
        {
            t_fit[i] = t[start_lim + i];
            y_fit[i] = y[start_lim + i];
        }

        double fit_coeffs[pol_ord + 1];
        polynomialFit(curr_win_size, pol_ord+1, t_fit, y_fit, fit_coeffs);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y_fit[j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs[k] * pow(t_fit[j], k);
            }
            f += pow(var_1, 2.0);
        }

        start_lim = v * curr_win_size + (N - N_s * curr_win_size);

        for(int i = 0; i < curr_win_size; i++)
        {
            t_fit[i] = t[start_lim + i];
            y_fit[i] = y[start_lim + i];
        }

        polynomialFit(curr_win_size, pol_ord+1, t_fit, y_fit, fit_coeffs);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_2 = y_fit[j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_2 -= fit_coeffs[k] * pow(t_fit[j], k);
            }
            f += pow(var_2, 2.0);
        }
    }

    f = sqrt(f / (2.0 * N_s * curr_win_size));

    free(t);

    return f;
}

//main loop for MFDFA (computes fluctuations starting from the beginning of the array y)
double flucMFDFAForwCompute(double *y, int curr_win_size, double q, int N, int pol_ord)
{
    double *t = malloc(N * sizeof(double));
    for(int i = 0; i < N; i++)
    {
        t[i] = (double)(i + 1);
    }

    int N_s = N / curr_win_size;
    double f = 0.0;

    #pragma omp parallel for reduction(+ : f)
    for(int v = 0; v < N_s; v++)
    {
        double rms = 0.0;
        int start_lim = v * curr_win_size;

        double t_fit[curr_win_size], y_fit[curr_win_size];
        for(int i = 0; i < curr_win_size; i++)
        {
            t_fit[i] = t[start_lim + i];
            y_fit[i] = y[start_lim + i];
        }

        double fit_coeffs[pol_ord + 1];
        polynomialFit(curr_win_size, pol_ord+1, t_fit, y_fit, fit_coeffs);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var = y_fit[j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var -= fit_coeffs[k] * pow(t_fit[j], k);
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
    }

    if(q == 0.0)
    {
        f = exp(f / (double)(2 * N_s));
    }
    else
    {
        f = pow(f / (double)N_s, 1 / (double)q);
    }

    free(t);

    return f;
}

//main loop for MFDFA (computes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the array y)
double flucMFDFAForwBackwCompute(double *y, int curr_win_size, double q, int N, int pol_ord)
{
    double *t = malloc(N * sizeof(double));
    for(int i = 0; i < N; i++)
    {
        t[i] = (double)(i + 1);
    }

    int N_s = N / curr_win_size;
    double f = 0.0;

    #pragma omp parallel for reduction(+ : f)
    for(int v = 0; v < N_s; v++)
    {
        double rms1 = 0.0;
        double rms2 = 0.0;
        int start_lim = v * curr_win_size;

        double t_fit[curr_win_size], y_fit[curr_win_size];
        for(int i = 0; i < curr_win_size; i++)
        {
            t_fit[i] = t[start_lim + i];
            y_fit[i] = y[start_lim + i];
        }

        double fit_coeffs[pol_ord + 1];
        polynomialFit(curr_win_size, pol_ord+1, t_fit, y_fit, fit_coeffs);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_1 = y_fit[j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs[k] * pow(t_fit[j], k);
            }
            rms1 += pow(var_1, 2.0);
        }

        start_lim = v * curr_win_size + (N - N_s * curr_win_size);

        for(int i = 0; i < curr_win_size; i++)
        {
            t_fit[i] = t[start_lim + i];
            y_fit[i] = y[start_lim + i];
        }

        polynomialFit(curr_win_size, pol_ord+1, t_fit, y_fit, fit_coeffs);

        for(int j = 0; j < curr_win_size; j++)
        {
            double var_2 = y_fit[j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_2 -= fit_coeffs[k] * pow(t_fit[j], k);
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
    }

    if(q == 0.0)
    {
        f = exp(f / (double)(4 * N_s));
    }
    else
    {
        f = pow(f / (double)(2 * N_s), 1 / (double)q);
    }

    free(t);

    return f;
}

//main loop for DCCA (computes fluctuations using absolute values)
double flucDCCAAbsCompute(double *y1, double *y2, int curr_win_size, int N, int pol_ord)
{
    double *t = malloc(N * sizeof(double));
    for(int i = 0; i < N; i++)
    {
        t[i] = (double)(i + 1);
    }

    int N_s = N - curr_win_size;
    double f = 0.0;

    #pragma omp parallel for reduction(+ : f)
    for(int v = 0; v < N_s; v++)
    {
        double t_fit[curr_win_size + 1], y_fit1[curr_win_size + 1], y_fit2[curr_win_size + 1];
        for(int i = 0; i <= curr_win_size; i++)
        {
            t_fit[i] = t[v + i];
            y_fit1[i] = y1[v + i];
            y_fit2[i] = y2[v + i];
        }

        double fit_coeffs1[pol_ord + 1], fit_coeffs2[pol_ord + 1];
        polynomialFit(curr_win_size+1, pol_ord+1, t_fit, y_fit1, fit_coeffs1);
        polynomialFit(curr_win_size+1, pol_ord+1, t_fit, y_fit2, fit_coeffs2);

        for(int j = 0; j <= curr_win_size; j++)
        {
            double var_1 = y_fit1[j];
            double var_2 = y_fit2[j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs1[k] * pow(t_fit[j], k);
                var_2 -= fit_coeffs2[k] * pow(t_fit[j], k);
            }
            f += fabs(var_1 * var_2);
        }
    }

    f = sqrt(f / (N_s * (curr_win_size - 1)));

    free(t);

    return f;
}

//main loop for DCCA (computes fluctuations without using absolute values)
double flucDCCANoAbsCompute(double *y1, double *y2, int curr_win_size, int N, int pol_ord)
{
    double *t = malloc(N * sizeof(double));
    for(int i = 0; i < N; i++)
    {
        t[i] = (double)(i + 1);
    }

    int N_s = N - curr_win_size;
    double f = 0.0;

    #pragma omp parallel for reduction(+ : f)
    for(int v = 0; v < N_s; v++)
    {
        double t_fit[curr_win_size + 1], y_fit1[curr_win_size + 1], y_fit2[curr_win_size + 1];
        for(int i = 0; i <= curr_win_size; i++)
        {
            t_fit[i] = t[v + i];
            y_fit1[i] = y1[v + i];
            y_fit2[i] = y2[v + i];
        }

        double fit_coeffs1[pol_ord + 1], fit_coeffs2[pol_ord + 1];
        polynomialFit(curr_win_size+1, pol_ord+1, t_fit, y_fit1, fit_coeffs1);
        polynomialFit(curr_win_size+1, pol_ord+1, t_fit, y_fit2, fit_coeffs2);

        for(int j = 0; j <= curr_win_size; j++)
        {
            double var_1 = y_fit1[j];
            double var_2 = y_fit2[j];
            for(int k = 0; k < (pol_ord + 1); k++)
            {
                var_1 -= fit_coeffs1[k] * pow(t_fit[j], k);
                var_2 -= fit_coeffs2[k] * pow(t_fit[j], k);
            }
            f += var_1 * var_2;
        }
    }

    f = f / (N_s * (curr_win_size - 1));

    free(t);

    return f;
}

//main loop for HT (computes fluctuations)
double HTCompute(double *y, int scale, int N, int pol_ord, int v)
{
    double *t = malloc(N * sizeof(double));
    for(int i = 0; i < N; i++)
    {
        t[i] = (double)(i + 1);
    }

    double f = 0.0;
    double t_fit[scale], y_fit[scale];
    for(int i = 0; i < scale; i++)
    {
        t_fit[i] = t[v + i];
        y_fit[i] = y[v + i];
    }

    double fit_coeffs[pol_ord + 1];
    polynomialFit(scale, pol_ord+1, t_fit, y_fit, fit_coeffs);

    for(int j = 0; j < scale; j++)
    {
        double var = y_fit[j];
        for(int k = 0; k < (pol_ord + 1); k++)
        {
            var -= fit_coeffs[k] * pow(t_fit[j], k);
        }
        f += pow(var, 2.0);
    }

    f = sqrt(f / (double)scale);

    free(t);

    return f;
}
