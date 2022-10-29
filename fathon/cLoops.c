//    cLoops.c - C loops of fathon package
//    Copyright (C) 2019-  Stefano Bianchi
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

//main loop for unbiased DFA
void flucUDFACompute(double *y_vec, double *t_vec, int y_len, int *wins_vec, int num_wins, int pol, double *f_vec)
{
#ifdef _WIN64
    int i = 0;
#endif

#pragma omp parallel for
#ifdef _WIN64
    for(i = 0; i < num_wins; i++)
#else
    for(int i = 0; i < num_wins; i++)
#endif
    {
        int s = wins_vec[i];
        int n_wins = y_len - s + 1;

        double *fit_coeffs = malloc((pol + 1) * sizeof(double));
        double *df = malloc(s * sizeof(double));

        double f = 0.0;
#ifdef _WIN64
        int start = 0;
        for(start = 0; start < n_wins; start++)
#else
        for(int start = 0; start < n_wins; start++)
#endif
        {
            polynomialFit(s, pol + 1, t_vec + start, y_vec + start, fit_coeffs);
            for(int j = 0; j < s; j++)
            {
                df[j] = y_vec[start + j];
                for(int k = 0; k < (pol + 1); k++)
                {
                    df[j] -= fit_coeffs[k] * pow(t_vec[start + j], k);
                }
            }
        
            double df_sum = 0.0, df_2_sum = 0.0, df_even_sum = 0.0, df_odd_sum = 0.0, df_shift_sum = 0.0;
            for(int j = 0; j < s; j++)
            {
                df_sum += df[j];
                df_2_sum += df[j] * df[j];
            }
            for(int j = 0; j < s; j += 2)
            {
                df_odd_sum += df[j];
            }
            for(int j = 1; j < s; j += 2)
            {
                df_even_sum += df[j];
            }
            for(int j = 0; j < (s - 1); j++)
            {
                df_shift_sum += (df[j] * df[j + 1]);
            }
    
            double df_neg_mean = (df_odd_sum - df_even_sum) / (double)s;
            double df_neg_var = df_2_sum / (double)s - df_neg_mean * df_neg_mean;
            double df_pos_mean = df_sum / (double)s;
            double df_pos_var = df_2_sum / (double)s - df_pos_mean * df_pos_mean;
        
            double df_pos_shift = (df_shift_sum + df_pos_mean * (df[0] + df[s - 1] - df_pos_mean * (s + 1))) / df_pos_var;
            double df_neg_shift = (-df_shift_sum + df_neg_mean * (df[0] + pow(-1.0, s + 1) * df[s - 1] - df_neg_mean * (s + 1))) / df_neg_var;
            double rho_A = (s + df_pos_shift) / (double)(2 * s - 1);
            double rho_B = (s + df_neg_shift) / (double)(2 * s - 1);
        
            double rho_A_star = rho_A + (1 + 3 * rho_A) / (double)(2 * s);
            double rho_B_star = rho_B + (1 + 3 * rho_B) / (double)(2 * s);
            f += ((rho_A_star + rho_B_star) * (1 - 1.0 / (double)(2 * s)) * df_pos_var);
        }
        f_vec[i] = sqrt(f * sqrt((s - 1) / (double)s) / (double)(n_wins));

        free(fit_coeffs);
        free(df);
    }
}

//main loop for DFA (computes fluctuations starting from the beginning of the array y)
void flucDFAForwCompute(double *y, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec)
{
#ifdef _WIN64
    int i = 0;
#endif

#pragma omp parallel for
#ifdef _WIN64
    for(i = 0; i < n_wins; i++)
#else
    for(int i = 0; i < n_wins; i++)
#endif
    {
        int curr_win_size = wins[i];
        int N_s = N / curr_win_size;
        double f = 0.0;
#ifdef _WIN64
        int v = 0;
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
    
        f_vec[i] = sqrt(f / (N_s * curr_win_size));
    }
}

//main loop for DFA (computes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the array y)
void flucDFAForwBackwCompute(double *y, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec)
{
#ifdef _WIN64
    int i = 0;
#endif

#pragma omp parallel for
#ifdef _WIN64
    for(i = 0; i < n_wins; i++)
#else
    for(int i = 0; i < n_wins; i++)
#endif
    {
        int curr_win_size = wins[i];
        int N_s = N / curr_win_size;
        double f = 0.0;
#ifdef _WIN64
        int v = 0;
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
    
        f_vec[i] = sqrt(f / (2.0 * N_s * curr_win_size));
    }
}

//main loop for MFDFA (computes fluctuations starting from the beginning of the array y)
void flucMFDFAForwCompute(double *y, double *t, int N, int *wins, int n_wins, double *qs, int n_q, int pol_ord, double *f_vec)
{
#ifdef _WIN64
    int iq = 0;
#endif

#ifdef _WIN64
#pragma omp parallel for
    for(iq = 0; iq < n_q; iq++)
    {
        int i = 0;
        for(i = 0; i < n_wins; i++)
#else
#pragma omp parallel for collapse(2)
    for(int iq = 0; iq < n_q; iq++)
    {
        for(int i = 0; i < n_wins; i++)
#endif
        {
            double q = qs[iq];
            int curr_win_size = wins[i];
            int N_s = N / curr_win_size;
            double f = 0.0;
#ifdef _WIN64
            int v = 0;
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
                f_vec[iq * n_wins + i] = exp(f / (double)(2 * N_s));
            }
            else
            {
                f_vec[iq * n_wins + i] = pow(f / (double)N_s, 1 / (double)q);
            }
        }
    }
}

//main loop for MFDFA (computes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the array y)
void flucMFDFAForwBackwCompute(double *y, double *t, int N, int *wins, int n_wins, double *qs, int n_q, int pol_ord, double *f_vec)
{
#ifdef _WIN64
    int iq = 0;
#endif

#ifdef _WIN64
#pragma omp parallel for
    for(iq = 0; iq < n_q; iq++)
    {
        int i = 0;
        for(i = 0; i < n_wins; i++)
#else
#pragma omp parallel for collapse(2)
    for(int iq = 0; iq < n_q; iq++)
    {
        for(int i = 0; i < n_wins; i++)
#endif
        {
            double q = qs[iq];
            int curr_win_size = wins[i];
            int N_s = N / curr_win_size;
            double f = 0.0;
#ifdef _WIN64
            int v = 0;
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
                f_vec[iq * n_wins + i] = exp(f / (double)(4 * N_s));
            }
            else
            {
                f_vec[iq * n_wins + i] = pow(f / (double)(2 * N_s), 1 / (double)q);
            }
        }
    }
}

//main loop for DCCA (computes fluctuations using absolute values)
void flucDCCAAbsCompute(double *y1, double *y2, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec)
{
#ifdef _WIN64
    int i = 0;
#endif

#pragma omp parallel for
#ifdef _WIN64
    for(i = 0; i < n_wins; i++)
#else
    for(int i = 0; i < n_wins; i++)
#endif
    {
        int curr_win_size = wins[i];
        int N_s = N - curr_win_size;
        double f = 0.0;
#ifdef _WIN64
        int v = 0;
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

        f_vec[i] = sqrt(f / (N_s * (curr_win_size - 1)));
    }
}

//main loop for DCCA (computes fluctuations without using absolute values)
void flucDCCANoAbsCompute(double *y1, double *y2, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec)
{
#ifdef _WIN64
    int i = 0;
#endif

#pragma omp parallel for
#ifdef _WIN64
    for(i = 0; i < n_wins; i++)
#else
    for(int i = 0; i < n_wins; i++)
#endif
    {
        int curr_win_size = wins[i];
        int N_s = N - curr_win_size;
        double f = 0.0;
#ifdef _WIN64
        int v = 0;
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

        f_vec[i] = f / (N_s * (curr_win_size - 1));
    }
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
void flucDCCAForwAbsComputeNoOverlap(double *y1, double *y2, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec)
{
#ifdef _WIN64
    int i = 0;
#endif

#pragma omp parallel for
#ifdef _WIN64
    for(i = 0; i < n_wins; i++)
#else
    for(int i = 0; i < n_wins; i++)
#endif
    {
        int curr_win_size = wins[i];
        int N_s = N / curr_win_size;
        double f = 0.0;
#ifdef _WIN64
        int v = 0;
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

        f_vec[i] = sqrt(f / (N_s * curr_win_size));
    }
}

//main loop for DCCA without overlap (ccomputes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the arrays y1 and y2, and using absolute values)
void flucDCCAForwBackwAbsComputeNoOverlap(double *y1, double *y2, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec)
{
#ifdef _WIN64
    int i = 0;
#endif

#pragma omp parallel for
#ifdef _WIN64
    for(i = 0; i < n_wins; i++)
#else
    for(int i = 0; i < n_wins; i++)
#endif
    {
        int curr_win_size = wins[i];
        int N_s = N / curr_win_size;
        double f = 0.0;
#ifdef _WIN64
        int v = 0;
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

        f_vec[i] = sqrt(f / (2.0 * N_s * curr_win_size));
    }
}

//main loop for DCCA without overlap (computes fluctuations starting from the beginning
// of the array y)
void flucDCCAForwNoAbsComputeNoOverlap(double *y1, double *y2, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec)
{
#ifdef _WIN64
    int i = 0;
#endif

#pragma omp parallel for
#ifdef _WIN64
    for(i = 0; i < n_wins; i++)
#else
    for(int i = 0; i < n_wins; i++)
#endif
    {
        int curr_win_size = wins[i];
        int N_s = N / curr_win_size;
        double f = 0.0;
#ifdef _WIN64
        int v = 0;
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

        f_vec[i] = f / (N_s * curr_win_size);
    }
}

//main loop for DCCA without overlap (computes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the arrays y1 and y2)
void flucDCCAForwBackwNoAbsComputeNoOverlap(double *y1, double *y2, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec)
{
#ifdef _WIN64
    int i = 0;
#endif

#pragma omp parallel for
#ifdef _WIN64
    for(i = 0; i < n_wins; i++)
#else
    for(int i = 0; i < n_wins; i++)
#endif
    {
        int curr_win_size = wins[i];
        int N_s = N / curr_win_size;
        double f = 0.0;
#ifdef _WIN64
        int v = 0;
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

        f_vec[i] = f / (2.0 * N_s * curr_win_size);
    }
}

//main loop for MFDCCA (computes fluctuations starting from the beginning of the array y)
void flucMFDCCAForwCompute(double *y1, double *y2, double *t, int N, int *wins, int n_wins, double *qs, int n_q, int pol_ord, double *f_vec)
{
#ifdef _WIN64
    int iq = 0;
#endif

#ifdef _WIN64
#pragma omp parallel for
    for(iq = 0; iq < n_q; iq++)
    {
        int i = 0;
        for(i = 0; i < n_wins; i++)
#else
#pragma omp parallel for collapse(2)
    for(int iq = 0; iq < n_q; iq++)
    {
        for(int i = 0; i < n_wins; i++)
#endif
        {
            double q = qs[iq];
            int curr_win_size = wins[i];
            int N_s = N / curr_win_size;
            double f = 0.0;
#ifdef _WIN64
            int v = 0;
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
                f_vec[iq * n_wins + i] = exp(f / (double)(2 * N_s));
            }
            else
            {
                f_vec[iq * n_wins + i] = pow(f / (double)N_s, 1 / (double)q);
            }
        }
    }
}

//main loop for MFDCCA (computes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the array y)
void flucMFDCCAForwBackwCompute(double *y1, double *y2, double *t, int N, int *wins, int n_wins, double *qs, int n_q, int pol_ord, double *f_vec)
{
#ifdef _WIN64
    int iq = 0;
#endif

#ifdef _WIN64
#pragma omp parallel for
    for(iq = 0; iq < n_q; iq++)
    {
        int i = 0;
        for(i = 0; i < n_wins; i++)
#else
#pragma omp parallel for collapse(2)
    for(int iq = 0; iq < n_q; iq++)
    {
        for(int i = 0; i < n_wins; i++)
#endif
        {
            double q = qs[iq];
            int curr_win_size = wins[i];
            int N_s = N / curr_win_size;
            double f = 0.0;
#ifdef _WIN64
            int v = 0;
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
                f_vec[iq * n_wins + i] = exp(f / (double)(4 * N_s));
            }
            else
            {
                f_vec[iq * n_wins + i] = pow(f / (double)(2 * N_s), 1 / (double)q);
            }
        }
    }
}
