//    cLoops.c - C loops of fathon package
//    Copyright (C) 2019  Stefano Bianchi
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

//main loop for DFA (computes fluctuations starting from the beginning of the array y)
double flucDFAForwCompute(double *y, int curr_win_size, int N, int pol_ord)
{		
	int N_s, start_lim, end_lim;
    double *t, *t_fit, *y_fit, *fit_coeffs;
	int v, j, k;
	double f = 0.0;

	fit_coeffs = malloc((pol_ord+1) * sizeof(double));
	t = malloc(N * sizeof(double));
	double_range(t, N, 1.0, 1.0);
    t_fit = malloc(curr_win_size * sizeof(double));
    y_fit = malloc(curr_win_size * sizeof(double));
    N_s = N / curr_win_size;
    for(v = 0; v < N_s; v++){
        start_lim = v * curr_win_size;
        end_lim = (v + 1) * curr_win_size - 1;
        slice_vec(t, t_fit, start_lim, end_lim);
        slice_vec(y, y_fit, start_lim, end_lim);
		polynomialFit(curr_win_size, pol_ord+1, t_fit, y_fit, fit_coeffs);
        for(j = 0; j < curr_win_size; j++){
			for(k = 0; k < pol_ord+1; k++){
				y_fit[j] -= fit_coeffs[k] * pow(t_fit[j], k);
			}
			f += pow(y_fit[j], 2.0);
		}
	}
    f = sqrt(f / (N_s*curr_win_size));
	free(fit_coeffs);
    free(t);
    free(t_fit);
    free(y_fit);
    return f;
}

//main loop for DFA (computes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the array y)
double flucDFAForwBackwCompute(double *y, int curr_win_size, int N, int pol_ord)
{		
	int N_s, start_lim, end_lim;
    double *t, *t_fit, *y_fit, *fit_coeffs;
	int v, j, k;
	double f = 0.0;

    fit_coeffs = malloc((pol_ord+1) * sizeof(double));
	t = malloc(N * sizeof(double));
	double_range(t, N, 1.0, 1.0);
    t_fit = malloc(curr_win_size * sizeof(double));
    y_fit = malloc(curr_win_size * sizeof(double));
    N_s = N / curr_win_size;
    for(v = 0; v < N_s; v++){
        start_lim = v * curr_win_size;
        end_lim = (v + 1) * curr_win_size - 1;
        slice_vec(t, t_fit, start_lim, end_lim);
        slice_vec(y, y_fit, start_lim, end_lim);
        polynomialFit(curr_win_size, pol_ord+1, t_fit, y_fit, fit_coeffs);
        for(j = 0; j < curr_win_size; j++){
            for(k = 0; k < pol_ord+1; k++){
				y_fit[j] -= fit_coeffs[k] * pow(t_fit[j], k);
			}
            f += pow(y_fit[j], 2.0);
		}
		start_lim = v * curr_win_size + (N - N_s * curr_win_size);
		end_lim = (v + 1) * curr_win_size + (N - N_s * curr_win_size) - 1;
		slice_vec(t, t_fit, start_lim, end_lim);
		slice_vec(y, y_fit, start_lim, end_lim);
		polynomialFit(curr_win_size, pol_ord+1, t_fit, y_fit, fit_coeffs);
		for(j = 0; j < curr_win_size; j++){
			for(k = 0; k < pol_ord+1; k++){
				y_fit[j] -= fit_coeffs[k] * pow(t_fit[j], k);
			}
			f += pow(y_fit[j], 2.0);
		}
	}
	f = sqrt(f/(2.0*N_s*curr_win_size));
    free(fit_coeffs);
    free(t);
    free(t_fit);
    free(y_fit);
    return f;
}

//main loop for MFDFA (computes fluctuations starting from the beginning of the array y)
double flucMFDFAForwCompute(double *y, int curr_win_size, double q, int N, int pol_ord)
{		
	int N_s, start_lim, end_lim;
	double rms;
    double *t, *t_fit, *y_fit, *fit_coeffs;
	int v, j, k;
	double f = 0.0;

    fit_coeffs = malloc((pol_ord+1) * sizeof(double));
	t = malloc(N * sizeof(double));
	double_range(t, N, 1.0, 1.0);
    t_fit = malloc(curr_win_size * sizeof(double));
    y_fit = malloc(curr_win_size * sizeof(double));
    N_s = N / curr_win_size;
    for(v = 0; v < N_s; v++){
		rms = 0.0;
        start_lim = v * curr_win_size;
        end_lim = (v + 1) * curr_win_size - 1;
        slice_vec(t, t_fit, start_lim, end_lim);
        slice_vec(y, y_fit, start_lim, end_lim);
        polynomialFit(curr_win_size, pol_ord+1, t_fit, y_fit, fit_coeffs);
        for(j = 0; j < curr_win_size; j++){
            for(k = 0; k < pol_ord+1; k++){
				y_fit[j] -= fit_coeffs[k] * pow(t_fit[j], k);
			}
            rms += pow(y_fit[j], 2.0);
		}
		if(q == 0.0){
			f += log(rms/(double)curr_win_size);
		}else{
			f += pow(rms/(double)curr_win_size, 0.5*q);
		}
	}
	if(q == 0.0){
		f = exp(f/(double)(2*N_s));
	}else{
		f = pow(f/(double)N_s, 1/(double)q);
	}
    free(fit_coeffs);
    free(t);
    free(t_fit);
    free(y_fit);
    return f;
}

//main loop for MFDFA (computes fluctuations starting from the beginning of the array y
//and then computes fluctuations again starting from the end of the array y)
double flucMFDFAForwBackwCompute(double *y, int curr_win_size, double q, int N, int pol_ord)
{		
	int N_s, start_lim, end_lim;
	double rms1, rms2;
    double *t, *t_fit, *y_fit, *fit_coeffs;
	int v, j, k;
	double f = 0.0;

    fit_coeffs = malloc((pol_ord+1) * sizeof(double));
	t = malloc(N * sizeof(double));
	double_range(t, N, 1.0, 1.0);
    t_fit = malloc(curr_win_size * sizeof(double));
    y_fit = malloc(curr_win_size * sizeof(double));
    N_s = N / curr_win_size;
    for(v = 0; v < N_s; v++){
		rms1 = 0.0;
		rms2 = 0.0;
        start_lim = v * curr_win_size;
        end_lim = (v + 1) * curr_win_size - 1;
        slice_vec(t, t_fit, start_lim, end_lim);
        slice_vec(y, y_fit, start_lim, end_lim);
        polynomialFit(curr_win_size, pol_ord+1, t_fit, y_fit, fit_coeffs);
        for(j = 0; j < curr_win_size; j++){
            for(k = 0; k < pol_ord+1; k++){
				y_fit[j] -= fit_coeffs[k] * pow(t_fit[j], k);
			}
            rms1 += pow(y_fit[j], 2.0);
		}
		start_lim = v * curr_win_size + (N - N_s * curr_win_size);
		end_lim = (v + 1) * curr_win_size + (N - N_s * curr_win_size) - 1;
		slice_vec(t, t_fit, start_lim, end_lim);
		slice_vec(y, y_fit, start_lim, end_lim);
		polynomialFit(curr_win_size, pol_ord+1, t_fit, y_fit, fit_coeffs);
		for(j = 0; j < curr_win_size; j++){
    		for(k = 0; k < pol_ord+1; k++){
				y_fit[j] -= fit_coeffs[k] * pow(t_fit[j], k);
			}
			rms2 += pow(y_fit[j], 2.0);
		}
		if(q == 0.0){
			f += (log(rms1/(double)curr_win_size) + log(rms2/(double)curr_win_size));
		}else{
			f += (pow(rms1/(double)curr_win_size, 0.5*q) + pow(rms2/(double)curr_win_size, 0.5*q));
		}
	}
	if(q == 0.0){
		f = exp(f/(double)(4*N_s));
	}else{
		f = pow(f/(double)(2*N_s), 1/(double)q);
	}
    free(fit_coeffs);
    free(t);
    free(t_fit);
    free(y_fit);
    return f;
}

//main loop for DCCA (computes fluctuations using absolute values)
double flucDCCAAbsCompute(double *y1, double *y2, int curr_win_size, int N, int pol_ord)
{
	int N_s, end_lim;
    double *t, *t_fit, *y_fit1, *y_fit2, *fit_coeffs1, *fit_coeffs2;
	int v, j, k;
	double f = 0.0;
	
	fit_coeffs1 = malloc((pol_ord+1) * sizeof(double));
	fit_coeffs2 = malloc((pol_ord+1) * sizeof(double));
    t = malloc(N * sizeof(double));
    double_range(t, N, 1.0, 1.0);
    t_fit = malloc((curr_win_size+1) * sizeof(double));
    y_fit1 = malloc((curr_win_size+1) * sizeof(double));
    y_fit2 = malloc((curr_win_size+1) * sizeof(double));
    N_s = N - curr_win_size;
    for(v = 0; v < N_s; v++){
        end_lim = v + curr_win_size;
        slice_vec(t, t_fit, v, end_lim);
        slice_vec(y1, y_fit1, v, end_lim);
        slice_vec(y2, y_fit2, v, end_lim);
        polynomialFit(curr_win_size+1, pol_ord+1, t_fit, y_fit1, fit_coeffs1);
        polynomialFit(curr_win_size+1, pol_ord+1, t_fit, y_fit2, fit_coeffs2);
		for(j = 0; j <= curr_win_size; j++){
    		for(k = 0; k < pol_ord+1; k++){
				y_fit1[j] -= fit_coeffs1[k] * pow(t_fit[j], k);
				y_fit2[j] -= fit_coeffs2[k] * pow(t_fit[j], k);
			}
			f += fabs(y_fit1[j] * y_fit2[j]);
        }
    }
    f = sqrt(f / (N_s*(curr_win_size-1)));
    free(fit_coeffs1);
    free(fit_coeffs2);
    free(t);
    free(t_fit);
    free(y_fit1);
    free(y_fit2);
    return f;
}

//main loop for DCCA (computes fluctuations without using absolute values)
double flucDCCANoAbsCompute(double *y1, double *y2, int curr_win_size, int N, int pol_ord)
{
	int N_s, end_lim;
    double *t, *t_fit, *y_fit1, *y_fit2, *fit_coeffs1, *fit_coeffs2;
	int v, j, k;
	double f = 0.0;
	
	fit_coeffs1 = malloc((pol_ord+1) * sizeof(double));
	fit_coeffs2 = malloc((pol_ord+1) * sizeof(double));
    t = malloc(N * sizeof(double));
    double_range(t, N, 1.0, 1.0);
    t_fit = malloc((curr_win_size+1) * sizeof(double));
    y_fit1 = malloc((curr_win_size+1) * sizeof(double));
    y_fit2 = malloc((curr_win_size+1) * sizeof(double));
    N_s = N - curr_win_size;
    for(v = 0; v < N_s; v++){
        end_lim = v + curr_win_size;
        slice_vec(t, t_fit, v, end_lim);
        slice_vec(y1, y_fit1, v, end_lim);
        slice_vec(y2, y_fit2, v, end_lim);
        polynomialFit(curr_win_size+1, pol_ord+1, t_fit, y_fit1, fit_coeffs1);
        polynomialFit(curr_win_size+1, pol_ord+1, t_fit, y_fit2, fit_coeffs2);
		for(j = 0; j <= curr_win_size; j++){
    		for(k = 0; k < pol_ord+1; k++){
				y_fit1[j] -= fit_coeffs1[k] * pow(t_fit[j], k);
				y_fit2[j] -= fit_coeffs2[k] * pow(t_fit[j], k);
			}
			f += y_fit1[j] * y_fit2[j];
		}
    }
    f = f / (N_s*(curr_win_size-1));
    free(fit_coeffs1);
    free(fit_coeffs2);
    free(t);
    free(t_fit);
    free(y_fit1);
    free(y_fit2);
    return f;
}

//main loop for HT (computes fluctuations)
double HTCompute(double *y, int scale, int N, int pol_ord, int v)
{
	int end_lim;
	double *t, *t_fit, *y_fit, *fit_coeffs;
	int j, k;
	double f = 0.0;

    fit_coeffs = malloc((pol_ord+1) * sizeof(double));
	t = malloc(N * sizeof(double));
	double_range(t, N, 1.0, 1.0);
    t_fit = malloc(scale * sizeof(double));
    y_fit = malloc(scale * sizeof(double));
	//for(v = 0; v < (N-scale+1); v++){
	end_lim = v + scale - 1;
	slice_vec(t, t_fit, v, end_lim);
	slice_vec(y, y_fit, v, end_lim);
	polynomialFit(scale, pol_ord+1, t_fit, y_fit, fit_coeffs);
	for(j = 0; j < scale; j++){
    		for(k = 0; k < pol_ord+1; k++){
			y_fit[j] -= fit_coeffs[k] * pow(t_fit[j], k);
		}
		f += pow(y_fit[j], 2.0);
	}
	f = sqrt(f/(double)scale);
	free(fit_coeffs);
    free(t);
    free(t_fit);
    free(y_fit);
    return f;
}
