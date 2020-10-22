//    cFuncs.h - array and mathematical operations C functions of fathon package
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multifit.h>

//polynomial fit
void polynomialFit(int obs, int degree, double *dx, double *dy, double *store)
{
    gsl_matrix *X = gsl_matrix_alloc(obs, degree);
    gsl_vector *y = gsl_vector_alloc(obs);
    gsl_vector *c = gsl_vector_alloc(degree);
    gsl_matrix *cov = gsl_matrix_alloc(degree, degree);

    for(int i = 0; i < obs; i++)
    {
        for(int j = 0; j < degree; j++)
        {
            gsl_matrix_set(X, i, j, pow(dx[i], j));
        }
        gsl_vector_set(y, i, dy[i]);
    }

    double chisq;
    gsl_multifit_linear_workspace *ws = gsl_multifit_linear_alloc(obs, degree);
    gsl_multifit_linear(X, y, c, cov, &chisq, ws);

    for(int i = 0; i < degree; i++)
    {
        store[i] = gsl_vector_get(c, i);
    }

    gsl_multifit_linear_free(ws);
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(c);
}
