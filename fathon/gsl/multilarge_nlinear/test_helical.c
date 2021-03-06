#define helical_N         3
#define helical_P         3

static double helical_x0[helical_P] = { -1.0, 0.0, 0.0 };
static double helical_x[helical_P] = { 1.0, 0.0, 0.0 };
static double helical_epsrel = 1.0e-12;

static double helical_J[helical_N * helical_P];

static void
helical_checksol(const double x[], const double sumsq,
                 const double epsrel, const char *sname,
                 const char *pname)
{
  size_t i;
  const double sumsq_exact = 0.0;

  gsl_test_rel(sumsq, sumsq_exact, epsrel, "%s/%s sumsq",
               sname, pname);

  for (i = 0; i < helical_P; ++i)
    {
      gsl_test_rel(x[i], helical_x[i], epsrel, "%s/%s i=%zu",
                   sname, pname, i);
    }
}

static int
helical_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double x3 = gsl_vector_get(x, 2);
  double theta = (x1 >= 0.0) ? 0.0 : 5.0;
  double nx = gsl_hypot(x1, x2);

  gsl_vector_set(f, 0, 10.0 * (x3 - 5.0/M_PI*atan(x2 / x1) - theta));
  gsl_vector_set(f, 1, 10.0*(nx - 1.0));
  gsl_vector_set(f, 2, x3);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
helical_df (CBLAS_TRANSPOSE_t TransJ, const gsl_vector * x,
            const gsl_vector * u, void * params, gsl_vector * v,
            gsl_matrix * JTJ)
{
  gsl_matrix_view J = gsl_matrix_view_array(helical_J, helical_N, helical_P);
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double nx = gsl_hypot(x1, x2);
  double nx_sq = nx * nx;
  double term1 = 50.0 / (M_PI * nx_sq);
  double term2 = 10.0 / nx;

  gsl_matrix_set(&J.matrix, 0, 0, term1*x2);
  gsl_matrix_set(&J.matrix, 0, 1, -term1*x1);
  gsl_matrix_set(&J.matrix, 0, 2, 10.0);

  gsl_matrix_set(&J.matrix, 1, 0, term2*x1);
  gsl_matrix_set(&J.matrix, 1, 1, term2*x2);
  gsl_matrix_set(&J.matrix, 1, 2, 0.0);

  gsl_matrix_set(&J.matrix, 2, 0, 0.0);
  gsl_matrix_set(&J.matrix, 2, 1, 0.0);
  gsl_matrix_set(&J.matrix, 2, 2, 1.0);

  if (v)
    gsl_blas_dgemv(TransJ, 1.0, &J.matrix, u, 0.0, v);

  if (JTJ)
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &J.matrix, 0.0, JTJ);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static int
helical_fvv (const gsl_vector * x, const gsl_vector * v,
             void *params, gsl_vector * fvv)
{
  double x1 = gsl_vector_get(x, 0);
  double x2 = gsl_vector_get(x, 1);
  double v1 = gsl_vector_get(v, 0);
  double v2 = gsl_vector_get(v, 1);
  double term1 = v2*x1 - v1*x2;
  double term2 = v1*x1 + v2*x2;
  double term3 = x1*x1 + x2*x2;

  gsl_vector_set(fvv, 0, 100.0 / M_PI * (term1 / term3) * (term2 / term3));
  gsl_vector_set(fvv, 1, 10.0 * (term1 * term1) / pow(term3, 1.5));
  gsl_vector_set(fvv, 2, 0.0);

  (void)params; /* avoid unused parameter warning */

  return GSL_SUCCESS;
}

static gsl_multilarge_nlinear_fdf helical_func =
{
  helical_f,
  helical_df,
  helical_fvv,
  helical_N,
  helical_P,
  NULL,
  0,
  0,
  0,
  0
};

static test_fdf_problem helical_problem =
{
  "helical",
  helical_x0,
  NULL,
  &helical_epsrel,
  &helical_checksol,
  &helical_func
};
