//    cLoops.h - header of C loops of fathon package
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


extern void flucDFAForwCompute(double *y, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec);
extern void flucDFAForwBackwCompute(double *y, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec);
extern void flucMFDFAForwCompute(double *y, double *t, int N, int *wins, int n_wins, double *qs, int n_q, int pol_ord, double *f_vec);
extern void flucMFDFAForwBackwCompute(double *y, double *t, int N, int *wins, int n_wins, double *qs, int n_q, int pol_ord, double *f_vec);
extern void flucDCCAAbsCompute(double *y1, double *y2, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec);
extern void flucDCCANoAbsCompute(double *y1, double *y2, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec);
extern double HTCompute(double *y, double *t, int scale, int N, int pol_ord, int v);
extern void flucDCCAForwAbsComputeNoOverlap(double *y1, double *y2, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec);
extern void flucDCCAForwBackwAbsComputeNoOverlap(double *y1, double *y2, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec);
extern void flucDCCAForwNoAbsComputeNoOverlap(double *y1, double *y2, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec);
extern void flucDCCAForwBackwNoAbsComputeNoOverlap(double *y1, double *y2, double *t, int N, int *wins, int n_wins, int pol_ord, double *f_vec);
extern void flucMFDCCAForwCompute(double *y1, double *y2, double *t, int N, int *wins, int n_wins, double *qs, int n_q, int pol_ord, double *f_vec);
extern void flucMFDCCAForwBackwCompute(double *y1, double *y2, double *t, int N, int *wins, int n_wins, double *qs, int n_q, int pol_ord, double *f_vec);
