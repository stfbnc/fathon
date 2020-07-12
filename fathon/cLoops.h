//    cLoops.h - header of C loops of fathon package
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

extern double flucDFAForwCompute(double *y, int curr_win_size, int N, int pol_ord);
extern double flucDFAForwBackwCompute(double *y, int curr_win_size, int N, int pol_ord);
extern double flucMFDFAForwCompute(double *y, int curr_win_size, double q, int N, int pol_ord);
extern double flucMFDFAForwBackwCompute(double *y, int curr_win_size, double q, int N, int pol_ord);
extern double flucDCCAAbsCompute(double *y1, double *y2, int curr_win_size, int N, int pol_ord);
extern double flucDCCANoAbsCompute(double *y1, double *y2, int curr_win_size, int N, int pol_ord);
extern double HTCompute(double *y, int scale, int N, int pol_ord, int v);
