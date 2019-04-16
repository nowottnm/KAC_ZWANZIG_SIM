// sammlung verwendeter Funktionen
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#ifndef HEADER_FILE_UTIL
#define HEADER_FILE_UTIL


double Norm_Diff(double *vec1,double *vec2, int dim, int p); //berechnet norm von || vec1-vec2 ||_p
void Lattice_Setup(double **positions,int  numbertoinf,int dim, double latticespacing); // Legt Koordianten auf mxn Gittermatrix der dimension numbertoinf**dim, dim
int ipow(int base, int exp);
int Kahan_Sum(int N, double input[], double *ans);
int Sign_Sum(int N, double input[], double *ans);
int Neumaier_Sum(int N, double input[N], double *ans);
void deleteSpaces(char src[], char dst[]);

void vec_zero(double *v,int  N);
void vec_zero_i(int *v,int  N);

int Update_Lattice_Position(int * pos, const double *vec, const double lattice_spacing);

// begin stepper with different external forces -------//

void VVerlet_Step(int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, 
					void (*deriv) (double *y, double *ans, double t,int N));

void VVerlet_Step_Target_Square(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, 
					void (*deriv) (double *y, double *ans, double t,int N));

void VVerlet_Step_deriv(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N));

void VVerlet_Step_Yukawa(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N));

void VVerlet_Step_Pore_Rectangle(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, // wie Verlet step, aber Teilchen reflektiert stehen falls in Target (Quadrat)
					void (*derivmethod) (double *y, double *ans, double t,int N));


void VVerlet_Step_Pore_Yukawa(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N));

// end stepper with different external forces -------//

void VVerlet(int N, double *y, double **ans, double *t,
			void (*deriv) (double *y, double *ans, double t, int N));

void VVerlet_parallel(const int N, const double *y, double **ans, double *t,							
				void (*derivmethod) (double *y, double *ans, double t, int N));

void VVerlet_parallel_burn(const int N, const  double *y, double **ans, double *t,							// Velocity Verlet für Start y, Ausgabe ans[LengthT][ORDER] zu Zeiten T
				void (*derivmethod) (double *y, double *ans, double t, int N));

void Bath_Setup(double *y, char *label, gsl_rng * r);

void Ommega_Setup(gsl_rng * r);

void Gamma_Setup(double * gamma);

void deriv(double *y, double *ans, double t, int N);
void deriv_parallel(const double *yin, double *ans, double t, int N);

void deriv_Yukawa(const double *y, double *ans, double t, int N);

#endif