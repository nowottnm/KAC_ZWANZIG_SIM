// sammlung Verlet steps

#ifndef HEADER_FILE_STEP
#define HEADER_FILE_STEP

void VVerlet_Step(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N));

void VVerlet_Step_deriv(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N));

void VVerlet_Step_hardsphere_reflect(const int N,  double  x[] , double v[], double a[], double *t, 	//mach einzelnen update Schritt nach Velocity-Verlets
					void (*derivmethod) (double *y, double *ans, double t,int N));					// stop at spehere, reflect


void VVerlet_Step_hardsphere_reflect_all(const int N,  double  x[] , double v[], double a[], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N));

void VVerlet_Step_deriv_Kupf(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N));

void VVerlet_Step_hardsphere_reflect_Kupf(const int N,  double  x[] , double v[], double a[], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N));
void VVerlet_Step_deriv_Box(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, 
					void (*derivmethod) (double *y, double *ans, double t,int N));
void VVerlet_Step_deriv_Box_Ramp(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, 
					void (*derivmethod) (double *y, double *ans, double t,int N));

#endif
