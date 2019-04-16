#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include "utility.h"
#include "gnuplot.h"
#include "constants.h"
#include <gsl/gsl_cblas.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_pow_int.h>

#ifdef EXTERN_VARIABLES
#define EXTERN
#else
#define EXTERN extern
#endif



#define OSSZI 14000
#define	SIMULATIONEN 10
#define KUPF_BINARY	1				// set to 1 for Kupferman/ fully mobile bath
#define VIRTUAL_FLAG	0			// set to 1 to transfer force one way onto virtual particle, no bath recouoling
#define SIMPSON_BINARY 1
#define REFLECTALL_FLAG 0
#define OMMEGA_READ_BINARY 0		// if true, don't regenerate ommega. use ommega.dat
#define DIM 1
#define ORDER ((2 * OSSZI + 2) * DIM)
#define LABELPREFIX "GIT_HUB"	//a string to add to your result folder
#define T_POWERS 0 							// no function, just compiler needs to know
#define T_MULT 1.0              			// same
#define LOGSCALE_BINARY 0       			// leave at 0, toggles time setup. WIP
#define TIME_STEPS  0.05					// snapshot intervall for saved values
#define LENGTH_T 	1000  +  9 				// Number N of saved steps (format must be N +9)
#define TIME_SHIFT 500
#define TIME_END (double)((double) (LENGTH_T-9) * TIME_STEPS) //leave macro as is
//#endif
#define THREADSNR 1 						// most efficient at dimension
#define	KBOLTZ  1.380641e-8
#define	TEMP  	100.0
#define	GAMMA  1.0 							// drag
#define LOWER_CUTOFF_FLAG 0 				// 1 for lower cutoff frequency when settinh up the ommega
#define	ALPHA  0.8							// exponent free random walk msqdplcmnt t^alpha
#define	mass  1.0 							// mass tracer particle
#define	TIME_ARRIVAL 3000.0 				// Diffusion time over lattice cell
#define MAXSTEPSIZE_BINARY	0 				// set to 1 to decouple from smallest period to fixed N * DT = MAXSTEPSIEMULT
#define MAXSTEPSIZEMULT    5.0e-4			// stepsize in terms of shortest period scale
#define	NUMBER_TO_INF  (500)				// muss gerade sein, hindernisse die Symetrisch um null verteilt sind in 1-D
#define	VOL_FRAC  0.0
#define LABEL "ZWANZIG_THERMALIZED"			// "ZWANZIG" or "ZWANZIG_THERMALISED"
#define SUMMATION_METHOD Neumaier_Sum  		// Sign_, Kahan oder Neumaier_ gefolgt von sum
#define STEPPER_METHOD VVerlet_Step_deriv_Box	//bis jetzt VVerlet_Step_Target_Square oder VVerlet_Step_deriv oder VVerlet_Step_Yukawa Verlet_Step_Target_Pore_Rectangle
#define POTI_HANDLE deriv_parallel
#define NORM_INT 0  						// für Norm_Diff die art der p-norm, 0 für infty
#define MAX_NR_PLOTS 50
#define ESCAPE_CELLS  0                		// anzahl der gitterzellen ohne hinderniss -1
#define REFLEC_BINARY 1 					// toggle bewteewn reflecting new(1) or old(0) velocity, or no velcotiy adjust (2)
#define YUKAWA_EXPONENT_FACTOR (double) 1.0	// factor a in V(r) = exp(- a r/d) * A /r^n
#define YUKAWA_POWER (double) 1.0			// factor n in V(r) = exp(- a r/d) * A /r^n
#define YUKAWA_TOL	(double) 1.0			// Yukawa wall set so thermal energy equals Potential resistance times TOL at target Length
#define BURN_INT 	0 						// number of short samples to approx DIFF_COEFF
#define BURN_END 300.0						// determine approx DIFF_COEFF in timespan (0, BURN_END)
#define BURN_START 100.0 					// determine start time for burn 
#define DIFF_ANALYTIC_BINARY 1 				// set 1 to use analytic instead of approximate coeccicient
#define FLAG_RAND_BATH 0 					// set 1 for unifrom random bath, 0 for uniform deteriministic

											//(so as to avoid sampling when particle is still equilibrating)
#define PORE_OPENING 0.01 					// Size of Pore Opening relativ to LATTICE SPACING
#define PORE_THICKNESS 0.1 
#define DENSITY_SAMPLING_LENGTH 10
#define DENSITY_SAMPLING_POINTS 10
#define DIFF_COEFF_MODIFY 0.01 				// modify for VVerlet_Step_drag_sphere
#define DENSITY_SNAPSHOTS 5
#define TARGET_METHOD	Get_Lattice_Targets_Corner_Square // method to determine free cells
#define PARALLEL  //nicht definieren falls ohne -fopenmp compiliert!
#ifdef PARALLEL
	#define THREADS THREADSNR
	#define SETNUMTHREADS omp_set_num_threads(THREADSNR);
#else
	#define SETNUMTHREADS
	#define THREADS 1  //wichtig zur zeitberechnung dass dann THREADS auf 1 gesetzt wird, wenn nicht parallel operiert
#endif

//#define OCL
#ifdef OCL
	#define VVERLET Opencl_VVerlet
#else
	#define VVERLET VVerlet_parallel
#endif
EXTERN int FLAG_DUMP_VALUE;			// if error occurs set 1, dump calculation, start fresh
EXTERN int NUMBER_PORE_OBSTACLES;  	//number of spheres in Yukawa pore wall per side
EXTERN double LATTICE_SPACING;
EXTERN double DIFF_COEFF;
EXTERN double TARGET_LENGTH;
EXTERN double VIRTUAL_X[DIM];
EXTERN double VIRTUAL_V[DIM];
EXTERN int OUTPUT_FLAG;				//set 1 before first output, after that 0 to hide messages
									// from Setuo functions for ommega and gamma
EXTERN int TARGET_CONTACT;			// counts number of hits on lattice
EXTERN int TARGET_FLAG;				// Toggles after first hit per run
EXTERN int lattice_position[DIM];	// aktuelle Gitterzelle	
EXTERN double massq[OSSZI];
EXTERN double time_first_contact;
EXTERN double coupling[OSSZI];  //Koppl. Konst. Gamma im Zwanzig
EXTERN double ommega[OSSZI]; 	// Osszi Kreisfreq
EXTERN double y[ORDER];			// Startvec
EXTERN double MAXSTEPSIZE;		//Schrittweite Integrator
EXTERN double LASTSTEPSIZE;
EXTERN int **POSITIONS; 	//enthält Koordinated der Gitterpunkte
EXTERN double L_BALL;		// ballistische Länge der Diffusion
EXTERN int (*Sum_Method)(int N, double input[], double *ans);
EXTERN void (*Poti_Handle)(const double *y, double *ans, double t,int N); 			// externe oder Hindernisskraft
EXTERN void (*Stepper_Method)(int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, 
					void (*deriv) (double *y, double *ans, double t,int N));			
					// Methode mit oder ohne Hindernisse für Verlet_schritt
EXTERN gsl_rng * RAND_GLOBAL;	// globally accessible random variable

