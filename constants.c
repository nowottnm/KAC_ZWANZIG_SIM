//setzt constanten auf
#define EXTERN_VARIABLES
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include "global.h"
#include "utility.h"
#include "steppmethods.h"
#include "gnuplot.h"
#include "Verlet_alphas.h"
#include "Opencl_Vverlet.h"
#include <gsl/gsl_cblas.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



void Set_Target_Length(double * TARGET_LENGTH, double volfrac,  void (*Stepper_Method)(int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, 
					void (*deriv) (double *y, double *ans, double t,int N)),
					void (*Poti_Handle) (double *y, double *ans, double t,int N))
{	int FLAG = 1;
	*TARGET_LENGTH = 0.0;
	if (Stepper_Method == VVerlet_Step_deriv)
	{  
		 *TARGET_LENGTH = 0.0;
		printf(" No  Targets\n");
		FLAG = 0;
	}
	//if (Stepper_Method == VVerlet_Step_Target_Circle)
	//{
	//	*TARGET_LENGTH = LATTICE_SPACING * sqrt(volfrac/ M_PI);
	//	printf(" Targets set to hard circles\n");
	//} 

	if (Stepper_Method == VVerlet_Step_hardsphere_reflect)
	{  
		 *TARGET_LENGTH = LATTICE_SPACING * sqrt(volfrac/ M_PI);
		printf(" Targets set to hard circles\n");
		FLAG = 0;
	}

	if (Stepper_Method == VVerlet_Step_hardsphere_reflect_Kupf)
	{  
		 *TARGET_LENGTH = LATTICE_SPACING * sqrt(volfrac/ M_PI);
		printf(" Targets set to hard circles Kupferman \n");
		FLAG = 0;
	}

	if (Stepper_Method == VVerlet_Step_hardsphere_reflect_all)
	{  
		 *TARGET_LENGTH = LATTICE_SPACING * sqrt(volfrac/ M_PI);
		printf(" Targets set to hard circles, reflect ALL particles\n");
		FLAG = 0;
	}

	if (Stepper_Method == VVerlet_Step_deriv_Kupf)
	{  
		*TARGET_LENGTH = 0.0;
		printf(" No  Targets\n, Kupfermann bath");
		FLAG = 0;
	}

	if (Stepper_Method == VVerlet_Step_deriv_Box)
	{  
		 //*TARGET_LENGTH = LATTICE_SPACING * sqrt(volfrac/ M_PI/5.0);
		//printf(" Targets set to soft Yukawa circles in 3 per pore corner\n");
		*TARGET_LENGTH = 0;
		FLAG = 0;
		if (DIM > 1)
		{
			printf("wrong Dimension, please switch to 1d");
			exit(1);
		}
		if (volfrac > 0.0)
		{
			printf("no need for volume fractions >0 \n");
			exit(1);
		}
		printf("box set up with hard wall \n");	
	}
	if(FLAG)
	{
		printf("Kein bekannter Stepper in global, cannot set tarhet length!\n");
		exit(1);
	}
}

void Constants(){
// call to setup the bathh parameters and give rudimentary output for the run
	int i;
	OUTPUT_FLAG = 1; 
	printf ("OSSZI = %d  DIM =%d  ORDER =%d  \n", OSSZI, DIM, ORDER);
	printf ("KBOLTZ = %e \nTIME_STEPS =%f  \nTIME_END =%f  \n", KBOLTZ, TIME_STEPS, TIME_END);
	for(i=0;i<OSSZI;i++){
		coupling[i] = 0.0;
	}
	for(i=0;i<OSSZI;i++){
		ommega[i] = 0.0;
	}
	for(i=0;i<OSSZI;i++){
		massq[i] = 1.0;
	}
	for(i=0;i<ORDER;i++){
		y[i] = 1.0;
	}
	printf("ALPHA = %1.2E  GAMMA = %1.2E  mass  = %1.2E \n", ALPHA, GAMMA, mass);
	printf("NUMBER_TO_INF = %d  \nTIME_ARRIVAL =%1.2E  \n ", NUMBER_TO_INF, TIME_ARRIVAL);
	printf("Startbedingungen = "); 
	printf("%s",LABEL); 
	printf("\n");
	if (!(OMMEGA_READ_BINARY))
	{
		// random Number Preparation for Taus GSL genrator
		gsl_rng * r;
		const gsl_rng_type * T;
		gsl_rng_env_setup();
		T = gsl_rng_taus2;
		r = gsl_rng_alloc (T);
		gsl_rng_set(r, time(NULL)); // Seed with time
		// Random number end // 
	    Ommega_Setup(r);
	    gsl_rng_free (r);
	}else	// read in from ommega.dat
	{
		FILE *fp;
		fp = fopen("ommega.dat", "r");	
		for (i=0;i<OSSZI;i++)
		{
			if (fscanf(fp, "%lf", &ommega[i]) != 1) 
			{
	            printf("ERROT READING ommega.dat");
	            break;
        	}
		}
		fclose(fp);
		printf("read in ommega.dat\n");
	}
	if(VIRTUAL_FLAG) printf("Virtuelles Teilchen benutzt, ohne Badrückkopplung");
  	Gamma_Setup(coupling);
	
  	double ommega_max = ommega[cblas_idamax(OSSZI, ommega, 1)];
  	double coupling_max = coupling[cblas_idamax(OSSZI, coupling, 1)];
  	printf("MAX coupling =  %1.2e\n",coupling_max);
  	if(!(MAXSTEPSIZE_BINARY))
	{
  		MAXSTEPSIZE = (2*M_PI/ommega[cblas_idamax(OSSZI, ommega, 1)] * MAXSTEPSIZEMULT);
  		printf("MAXSTEPSIZE =  %1.2e set according to smallest period times %3.3e\n",MAXSTEPSIZE, MAXSTEPSIZEMULT);
  	}else
  	{
  		MAXSTEPSIZE = (MAXSTEPSIZEMULT / (double) OSSZI);
  		printf("MAXSTEPSIZE =  %1.2e set according to OSSZInr times Stepsize dt = %3.3e\n",MAXSTEPSIZE, MAXSTEPSIZEMULT);
  	}
  	Sum_Method =  SUMMATION_METHOD;
  	Stepper_Method = STEPPER_METHOD;
  	Poti_Handle = POTI_HANDLE;

  	int m = ipow(NUMBER_TO_INF,DIM); int n = DIM;
  	// Setze Gitter auf auf globale Matrix POSITIONS
	POSITIONS = (double **) malloc(m * sizeof(int *));
	POSITIONS[0] = (double *) malloc(m* n * sizeof(int));
	for (i=1; i<m; i++) POSITIONS[i] = POSITIONS[0] + n *i;
	vec_zero_i(lattice_position,DIM);
	printf("lattice allocated, current lattice cell set to zero\n");
	TARGET_CONTACT = 0;
	OUTPUT_FLAG = 0; // subsequent setup calls should not display text
	printf("Running on %d threads, optimally DIM = %d\n", THREADS, DIM);

	// set masses for chosen type of coupling
	if( !(KUPF_BINARY) )
	{
		for (int osi = 0; osi < OSSZI; osi ++)			
		{
			massq[osi] 		= 1.0;
		}
		printf("Set all masses équal\n");
	}else
	{
		for (int osi = 0; osi < OSSZI; osi ++)			
		{
			massq[osi] 		= coupling[osi] / ommega[osi] / ommega[osi];
			
		}
		printf("Set all masses according to k = m w^2l\n");
	}
	// Setup masses of bath for simpsons rule, split ends, odd and even parts

	if (SIMPSON_BINARY)
	{
		for (int osi = 1; osi < OSSZI/2-1; osi ++)
		{			
			massq[osi * 2] 		*= 2.0/3.0;
			massq[osi * 2 + 1] 	*= 4.0/3.0;
		}
		massq[0] 		*= 1.0/3.0;;
		massq[OSSZI-1] 	*= 1.0/3.0;

		if (KUPF_BINARY)
		{
			for (int osi = 0; osi < OSSZI; osi ++)			
			{
				coupling[osi] = massq[osi] * ommega[osi] * ommega[osi];
			}	
		}
	}
}


