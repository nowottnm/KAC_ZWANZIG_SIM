// Sammlung verwendeter Funktionen
#include "global.h"
#include "targets.h"
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <stdio.h>
#include "utility.h"
#include <string.h>
#include <time.h>
#include <omp.h>


void VVerlet_Step(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N)) {
	double temp_vec[N], abl[N];
	int i,j;
	for (j = 0; j < DIM; j++ )		// heavy particle outside parallel region
	{	

		v[j] = v[j] + 0.5 * a[j] * MAXSTEPSIZE;
		x[j] = x[j] + v[j] * MAXSTEPSIZE;
		temp_vec[j] = x[j];					//temp_vec zur abl- berechnung
		temp_vec[DIM + j] = v[j] * mass;
	}
	#pragma omp parallel for 			//parallelisiere über Osszillatoren
	for (i= 0 ; i < OSSZI; i++)
	{
		for (j = 0; j < DIM; j++ )
		{
			v[DIM + j + i*DIM] = v[DIM + j + i*DIM] + 0.5 * a[DIM + j + i*DIM] * MAXSTEPSIZE;
			x[DIM + j + i*DIM] = x[DIM + j + i*DIM] + v[DIM + j + i*DIM] * MAXSTEPSIZE;
			temp_vec[2*DIM + i * 2*DIM + j] = x[DIM + j + i*DIM];	//temp_vec zur abl- berechnung
			temp_vec[3*DIM + i * 2*DIM + j] = v[DIM + j + i*DIM] ;
		}
	}

	derivmethod(temp_vec, abl, *t, N); //bilde Ableitung zu t
										// gehe wieder in x-v Form für a
	for (j = 0; j < DIM; j++ )
	{
		a[j] = abl[DIM + j] / mass;
		v[j] = v[j] + 0.5 * a[j] * MAXSTEPSIZE;
	}

	#pragma omp parallel for 
	for (i= 0 ; i < OSSZI; i++)
	{
		for (j = 0; j < DIM; j++ )
		{
			a[DIM +j + i*DIM] = abl[3*DIM + j + i * 2*DIM]/massq[i];
			v[DIM +j + i*DIM] = v[DIM +j + i*DIM]\
									+ 0.5 * a[DIM +j + i*DIM] * MAXSTEPSIZE;
		}
	}
	//#pragma omp parallel for 
	//for (i= 0 ; i < N/2; i++)
	//{
	//	v[i] = v[i] + 0.5 * a[i] * MAXSTEPSIZE;
	//}	
	Update_Lattice_Position(lattice_position, x, LATTICE_SPACING);	
	*t = *t + MAXSTEPSIZE;	
}

void VVerlet_Step_deriv(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N)) {

	int i,j;
	SETNUMTHREADS
	double sum[DIM];
	#pragma omp parallel 			// Values for Bath particles parralized///////////////////
	{
		int ID = omp_get_thread_num();
		int MAX_THREADS = omp_get_num_threads();
		
		if (MAX_THREADS < DIM)  // Catch case : not enough threads
		{	for (j= 0; j<DIM; j++)
			{
				sum[j] = 0.0f;				// prepare sum for heavy particle
			}
		}else if (ID < DIM)  // heavy particle dimensions go to threads
		{
			sum[ID] = 0.0f;				// prepare sum for heavy particle
		}
		#pragma omp for
		for (i= 0 ; i < N/2; i++)	// Erstes Halbupdate ////
		{
			double v_temp = v[i] = v[i] + 0.5 * a[i] * MAXSTEPSIZE;		//vnext = v + a/2 *dt
			x[i] = x[i] + v_temp * MAXSTEPSIZE;			//xnext = x + vnext *dt
		}
		#pragma omp for 		// Values for Bath particles parralized
		for (i = 0 ; i < OSSZI; i++)
		{ 
			double coup = coupling[i];
			double om = ommega[i];
			for(j = 0; j < DIM; j++)
			{		// Ableitung /
				double a_temp = a[DIM + j + i*DIM] = - om * om * x[DIM + j + i*DIM]\
											  + coup * x[j]; 		// dp_q/dt = -w^2 * q+ gamma * x
				v[DIM +j + i*DIM] = v[DIM +j + i*DIM]\
										+ 0.5 * a_temp * MAXSTEPSIZE; // Zweites Halbupdate ////

			}
		}
		double private_sum[DIM];
		double private_c[DIM] = {0.0,0.0};			// Kahan correction
		double private_t[DIM] = {0.0,0.0};
		double private_y[DIM] = {0.0,0.0};
		for(j = 0; j < DIM; j++)		// Values for heavy particle
		{
			private_sum[j] = 0.0; 			// JEder thread berechnet Partialsumme
		}
		#pragma omp for 			// berechne die Summe gamma_i *q_i - gamma_i^2/ommega_i^2 * x
		for (i = 0 ; i < OSSZI; i++)
		{ 	
			double coup = coupling[i];
			double om = ommega[i];
			for (j=0; j<DIM; j++)
			{	// Kahan summation scheme
				private_y[j] = coup * x[DIM + j + i*DIM]\
						 - pow(coup,2.0)/pow(om,2.0) * x[j]  - private_c[j]; 
				private_t[j] = private_sum[j] + private_y[j];
				private_c[j] = (private_t[j] - private_sum[j]) - private_y[j];
				private_sum[j] = private_t[j];
			}
		}
        for(j=0; j<DIM; j++) 
        {
        	#pragma omp atomic		//addiere threadsummen auf
            sum[j] =sum[j] + private_sum[j];
        }
	    #pragma omp barrier
	    if ( MAX_THREADS < DIM)  // Catch case : not enough threads
		{	for (j= 0; j<DIM; j++)
			{
				a[j] = sum[j] /mass ;   	//p-Ableitung
			}
		}else if (ID < DIM)  // heavy particle dimensions go to threads
		{
			a[ID] = sum[ID] /mass ;   	//p-Ableitung
		}
	}
	
	for(j = 0; j< DIM; j++)
	{
		v[j] = v[j] + 0.5 * a[j] * MAXSTEPSIZE;
	}	
	*t = *t + MAXSTEPSIZE;	
}

void VVerlet_Step_hardsphere_reflect(const int N,  double  x[] , double v[], double a[], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N)) {

	int i,j;
	SETNUMTHREADS
	double sum[DIM];						// variable stepsize maximum MAXSTEPSIZE
	double target_vec[20];
	int reflec_flag = 0; 					// set to 1of target if root found
	int target_nr_hit = 0;					// number of actual target hit
	double contact_stepsize;
	// upper search bound
	int number_targets;
	double lower_t_global[2*DIM]; 			// array containing contact times for all targets in cell
	int reflec_flag_global[2*DIM];			// global array tro indicate reflection on target
	for (i=0; i < 2*DIM; i++)
	{
		lower_t_global[i] = MAXSTEPSIZE;	// all times a equal or smaller stepsize, so search for minimum later
		reflec_flag_global[i] = 0;

	}
	// determine time to collsion --------------------------------------------------------------------------------------------
	contact_stepsize = MAXSTEPSIZE;
	Update_Lattice_Position(lattice_position, x, LATTICE_SPACING);
	if 	(	(number_targets =  TARGET_METHOD(lattice_position, target_vec, 	// Find Targets
		LATTICE_SPACING, ESCAPE_CELLS)\
		)>0
		)
	{	
		#pragma omp parallel for
		for (int target_i = 0; target_i < number_targets; target_i++)
		{
			double lower_t; 						// lower search bound
			double upper_t;
			int signold; int signnew;
			lower_t = 0.0; 							// lower search bound
			upper_t = MAXSTEPSIZE;					// upper search bound
			double steps = 500.0;
			int root_flag = 1;						// set to one if root found and also one at fitst loop
			int root_loop = 0;						// number of current loop
			int root_loop_max = 3;
			while ((root_flag == 1)	&&	(root_loop < root_loop_max))
			{
				root_flag = 0;
				root_loop += 1;
			// get sign lower bound before loop
				double eval_value = 0.0;
				for(i = 0; i< DIM; i++)
				{
					double v_temp = v[i] + 0.5 * a[i] * lower_t;
					double x_temp = x[i] + v_temp * lower_t;
					eval_value += pow( (x_temp - target_vec[i + DIM * target_i]) , 2.0);		// f = sum (x-R)^2
				}
				eval_value = eval_value - pow(TARGET_LENGTH , 2.0);
				if (eval_value > 0)
				{
					signold = 1;
				}else
				{
					signold = -1;
				}
				double dt = (upper_t - lower_t)/ steps;
			// check all signs in intervall till first change
				for (double t_check = lower_t; t_check < upper_t; t_check += dt)
				{
					eval_value = 0.0;
					for(i = 0; i< DIM; i++)
					{
						double v_temp = v[i] + 0.5 * a[i] * t_check;
						double x_temp = x[i] + v_temp * t_check;
						eval_value += pow( (x_temp - target_vec[i + DIM * target_i]) , 2.0);		// f = sum (x-R)^2
					}
					eval_value = eval_value - pow(TARGET_LENGTH , 2.0);
					//printf("%2.2e\n", eval_value);
					if (eval_value > 0)
					{
						signnew = 1;
					}else
					{
						signnew = -1;
					}
					if (signnew != signold)				// at first sign change adjust intervall
					{
						lower_t = t_check - dt;
						steps = 10000.0;				// refine only 4 orders of magnitude per loop
						upper_t = t_check;
						root_flag = 1;					// redo loop to refine root
						lower_t_global[target_i] = lower_t;	
						break;
					}
					signold = signnew;
				}
			}
		}
	}// collision detection end--------------------------------------------------------------------------------------------------
	double minimum_t = MAXSTEPSIZE;
	for (i=0; i < number_targets; i++)				//find smallest reflection time
	{
		if(lower_t_global[i] < minimum_t)
		{
			reflec_flag = 1;
			target_nr_hit = i;
			minimum_t = lower_t_global[i];
		}
	}
	if (reflec_flag ==1 )
	{
		if (TARGET_FLAG)				// check for first contact
		{
			time_first_contact += *t + minimum_t;
			TARGET_FLAG = 0;// no need to test again after 1st contact
		}
		contact_stepsize =  minimum_t;				// pick lower boud so particle is JUST outside the boundary - Else zero search converges badly
	}else
	{
		contact_stepsize = MAXSTEPSIZE;
	}

	#pragma omp parallel 							// Values for Bath particles parralized//
	{
		int ID = omp_get_thread_num();
		int MAX_THREADS = omp_get_num_threads();
		
		if (MAX_THREADS < DIM) 						// Catch case : not enough threads
		{	for (j= 0; j<DIM; j++)
			{
				sum[j] = 0.0f;						// prepare sum for heavy particle
			}
		}else if (ID < DIM)  						// heavy particle dimensions go to threads
		{
			sum[ID] = 0.0f;							// prepare sum for heavy particle
		}
		#pragma omp for
		for (i= 0 ; i < N/2; i++)					// Erstes Halbupdate ////
		{
			double v_temp = v[i] = v[i] + 0.5 * a[i] * contact_stepsize;		//vnext = v + a/2 *dt
			x[i] = x[i] + v_temp * contact_stepsize;			//xnext = x + vnext *dt
		}
		#pragma omp for 		// Values for Bath particles parralized
		for (i = 0 ; i < OSSZI; i++)
		{ 
			double coup = coupling[i];
			double om = ommega[i];
			for(j = 0; j < DIM; j++)
			{		// Ableitung /
				double a_temp = a[DIM + j + i*DIM] = - om * om * x[DIM + j + i*DIM]\
											  + coup * x[j]; 		// dp_q/dt = -w^2 * q+ gamma * x
				v[DIM +j + i*DIM] = v[DIM +j + i*DIM]\
										+ 0.5 * a_temp * contact_stepsize; // Zweites Halbupdate ////

			}
		}
		double private_sum[DIM];
		double private_c[DIM] = {0.0,0.0};			// Kahan correction
		double private_t[DIM] = {0.0,0.0};
		double private_y[DIM] = {0.0,0.0};
		for(j = 0; j < DIM; j++)					// Values for heavy particle
		{
			private_sum[j] = 0.0; 					// ech thread calcs partial sum
		}
		#pragma omp for 							// calc sum gamma_i *q_i - gamma_i^2/ommega_i^2 * x
		for (i = 0 ; i < OSSZI; i++)
		{ 	
			double coup = coupling[i];
			double om = ommega[i];
			for (j=0; j<DIM; j++)
			{										// Kahan summation scheme
				private_y[j] = coup * x[DIM + j + i*DIM]\
						 - pow(coup,2.0)/pow(om,2.0) * x[j]  - private_c[j]; 
				private_t[j] = private_sum[j] + private_y[j];
				private_c[j] = (private_t[j] - private_sum[j]) - private_y[j];
				private_sum[j] = private_t[j];
			}
		}
        for(j=0; j<DIM; j++) 
        {
        	#pragma omp atomic		//addiere threadsummen auf
            sum[j] =sum[j] + private_sum[j];
        }
	    #pragma omp barrier
	    if ( MAX_THREADS < DIM)  // Catch case : not enough threads
		{	for (j= 0; j<DIM; j++)
			{
				a[j] = sum[j] /mass ;   	//p-Ableitung
			}
		}else if (ID < DIM)  // heavy particle dimensions go to threads
		{
			a[ID] = sum[ID] /mass ;   	//p-Ableitung
		}
	}
	
	for(j = 0; j< DIM; j++)
	{
		v[j] = v[j] + 0.5 * a[j] * contact_stepsize;
	}	
	*t = *t + contact_stepsize;	
	double norm = 0.0;
	for(j = 0; j< DIM; j++)
	{
			norm  += pow(target_vec[j + DIM * target_nr_hit] - x[j],2.0);
	}
	norm = sqrt(norm);
	if (number_targets)													
	{
		if(norm < TARGET_LENGTH * 0.98)									// catch case particle entered sphere, only nonempty cells
		{
			FLAG_DUMP_VALUE  = 1;
			printf("\nsphere entered at t=%2.2e, lattice point %d %d", *t, lattice_position[0], lattice_position[1]);
		}	
	}
	
	if ((reflec_flag))												// reflect on surface normal. Assume particle is on the surface
																	// prevent case of inner reflection in case particle still got in
	{
		double surface_normal[DIM];
		double vec_norm = 0.0;
		double scalar_prod = 0.0;
		// find vector towards target
		for(j = 0; j< DIM; j++)
		{
			surface_normal[j] = target_vec[j + DIM * target_nr_hit] - x[j];
		}
		for(j = 0; j< DIM; j++)
		{
			vec_norm += surface_normal[j] * surface_normal[j];
		}
		// point surface vector outwards and normalize
		vec_norm = sqrt(vec_norm);
		for(j = 0; j< DIM; j++)
		{
			surface_normal[j] = surface_normal[j]/vec_norm;
		}
		// find part of v in surface direction
		for(j = 0; j< DIM; j++)
		{
			scalar_prod += surface_normal[j] * v[j];
		}
		// reflect that part
		for(j = 0; j< DIM; j++)
		{
			v[j] = v[j] - 2.0 * scalar_prod * surface_normal[j];
		}

	}
}

void VVerlet_Step_hardsphere_reflect_all(const int N,  double  x[] , double v[], double a[], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N)) {

	int i,j;
	SETNUMTHREADS
	double sum[DIM];						// variable stepsize maximum MAXSTEPSIZE
	double target_vec[20];
	int reflec_flag = 0; 					// set to 1of target if root found
	int target_nr_hit = 0;					// number of actual target hit
	double contact_stepsize;
	// upper search bound
	int number_targets;
	double lower_t_global[2*DIM]; 			// array containing contact times for all targets in cell
	int reflec_flag_global[2*DIM];			// global array tro indicate reflection on target
	for (i=0; i < 2*DIM; i++)
	{
		lower_t_global[i] = MAXSTEPSIZE;	// all times a equal or smaller stepsize, so search for minimum later
		reflec_flag_global[i] = 0;

	}
	// determine time to collsion --------------------------------------------------------------------------------------------
	contact_stepsize = MAXSTEPSIZE;
	Update_Lattice_Position(lattice_position, x, LATTICE_SPACING);
	if 	(	(number_targets =  TARGET_METHOD(lattice_position, target_vec, 	// Find Targets
		LATTICE_SPACING, ESCAPE_CELLS)\
		)>0
		)
	{	
		#pragma omp parallel for
		for (int target_i = 0; target_i < number_targets; target_i++)
		{
			double lower_t; 						// lower search bound
			double upper_t;
			int signold; int signnew;
			lower_t = 0.0; 							// lower search bound
			upper_t = MAXSTEPSIZE;					// upper search bound
			double steps = 500.0;
			int root_flag = 1;						// set to one if root found and also one at fitst loop
			int root_loop = 0;						// number of current loop
			int root_loop_max = 3;
			while ((root_flag == 1)	&&	(root_loop < root_loop_max))
			{
				root_flag = 0;
				root_loop += 1;
			// get sign lower bound before loop
				double eval_value = 0.0;
				for(i = 0; i< DIM; i++)
				{
					double v_temp = v[i] + 0.5 * a[i] * lower_t;
					double x_temp = x[i] + v_temp * lower_t;
					eval_value += pow( (x_temp - target_vec[i + DIM * target_i]) , 2.0);		// f = sum (x-R)^2
				}
				eval_value = eval_value - pow(TARGET_LENGTH , 2.0);
				if (eval_value > 0)
				{
					signold = 1;
				}else
				{
					signold = -1;
				}
				double dt = (upper_t - lower_t)/ steps;
			// check all signs in intervall till first change
				for (double t_check = lower_t; t_check < upper_t; t_check += dt)
				{
					eval_value = 0.0;
					for(i = 0; i< DIM; i++)
					{
						double v_temp = v[i] + 0.5 * a[i] * t_check;
						double x_temp = x[i] + v_temp * t_check;
						eval_value += pow( (x_temp - target_vec[i + DIM * target_i]) , 2.0);		// f = sum (x-R)^2
					}
					eval_value = eval_value - pow(TARGET_LENGTH , 2.0);
					//printf("%2.2e\n", eval_value);
					if (eval_value > 0)
					{
						signnew = 1;
					}else
					{
						signnew = -1;
					}
					if (signnew != signold)				// at first sign change adjust intervall
					{
						lower_t = t_check - dt;
						steps = 10000.0;				// refine only 4 orders of magnitude per loop
						upper_t = t_check;
						root_flag = 1;					// redo loop to refine root
						lower_t_global[target_i] = lower_t;	
						break;
					}
					signold = signnew;
				}
			}
		}
	}// collision detection end--------------------------------------------------------------------------------------------------
	double minimum_t = MAXSTEPSIZE;
	for (i=0; i < number_targets; i++)				//find smallest reflection time
	{
		if(lower_t_global[i] < minimum_t)
		{
			reflec_flag = 1;
			target_nr_hit = i;
			minimum_t = lower_t_global[i];
		}
	}
	if (reflec_flag ==1 )
	{
		if (TARGET_FLAG)				// check for first contact
		{
			time_first_contact += *t + minimum_t;
			TARGET_FLAG = 0;// no need to test again after 1st contact
		}
		contact_stepsize =  minimum_t;				// pick lower boud so particle is JUST outside the boundary - Else zero search converges badly
	}else
	{
		contact_stepsize = MAXSTEPSIZE;
	}

	#pragma omp parallel 							// Values for Bath particles parralized//
	{
		int ID = omp_get_thread_num();
		int MAX_THREADS = omp_get_num_threads();
		
		if (MAX_THREADS < DIM) 						// Catch case : not enough threads
		{	for (j= 0; j<DIM; j++)
			{
				sum[j] = 0.0f;						// prepare sum for heavy particle
			}
		}else if (ID < DIM)  						// heavy particle dimensions go to threads
		{
			sum[ID] = 0.0f;							// prepare sum for heavy particle
		}
		#pragma omp for
		for (i= 0 ; i < N/2; i++)					// Erstes Halbupdate ////
		{
			double v_temp = v[i] = v[i] + 0.5 * a[i] * contact_stepsize;		//vnext = v + a/2 *dt
			x[i] = x[i] + v_temp * contact_stepsize;			//xnext = x + vnext *dt
		}
		#pragma omp for 		// Values for Bath particles parralized
		for (i = 0 ; i < OSSZI; i++)
		{ 
			double coup = coupling[i];
			double om = ommega[i];
			for(j = 0; j < DIM; j++)
			{		// Ableitung /
				double a_temp = a[DIM + j + i*DIM] = - om * om * x[DIM + j + i*DIM]\
											  + coup * x[j]; 		// dp_q/dt = -w^2 * q+ gamma * x
				v[DIM +j + i*DIM] = v[DIM +j + i*DIM]\
										+ 0.5 * a_temp * contact_stepsize; // Zweites Halbupdate ////

			}
		}
		double private_sum[DIM];
		double private_c[DIM] = {0.0,0.0};			// Kahan correction
		double private_t[DIM] = {0.0,0.0};
		double private_y[DIM] = {0.0,0.0};
		for(j = 0; j < DIM; j++)					// Values for heavy particle
		{
			private_sum[j] = 0.0; 					// ech thread calcs partial sum
		}
		#pragma omp for 							// calc sum gamma_i *q_i - gamma_i^2/ommega_i^2 * x
		for (i = 0 ; i < OSSZI; i++)
		{ 	
			double coup = coupling[i];
			double om = ommega[i];
			for (j=0; j<DIM; j++)
			{										// Kahan summation scheme
				private_y[j] = coup * x[DIM + j + i*DIM]\
						 - pow(coup,2.0)/pow(om,2.0) * x[j]  - private_c[j]; 
				private_t[j] = private_sum[j] + private_y[j];
				private_c[j] = (private_t[j] - private_sum[j]) - private_y[j];
				private_sum[j] = private_t[j];
			}
		}
        for(j=0; j<DIM; j++) 
        {
        	#pragma omp atomic		//addiere threadsummen auf
            sum[j] =sum[j] + private_sum[j];
        }
	    #pragma omp barrier
	    if ( MAX_THREADS < DIM)  // Catch case : not enough threads
		{	for (j= 0; j<DIM; j++)
			{
				a[j] = sum[j] /mass ;   	//p-Ableitung
			}
		}else if (ID < DIM)  // heavy particle dimensions go to threads
		{
			a[ID] = sum[ID] /mass ;   	//p-Ableitung
		}
	}
	
	for(j = 0; j< DIM; j++)
	{
		v[j] = v[j] + 0.5 * a[j] * contact_stepsize;
	}	
	*t = *t + contact_stepsize;	
	double norm = 0.0;
	for(j = 0; j< DIM; j++)
	{
			norm  += pow(target_vec[j + DIM * target_nr_hit] - x[j],2.0);
	}
	norm = sqrt(norm);
	if (number_targets)													
	{
		if(norm < TARGET_LENGTH * 0.98)									// catch case particle entered sphere, only nonempty cells
		{
			FLAG_DUMP_VALUE  = 1;
			printf("\nsphere entered at t=%2.2e, lattice point %d %d", *t, lattice_position[0], lattice_position[1]);
		}	
	}
	
	if ((reflec_flag))												// reflect on surface normal. Assume particle is on the surface
																	// prevent case of inner reflection in case particle still got in
	{
		double surface_normal[DIM];
		double vec_norm = 0.0;
		double scalar_prod = 0.0;
		// find vector towards target
		for(j = 0; j< DIM; j++)
		{
			surface_normal[j] = target_vec[j + DIM * target_nr_hit] - x[j];
		}
		for(j = 0; j< DIM; j++)
		{
			vec_norm += surface_normal[j] * surface_normal[j];
		}
		// point surface vector outwards and normalize
		vec_norm = sqrt(vec_norm);
		for(j = 0; j< DIM; j++)
		{
			surface_normal[j] = surface_normal[j]/vec_norm;
		}
		// find part of v in surface direction
		for(j = 0; j< DIM; j++)
		{
			scalar_prod += surface_normal[j] * v[j];
		}
		// reflect that part for ALL
		#pragma omp parallel for
		for (i = 0; i < OSSZI; i++ )
		{
			for(j = 0; j< DIM; j++)
			{
				v[i*DIM + j] = v[i*DIM + j] - 2.0 * scalar_prod * surface_normal[j];
			}	
		}
	}
}

void VVerlet_Step_deriv_Kupf(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, 
					void (*derivmethod) (double *y, double *ans, double t,int N)) {
	//mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
	// use fully mobile Kupferman bath that is coupling U = sum gamma_i (q_i-X)^2

	int i,j;
	SETNUMTHREADS
	double sum[DIM];
	#pragma omp parallel 			// Values for Bath particles parralized///////////////////
	{
		int ID = omp_get_thread_num();
		int MAX_THREADS = omp_get_num_threads();
		
		if (MAX_THREADS < DIM)  // Catch case : not enough threads
		{	for (j= 0; j<DIM; j++)
			{
				sum[j] = 0.0f;				// prepare sum for heavy particle
			}
		}else if (ID < DIM)  // heavy particle dimensions go to threads
		{
			sum[ID] = 0.0f;				// prepare sum for heavy particle
		}
		#pragma omp for
		for (i= 0 ; i < N/2; i++)	// Erstes Halbupdate ////
		{
			double v_temp = v[i] = v[i] + 0.5 * a[i] * MAXSTEPSIZE;		//vnext = v + a/2 *dt
			x[i] = x[i] + v_temp * MAXSTEPSIZE;			//xnext = x + vnext *dt
		}
		#pragma omp for 		// Values for Bath particles parralized
		for (i = 0 ; i < OSSZI; i++)
		{ 
			double coup = coupling[i];
			double om = ommega[i];
			for(j = 0; j < DIM; j++)
			{		// Ableitung /
				double a_temp = a[DIM + j + i*DIM] = -coup * (x[DIM + j + i*DIM] - x[j])/massq[i]; 	
				v[DIM +j + i*DIM] = v[DIM +j + i*DIM]\
										+ 0.5 * a_temp * MAXSTEPSIZE; // Zweites Halbupdate ////

			}
		}
		double private_sum[DIM];
		double private_c[DIM] = {0.0,0.0};			// Kahan correction
		double private_t[DIM] = {0.0,0.0};
		double private_y[DIM] = {0.0,0.0};
		for(j = 0; j < DIM; j++)		// Values for heavy particle
		{
			private_sum[j] = 0.0; 			// JEder thread berechnet Partialsumme
		}
		#pragma omp for 			// berechne die Summe gamma_i *q_i - gamma_i^2/ommega_i^2 * x
		for (i = 0 ; i < OSSZI; i++)
		{ 	
			double coup = coupling[i];
			double om = ommega[i];
			for (j=0; j<DIM; j++)
			{	// Kahan summation scheme
				private_y[j] = coup * (x[DIM + j + i*DIM] - x[j]) - private_c[j]; 
				private_t[j] = private_sum[j] + private_y[j];
				private_c[j] = (private_t[j] - private_sum[j]) - private_y[j];
				private_sum[j] = private_t[j];
			}
		}
        for(j=0; j<DIM; j++) 
        {
        	#pragma omp atomic		//addiere threadsummen auf
            sum[j] =sum[j] + private_sum[j];
        }
	    #pragma omp barrier
	    if ( MAX_THREADS < DIM)  // Catch case : not enough threads
		{	for (j= 0; j<DIM; j++)
			{
				a[j] = sum[j] /mass ;   	//p-Ableitung
			}
		}else if (ID < DIM)  // heavy particle dimensions go to threads
		{
			a[ID] = sum[ID] /mass ;   	//p-Ableitung
		}
	}
	
	for(j = 0; j< DIM; j++)
	{
		v[j] = v[j] + 0.5 * a[j] * MAXSTEPSIZE;
		VIRTUAL_V[j] = v[j];
		VIRTUAL_X[j] = x[j];
	}	
	*t = *t + MAXSTEPSIZE;	
}

void VVerlet_Step_hardsphere_reflect_Kupf(const int N,  double  x[] , double v[], double a[], double *t, //mach einzelnen update Schritt nach Velocity-Verlet ohne Hinderniss
					void (*derivmethod) (double *y, double *ans, double t,int N)) {

	int i,j;
	SETNUMTHREADS
	double sum[DIM];						// variable stepsize maximum MAXSTEPSIZE
	double target_vec[20];
	int reflec_flag = 0; 					// set to 1of target if root found
	int target_nr_hit = 0;					// number of actual target hit
	double contact_stepsize = LASTSTEPSIZE;
	// upper search bound
	int number_targets;
	double lower_t_global[2*DIM]; 			// array containing contact times for all targets in cell
	int reflec_flag_global[2*DIM];			// global array tro indicate reflection on target
	for (i=0; i < 2*DIM; i++)
	{
		lower_t_global[i] = LASTSTEPSIZE;	// all times a equal or smaller stepsize, so search for minimum later
		reflec_flag_global[i] = 0;

	}
	// determine time to collsion --------------------------------------------------------------------------------------------
	contact_stepsize = LASTSTEPSIZE;
	if(VIRTUAL_FLAG)
	{
		Update_Lattice_Position(lattice_position, VIRTUAL_X, LATTICE_SPACING);
	}else
	{
		Update_Lattice_Position(lattice_position, x, LATTICE_SPACING);
	}
	
	if 	(	(number_targets =  TARGET_METHOD(lattice_position, target_vec, 	// Find Targets
		LATTICE_SPACING, ESCAPE_CELLS)\
		)>0
		)
	{
		double N_CHECK = 1000.0;		
		#pragma omp parallel for
		for (int target_i = 0; target_i < number_targets; target_i++)	// check for all targets in increments of LASTSTEPSIZE/N_CHECK for contact
		{
			double dh = LASTSTEPSIZE /  N_CHECK;
			for (double h_check = 0.0; h_check < LASTSTEPSIZE; h_check += dh)
			{
				double r_check = 0.0;
				for (int j = 0; j < DIM; j++)
				{
					double v_temp = v[j] + 0.5 * a[j] * h_check;
					double x_temp = x[j] + v_temp * h_check;
					if(VIRTUAL_FLAG)
					{
						v_temp = VIRTUAL_V[j] + 0.5 * a[j] * h_check;
						x_temp = VIRTUAL_X[j] + v_temp * h_check;
					} 
					r_check += pow( (x_temp - target_vec[j + target_i * DIM] ) ,2.0 );
				}
				r_check = sqrt(r_check);
				if (r_check < TARGET_LENGTH)
				{
					reflec_flag = 1;
					target_nr_hit = target_i;
					contact_stepsize = h_check - dh;
					if (TARGET_FLAG)				// check for first contact
					{
						time_first_contact += *t + contact_stepsize;
						TARGET_FLAG = 0;// no need to test again after 1st contact
					}
					break;
					break;
				} 
			}	
		}
	}
	//-------------------------------------------------------------------------------------------------------------------
	#pragma omp parallel 							// Values for Bath particles parralized//
	{
		int ID = omp_get_thread_num();
		int MAX_THREADS = omp_get_num_threads();
		
		if (MAX_THREADS < DIM) 						// Catch case : not enough threads
		{	for (j= 0; j<DIM; j++)
			{
				sum[j] = 0.0f;						// prepare sum for heavy particle
			}
		}else if (ID < DIM)  						// heavy particle dimensions go to threads
		{
			sum[ID] = 0.0f;							// prepare sum for heavy particle
		}
		#pragma omp for
		for (i= 0 ; i < N/2; i++)					// Erstes Halbupdate ////
		{
			double v_temp = v[i] = v[i] + 0.5 * a[i] * contact_stepsize;		//vnext = v + a/2 *dt
			x[i] = x[i] + v_temp * contact_stepsize;			//xnext = x + vnext *dt
		}
		if(VIRTUAL_FLAG)
		{
			#pragma omp for
			for (i= 0 ; i < DIM; i++)					// Erstes Halbupdate ////
			{
				double v_temp = VIRTUAL_V[i] = VIRTUAL_V[i] + 0.5 * a[i] * contact_stepsize;		//vnext = v + a/2 *dt
				VIRTUAL_X[i] = VIRTUAL_X[i] + v_temp * contact_stepsize;			//xnext = x + vnext *dt
			}
		}
		#pragma omp for 		// Values for Bath particles parralized
		for (i = 0 ; i < OSSZI; i++)
		{ 
			double coup = coupling[i];
			double om = ommega[i];
			for(j = 0; j < DIM; j++)
			{		// Ableitung /
				double a_temp = a[DIM + j + i*DIM] = -coup * (x[DIM + j + i*DIM] - x[j])/massq[i]; 		// dp_q/dt = -w^2 * q+ gamma * x
				v[DIM +j + i*DIM] = v[DIM +j + i*DIM]\
										+ 0.5 * a_temp * contact_stepsize; // Zweites Halbupdate ////

			}
		}
		double private_sum[DIM];
		double private_c[DIM] = {0.0,0.0};			// Kahan correction
		double private_t[DIM] = {0.0,0.0};
		double private_y[DIM] = {0.0,0.0};
		for(j = 0; j < DIM; j++)					// Values for heavy particle
		{
			private_sum[j] = 0.0; 					// ech thread calcs partial sum
		}
		#pragma omp for 							// calc sum gamma_i *q_i - gamma_i^2/ommega_i^2 * x
		for (i = 0 ; i < OSSZI; i++)
		{ 	
			double coup = coupling[i];
			double om = ommega[i];
			for (j=0; j<DIM; j++)
			{										// Kahan summation scheme
				private_y[j] = coup * (x[DIM + j + i*DIM] - x[j]) - private_c[j]; 
				private_t[j] = private_sum[j] + private_y[j];
				private_c[j] = (private_t[j] - private_sum[j]) - private_y[j];
				private_sum[j] = private_t[j];
			}
		}
        for(j=0; j<DIM; j++) 
        {
        	#pragma omp atomic		//addiere threadsummen auf
            sum[j] =sum[j] + private_sum[j];
        }
	    #pragma omp barrier
	    if ( MAX_THREADS < DIM)  // Catch case : not enough threads
		{	for (j= 0; j<DIM; j++)
			{
				a[j] = sum[j] /mass ;   	//p-Ableitung
			}
		}else if (ID < DIM)  // heavy particle dimensions go to threads
		{
			a[ID] = sum[ID] /mass ;   	//p-Ableitung
		}
	}
	
	for(j = 0; j< DIM; j++)
	{
		v[j] = v[j] + 0.5 * a[j] * contact_stepsize;
	}
	if(VIRTUAL_FLAG)
	{
		for(j = 0; j< DIM; j++)
		{
			VIRTUAL_V[j] = VIRTUAL_V[j] + 0.5 * a[j] * contact_stepsize;
		}
	}	
	*t = *t + contact_stepsize;	
	double norm = 0.0;
	for(j = 0; j< DIM; j++)
	{
		if(VIRTUAL_FLAG)
		{
			norm  += pow(target_vec[j + DIM * target_nr_hit] - VIRTUAL_X[j],2.0);
		}else
		{
			norm  += pow(target_vec[j + DIM * target_nr_hit] - x[j],2.0);
		}	
	}
	norm = sqrt(norm);
	if (number_targets)													
	{
		if(norm < TARGET_LENGTH * 0.98)									// catch case particle entered sphere, only nonempty cells
		{
			FLAG_DUMP_VALUE  = 1;
			printf("\nsphere entered at t=%2.2e, lattice point %d %d", *t, lattice_position[0], lattice_position[1]);
		}	
	}
	
	if ((reflec_flag))												// reflect on surface normal. Assume particle is on the surface
																	// prevent case of inner reflection in case particle still got in
	{
		double surface_normal[DIM];
		double vec_norm = 0.0;
		double scalar_prod = 0.0;
		// find vector towards target
		for(j = 0; j< DIM; j++)
		{

			if(VIRTUAL_FLAG)
			{
				surface_normal[j] = target_vec[j + DIM * target_nr_hit] - VIRTUAL_X[j];
			}else
			{
				surface_normal[j] = target_vec[j + DIM * target_nr_hit] - x[j];
			}
		}
		for(j = 0; j< DIM; j++)
		{
			vec_norm += surface_normal[j] * surface_normal[j];
		}
		// point surface vector outwards and normalize
		vec_norm = sqrt(vec_norm);
		for(j = 0; j< DIM; j++)
		{
			surface_normal[j] = surface_normal[j]/vec_norm;
		}
		// find part of v in surface direction
		for(j = 0; j< DIM; j++)
		{
			if(VIRTUAL_FLAG)
			{
				scalar_prod += surface_normal[j] * VIRTUAL_V[j];
			}else
			{
				scalar_prod += surface_normal[j] * v[j];
			}
		}
		// reflect that part
		for(j = 0; j< DIM; j++)
		{
			if(VIRTUAL_FLAG)
			{
				VIRTUAL_V[j] = VIRTUAL_V[j] - 2.0 * scalar_prod * surface_normal[j];
			}else
			{
				v[j] = v[j] - 2.0 * scalar_prod * surface_normal[j];
			}
		}
	}
	if ((reflec_flag))
	{
		LASTSTEPSIZE = LASTSTEPSIZE - contact_stepsize; 		// do remaining step time
	}
	else
	{
		LASTSTEPSIZE = MAXSTEPSIZE;
	}
}

void VVerlet_Step_deriv_Box(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, 
					void (*derivmethod) (double *y, double *ans, double t,int N)) {
	//macht einzelnen update schritt in Box der Länge Lattice_Spacting
	// use fully mobile Kupferman bath that is coupling U = sum gamma_i (q_i-X)^2

	int i,j;
	SETNUMTHREADS
	double sum = 0.0;
	double contact_stepsize = LASTSTEPSIZE;
	int reflec_flag = 0;
	// find contact time ------------------------------------
	double aa,b,c;
	aa=b=c=0;
	
	int REFLEC_WALL = 1;
	if (REFLEC_WALL)
	{											// debugging purposes, free evolution if false
		for (int i_bound = -1; i_bound < 2; i_bound += 2)
		{
			double root1 = 0.0; double  root2 = 0.0;
			if(VIRTUAL_FLAG)
			{
				c = VIRTUAL_X[0] + i_bound * LATTICE_SPACING/2.0;					
				b = VIRTUAL_V[0];
				aa = 0.5 * a[0];
			}else
			{
				c = x[0] + i_bound * LATTICE_SPACING/2.0;					
				b = v[0];
				aa = 0.5 * a[0];
			}
			
			if(b >= 0)													//  minimize rounoff err by case split
			{
				root1 = (-b - sqrt(b*b - 4.0 * aa * c)) / (2.0 * aa);
				root2 = 2.0 * c /(-b - sqrt(b*b - 4.0 * aa * c) );
			}else
			{
				root1 =  2.0 * c /(-b + sqrt(b*b - 4.0 * aa * c) );
				root2 = (-b + sqrt(b*b - 4.0 * aa * c)) / (2.0 * aa);
			}
			if( (root1 > 0) && (root1 < contact_stepsize) )
			{
				contact_stepsize = 0.995 * root1;
				reflec_flag = 1;
			}
			if( (root2 > 0) && (root2 < contact_stepsize) )
			{
				contact_stepsize = 0.995 * root2;
				reflec_flag = 1;
			}
		}
	}
	// contact time end ||||||||||||||||||||||||||||||||||||<

	#pragma omp parallel for
	for (i= 0 ; i < N/2; i++)	// Erstes Halbupdate ////
	{
		double v_temp = v[i] = v[i] + 0.5 * a[i] * contact_stepsize;		//vnext = v + a/2 *dt
		x[i] = x[i] + v_temp * contact_stepsize;							//xnext = x + vnext *dt

	}
	if(VIRTUAL_FLAG)
	{
		double v_temp = VIRTUAL_V[0] = VIRTUAL_V[0] + 0.5 * a[0] * contact_stepsize;	// virtual particle
		VIRTUAL_X[0] = VIRTUAL_X[0] + v_temp * contact_stepsize;	
	}
	for (i = 0 ; i < OSSZI; i++)
	{ 
		double coup = coupling[i];
		double om = ommega[i];
		double a_temp = a[1 + i] = -coup * (x[1 + i] - x[0])/massq[i]; 	
		v[1 + i] = v[1 + i]\
					+ 0.5 * a_temp * contact_stepsize; // Zweites Halbupdate ////
	}
	#pragma omp parallel for reduction(+:sum) 			// berechne die Summe gamma_i *q_i - gamma_i^2/ommega_i^2 * x
	for (i = 0 ; i < OSSZI; i++)
	{ 
		double coup = coupling[i];	
		sum += coup * (x[1 + i] - x[0]); 
	}
	a[0] = sum /mass;
	v[0] = v[0] + 0.5 * a[0] * contact_stepsize;
	if(VIRTUAL_FLAG)
	{
		VIRTUAL_V[0] = VIRTUAL_V[0] + 0.5 * a[0] * contact_stepsize;	
	}
	*t = *t + contact_stepsize;
	if(reflec_flag)
	{
		LASTSTEPSIZE = LASTSTEPSIZE - contact_stepsize;
		if(REFLECTALL_FLAG)								// reflect whole bath if flag is set
		{
			for (i = 0 ; i < OSSZI; i++)
			{ 
				double v_rel = v[1 + i] - v[0];
				v[1 + i] = - v[0] + v_rel;						//keep relative velocity after reflection of heavy particle the same
			}
		}
		if(VIRTUAL_FLAG)
		{
			VIRTUAL_V[0] = -VIRTUAL_V[0];										// reflect heavy
		}else
		{
			v[0] = -v[0];										// reflect heavy
		}
	}
	else
	{
		LASTSTEPSIZE = MAXSTEPSIZE;
	}	
}

void VVerlet_Step_deriv_Box_Ramp(const int N,  double  x[N/2] , double v[N/2], double a[N/2], double *t, 
					void (*derivmethod) (double *y, double *ans, double t,int N)) {
	//macht einzelnen update schritt in Box der Länge Lattice_Spacting
	// use fully mobile Kupferman bath that is coupling U = sum gamma_i (q_i-X)^2

	int i,j;
	SETNUMTHREADS
	double sum = 0.0;
	double contact_stepsize = LASTSTEPSIZE;
	int reflec_flag = 0;
	


	#pragma omp parallel for
	for (i= 0 ; i < N/2; i++)	// Erstes Halbupdate ////
	{
		double v_temp = v[i] = v[i] + 0.5 * a[i] * contact_stepsize;		//vnext = v + a/2 *dt
		x[i] = x[i] + v_temp * contact_stepsize;							//xnext = x + vnext *dt
	}
	for (i = 0 ; i < OSSZI; i++)
	{ 
		double coup = coupling[i];
		double om = ommega[i];
		double a_temp = a[1 + i] = -coup * (x[1 + i] - x[0])/massq[i]; 	
		v[1 + i] = v[1 + i]\
					+ 0.5 * a_temp * contact_stepsize; // Zweites Halbupdate ////
	}
	#pragma omp parallel for reduction(+:sum) 			// berechne die Summe gamma_i *q_i - gamma_i^2/ommega_i^2 * x
	for (i = 0 ; i < OSSZI; i++)
	{ 
		double coup = coupling[i];	
		sum += coup * (x[1 + i] - x[0]); 
	}
	a[0] = sum /mass;
	
	if(x[0] < -(LATTICE_SPACING/2.0 - L_BALL) )			//left wall;
	{
		double height_ramp = 2*KBOLTZ*TEMP;
		a[0] += height_ramp/L_BALL/mass;
	}
	if(x[0] > (LATTICE_SPACING/2.0 - L_BALL) )			//right wall;
	{
		double height_ramp = 2*KBOLTZ*TEMP;
		a[0] += -height_ramp/L_BALL/mass;
	}
	v[0] = v[0] + 0.5 * a[0] * contact_stepsize;
	
	*t = *t + contact_stepsize;	
}

