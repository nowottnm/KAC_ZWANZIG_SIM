// Sammlung verwendeter Funktionen
#include "global.h"
#include "steppmethods.h"
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








double Norm_Diff(double *vec1,double *vec2, int dim, int p) //berechnet norm von || vec1-vec2 ||_p
{
	double result = 0.0d;
	double temp;
	if (p == 1)
	{
		for (int i = 0; i < DIM; i++)
		{
			result += fabs(vec1[i] - vec2[i]); 
		}
		return(result);	
	}
	if (p == 2)
	{
		for (int i = 0; i < DIM; i++)
		{
			result += pow(vec1[i] - vec2[i],2.0); 
		}
		return(sqrt(result));	
	}
	if (p == 0)
	{
		for (int i = 0; i < DIM; i++)
		{
			temp = fabs(vec1[i] - vec2[i]);
			if (result < temp)
			{
				result = temp;
			} 
		}
		return(result);	
	}
	return (0.0f);
}


void Lattice_Setup(double **positions,int  numbertoinf,int dim, double latticespacing) // puts coordinates into lattice matrix of dimension numbertoinf**dim, dim
{
	printf("setze Gitter auf \n");
	int m = -numbertoinf + 1;
	int l = -numbertoinf + 1;
	int count = 0;
	while(count<ipow(numbertoinf,dim))
	{
		if (m> (numbertoinf - 1) ) // if over line, restart
		{
			m = -numbertoinf +1;
			l = l + 2; 
		}
		positions[count][0] = m * latticespacing/2;
		positions[count][1] = l * latticespacing/2;
		m = m + 2;
		count = count +1;
	}
	printf("Gitter aufgesetzt \n");
}

int ipow(int base, int exp) // I^p for integer i
{
	int result = 1;
	if(		(base < 0)\
		||	(exp < 0))
	{
		printf("wrong base or exponent, smaller zero!\n");
		exit;
	}
	while (exp)
	{
		if (exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
	}
	return result;
}

int Kahan_Sum(int N, double input[], double *ans){  // Add elements via Kahan summation, code is pretty much Wikipedia copy paste
	int i;
	double y,t,sum = 0.0d;
    double c = 0.0d;                 // A running compensation for lost low-order bits.
    for (i=0; i<N; i++){
        y = input[i] - c;	    // So far, so good: c is zero.
        t = sum + y;         	// Alas, sum is big, y small, so low-order digits of y are lost.
        c = (t - sum) - y;       // (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
        sum = t;                 // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
    }                      		// Next time around, the lost low part will be added to y in a fresh attempt.
    *ans  = (sum);
    return 0;
}

int Sign_Sum(int N, double *input, double *ans){  // Add seperated by sign
	int i;
	double sum = 0.0, negsum = 0.0;
	for (i=0; i<N; i++){
		if(input[i]<0.0) negsum += input[i];
		if(input[i]>=0.0) sum += input[i];
	}                      		
	*ans = sum + negsum;
	return 0;
}

int Neumaier_Sum(int N, double input[N], double *ans){  // Add elements of input by Neumaier summation
	int i;
	double t,sum = input[0];
    double c = 0.0;                 // A running compensation for lost low-order bits.
    for (i=1; i<N; i++){
    	t = sum + input[i];
    	if (fabs(sum)>=fabs(input[i]))
    	{
        	c += (sum - t) + input[i]; // If sum is bigger, low-order digits of input[i] are lost.	
        } else
        {
            c += (input[i] - t) + sum; // Else low-order digits of sum are lost
        }
        sum = t;
    }                      // Next time around, the lost low part will be added to y in a fresh attempt. 
    *ans = (sum + c);         // Correction only applied once in the very end
    return 0;
}


void deleteSpaces(char src[], char dst[]){      // deletes spaces from string
   // src is supposed to be zero ended
   // dst is supposed to be large enough to hold src
	int s, d=0;
	for (s=0; src[s] != 0; s++)
		if (src[s] != ' ') 
		{
			dst[d] = src[s];
			if (dst[d] == ':') 
			{
				dst[d] = '_';
			}
			d++;
		}
		
		dst[d] = 0;
	}


void vec_zero(double *v,int  N){ // zero target vector of double
	int i;
	for (i= 0 ; i < N; i++)
	{
		v[i] = 0.0;	
	}	

}
void vec_zero_i(int *v,int  N){ // zero target vector of ints
	int i;
	for (i= 0 ; i < N; i++)
	{
		v[i] = 0;	
	}	

}
void zero_Impuls(double *v)
{ 								// zero target vector of momenta
	double P[OSSZI + 1];
	double sum;
	double OS = OSSZI * 1.0;
	for (int j = 0; j < DIM; j++)
	{	
		for (int i = 0; i < OSSZI; i++)
		{
			P[i] = v[DIM + i * DIM + j]; 
		}
		P[OSSZI] = v[j] * mass;	
		Sum_Method(OSSZI+1,P,&sum);
		for (int i = 0; i < OSSZI; i++)
		{
			v[DIM + i * DIM + j] -= sum/( OS  + 1.0); 
		}
		v[j] -= sum/( OS  + 1.0) /mass;
	}	
}

int Update_Lattice_Position(int * pos, const  double *vec, const double lattice_spacing){ // updatet Int-Vector pos when leaving square lattice cell
	
	int ret_val = 0;			// return zero if no update neccessary, -1 if out of bounds
	for (int j = 0; j< DIM; j++)
	{
		int trigger = 1;
		while(trigger)
		{
			if (vec[j]< pos[j]*lattice_spacing - lattice_spacing/2.0)	
			{
				pos[j]-=1;
				trigger = 1;
				ret_val = 1;
			}else 
			{
				trigger = 0;
			}
			if (vec[j]> pos[j]*lattice_spacing + lattice_spacing/2.0)
			{
				pos[j]+=1;
				trigger = 1;
				ret_val = 1;
			}else
			{
				trigger = 0;
			}
		}
	}
	if ((abs(pos[0])> NUMBER_TO_INF) || (abs(pos[1])> NUMBER_TO_INF))
	{
		printf(" \n particle out of bounds");
		FLAG_DUMP_VALUE = 1;											// error occured
		return(-1);
	}
	return(ret_val);
}


void VVerlet(int N, double *y, double **ans, double *t,							// Velocity Verlet for Start y, return ans[LengthT][ORDER] to times T
				void (*derivmethod) (double *y, double *ans, double t, int N)){
	int i,j;
	for(i=0; i< N; i++){
	ans[0][i] = y[i];  				// give initial vec to answer
	}
	double x[N/2], v[N/2], a[N/2], z[N];
	vec_zero(x,N/2);
	vec_zero(v,N/2);
	vec_zero(a,N/2);
	vec_zero(z,N);
	int current_index = 0;
	double current_t = t[current_index];
	derivmethod(y, z, current_t, N); 				// get first derivative							
	for (j = 0; j < DIM; j++ )						// bring y-p  Form to x-v-a
	{
		x[j] = y[j];
		v[j] = y[DIM + j]/mass;
		a[j] = z[DIM + j]/mass;
		#pragma omp parallel for
		for (i= 0 ; i < OSSZI; i++)
		{
					x[DIM + j + i*DIM] = y[2*DIM + j + i * 2*DIM];
					v[DIM + j + i*DIM] = y[3*DIM + j + i * 2*DIM]/massq[i];
					a[DIM + j + i*DIM] = z[3*DIM + j + i * 2*DIM]/massq[i];
		}
	}

	while(current_index < LENGTH_T-1)
	{
		while(fabs(current_t - t[current_index + 1])> MAXSTEPSIZE)
		{ 	//integrate till space till next time smaller than stepsize
			Stepper_Method(N, x, v, a, &current_t, derivmethod);
		}
												// get to y-Form
		current_index = current_index +1;
		for (j = 0; j < DIM; j++ )
		{
			ans[current_index][DIM + j] = mass * v[j];
			ans[current_index][j] 	    = x[j];
			#pragma omp parallel for 
			for (i= 0 ; i < OSSZI; i++)
			{
						ans[current_index][2*DIM + j + i * 2*DIM] = x[DIM +j+  i*DIM];
						ans[current_index][3*DIM + j + i * 2*DIM] = v[DIM +j+  i*DIM] * massq[i]; // übergib SChritt an antwort
			}
		}		
	}

} 

void VVerlet_parallel(const int N, const  double *y, double **ans, double *t,							// Velocity Verlet for Start y, return ans[LengthT][ORDER] to times T
				void (*derivmethod) (double *y, double *ans, double t, int N)){
	int i,j;
	for(i=0; i< N; i++){
	ans[0][i] = y[i];  				// give initial vec to answer
	}
	TARGET_FLAG = 1;
	double x[N/2], v[N/2], a[N/2], z[N];
	vec_zero(x,N/2);
	vec_zero(v,N/2);
	vec_zero(a,N/2);
	vec_zero(z,N);
	int current_index = 0;
	double current_t = t[current_index];
	derivmethod(y, z, current_t, N); 						// giv f'
	int steps_in_intervall = 0;						// counts steps done in intervall  t, t+ dt
	int max_steps_in_intervall = 100000;			// after X steps, assume error occured
	for (j = 0; j < DIM; j++ )						//  y-p  --> x-v-a form
	{
		x[j] = y[j];
		VIRTUAL_X[j] = x[j];
		v[j] = y[DIM + j]/mass;
		VIRTUAL_V[j] = v[j];
		a[j] = z[DIM + j]/mass;
	}

	#pragma omp parallel for 
	for (i= 0 ; i < OSSZI; i++)
	{
		for (j = 0; j < DIM; j++ )
		{
			x[DIM + j + i*DIM] = y[2*DIM + j + i * 2*DIM];
			v[DIM + j + i*DIM] = y[3*DIM + j + i * 2*DIM]/massq[i];
			a[DIM + j + i*DIM] = z[3*DIM + j + i * 2*DIM]/massq[i];
		}
	}
	LASTSTEPSIZE = MAXSTEPSIZE;
	while(current_index < LENGTH_T-1)
	{
		while(fabs(current_t - t[current_index + 1])> MAXSTEPSIZE)			//integrate till space till next time smaller than stepsize
		{ 	
			Stepper_Method(N, x, v, a, &current_t, derivmethod);
			steps_in_intervall++;
			if(steps_in_intervall > max_steps_in_intervall)
			{
				FLAG_DUMP_VALUE  = 1;										// too many steps, error occured
				printf("too many steps taken to advance one time point \n");
			}
			if (FLAG_DUMP_VALUE) 											// if error occured somewhere leave calculation
			{
				return;
			}
		}
		steps_in_intervall = 0; 				// reset stepcount								
		current_index = current_index +1;		// get  back to y-form
		for (j = 0; j < DIM; j++ )
		{
			if(VIRTUAL_FLAG)
			{
				ans[current_index][DIM + j] = mass * VIRTUAL_V[j];
				ans[current_index][j] 	    = VIRTUAL_X[j];
			}else
			{
				ans[current_index][DIM + j] = mass * v[j];
				ans[current_index][j] 	    = x[j];
			}
		}
		#pragma omp parallel for 
		for (i= 0 ; i < OSSZI; i++)
		{
			for (j = 0; j < DIM; j++ )
			{
				ans[current_index][2*DIM + j + i * 2*DIM] = x[DIM +j+  i*DIM];
				ans[current_index][3*DIM + j + i * 2*DIM] = v[DIM +j+  i*DIM] * massq[i] ; // give step to answer
			}
		}		
	}
	if (TARGET_FLAG)		// if no contact, give simulation length as contact time
	{
		time_first_contact += TIME_END;
		TARGET_FLAG = 0;
	}
} 

void VVerlet_parallel_burn(const int N, const  double *y, double **ans, double *t,							// Velocity Verlet für Start y, Ausgabe ans[LengthT][ORDER] zu Zeiten T
				void (*derivmethod) (double *y, double *ans, double t, int N))
// always involes free particle for burning samples
{
	int i,j;
	for(i=0; i< N; i++){
	ans[0][i] = y[i];  				// übergib startvektor an antwort
	}
	TARGET_FLAG = 1;
	double x[N/2], v[N/2], a[N/2], z[N];
	vec_zero(x,N/2);
	vec_zero(v,N/2);
	vec_zero(a,N/2);
	vec_zero(z,N);
	int current_index = 0;
	double current_t = t[current_index];
	derivmethod(y, z, current_t, N); 				// übergib erste Ableitung
	// bringe aus y-p y Form auf x-v-a
	for (j = 0; j < DIM; j++ )
	{
		x[j] = y[j];
		v[j] = y[DIM + j]/mass;
		a[j] = z[DIM + j]/mass;
	}
	#pragma omp parallel for 
	for (i= 0 ; i < OSSZI; i++)
	{
		for (j = 0; j < DIM; j++ )
		{
			x[DIM + j + i*DIM] = y[2*DIM + j + i * 2*DIM];
			v[DIM + j + i*DIM] = y[3*DIM + j + i * 2*DIM]/massq[i];
			a[DIM + j + i*DIM] = z[3*DIM + j + i * 2*DIM]/massq[i];
		}
	}
	

	while(current_index < LENGTH_T-1)
	{
		while(fabs(current_t - t[current_index + 1])> MAXSTEPSIZE)
		{ 	//integriere bis abstand zu  nächstem t kleiner schrittweite
			VVerlet_Step_deriv(N, x, v, a, &current_t, derivmethod);
			//printf ("%fx      %ft  \n", x[1], current_t);
		}
												// rechne wieder auf y-Form
		current_index = current_index +1;
		for (j = 0; j < DIM; j++ )
		{
			ans[current_index][DIM + j] = mass * v[j];
			ans[current_index][j] 	    = x[j];
		}
		#pragma omp parallel for 
		for (i= 0 ; i < OSSZI; i++)
		{
			for (j = 0; j < DIM; j++ )
			{
				ans[current_index][2*DIM + j + i * 2*DIM] = x[DIM +j+  i*DIM];
				ans[current_index][3*DIM + j + i * 2*DIM] = v[DIM +j+  i*DIM] * massq[i]; // übergib SChritt an antwort
			}
		}		
	}
	if (TARGET_FLAG)		// FAlls kein erstkontakt addiere gesamtzeit
	{
		time_first_contact += TIME_END;
		TARGET_FLAG = 0;
	}
} 





void Bath_Setup(double *y, char *inputlabel, gsl_rng * r){
	int a,b,i;
	double u,P[OSSZI], sum=0.0;

	if (strcmp(inputlabel,"ZWANZIG")==0)
	{
		for (b = 0; b < DIM; b++)
		{
			y[b] = 0.0;					//massives particle seperatly
			y[DIM+b] = 0.0;
			for (a = 0; a < OSSZI; a++)
			{
				u = gsl_ran_gaussian(r, 1.0);
				y[2*DIM + a * 2*DIM + b] = sqrt(KBOLTZ * TEMP / massq[a]) / ommega[a] * u; 	//q = (k*T/m)^0.5/ommega * N(0,1)
				u = gsl_ran_gaussian(r, 1.0);
				y[3*DIM + a * 2*DIM + b] = sqrt(KBOLTZ * TEMP * massq[a]) * u;			//p_q = (k*T*m)^0.5 * N(0,1)
			}
		}
		return;	
	}
	if (strcmp(inputlabel,"ZWANZIG_THERMALIZED")==0)
	{
		for (b = 0; b < DIM; b++)
		{
			y[b] = 0.0;					//massives particle seperatly
			u = gsl_ran_gaussian(r, 1.0);
			y[DIM+b] =  sqrt(KBOLTZ * TEMP * mass) * u;			//p_q = (k*T*m)^0.5 * N(0,1)
			for (a = 0; a < OSSZI; a++)
			{
				u = gsl_ran_gaussian(r, 1.0);
				y[2*DIM + a * 2*DIM + b] = sqrt(KBOLTZ * TEMP / massq[a]) / ommega[a] * u; 	//q = (k*T/m)^0.5/ommega * N(0,1)
				u = gsl_ran_gaussian(r, 1.0);
				y[3*DIM + a * 2*DIM + b] = sqrt(KBOLTZ * TEMP * massq[a]) * u;			//p_q = (k*T*m)^0.5 * N(0,1)
			}
		}
		return;	
	}
	printf("FALSCHES BADLABEL");
	exit(1);
}

void Ommega_Setup(gsl_rng * r)
{

 															// draw uniformly nonrandom
	double nr_osszillators = 1.0 * OSSZI; 					// cast Osszi to double
	double higher_cutoff = 1.0/(0.5 * TIME_STEPS); 			//cleanly resolve smallest steps
	double lower_cutoff = 1.0/(1.6 * TIME_END); 			//cleanly resolve till end time


	higher_cutoff = 5.0* pow( (GAMMA/mass), 1.0/(2.0-ALPHA));
	lower_cutoff = 5E-5; 
	double ommega_step = (higher_cutoff - lower_cutoff)/(nr_osszillators);			
	for (int i = 0; i < OSSZI; i++)
	{
		ommega[i] = lower_cutoff + i * ommega_step;
	}
	printf("\nnonrandom bath ");
	printf ("Ommega_Max =   %1.3e  ", ommega[OSSZI-1] );
	printf ("Ommega_Min =   %1.3e\n", ommega[0] );
}

void Gamma_Setup(double *gamma){ 						//Set coupling for Exponent Alpha so that MSD s^2 prop to
														//t^ALPHA,  ALPHA = 1 (classical diff) seperately
														// put constants in gamma
	int i;
	double TOL = 1e-3;
	double nr_osz = (double) OSSZI;
	double delta_ommega = ommega[1]-ommega[0];
	if (! (KUPF_BINARY))
	{
		if (fabs(ALPHA-1.0)<TOL){  // Alpha = 1 
		for(i = 0; i < OSSZI; i++)
		{
			gamma[i] = sqrt(GAMMA) * ommega[i] * sqrt(delta_ommega);
		}
		if (OUTPUT_FLAG)
		{
			printf("ALPHA rounded to 1 \n");
		}
		}
		if (fabs(ALPHA-1.0)>TOL){
			for(i = 0; i < OSSZI; i++)
			{
				gamma[i] = sqrt(GAMMA) * pow(ommega[i], 0.5 + ALPHA * 0.5)\
				* sqrt( 2/sqrt(M_PI * 2) * sin(M_PI * ALPHA/2.0) * gsl_sf_gamma(1 - ALPHA)) * sqrt(delta_ommega);
			}
		}
		if (OUTPUT_FLAG)
		{
			printf("Konstanten aufgesetzt für Alpha =  %f \n " , ALPHA);
		}
	}else
	{
		if (fabs(ALPHA-1.0)<TOL)	// Alpha = 1 
		{  
			for(i = 0; i < OSSZI; i++)
			{
			
				gamma[i] = 1.0 / sqrt(2.0 * M_PI) * GAMMA * (delta_ommega);
			}
			printf("Konstanten aufgesetzt für KUPF_A rounded to 1  \n " );
		}else
		{
			for(i = 0; i < OSSZI; i++)
			{
			
				gamma[i] = 1.0 / sqrt(2.0 * M_PI) * GAMMA  * gsl_sf_gamma(1 - ALPHA)\
				* sin(M_PI * ALPHA/2.0) * pow(ommega[i], ALPHA - 1.0) * (delta_ommega);
			}
			printf("Konstanten aufgesetzt für KUPF_A = %2.2f  \n ", ALPHA );
		}

		
	}
}


void deriv(double *yin, double *ans, double t, int N){ 		// Hamiltonian derivative
															
	int i,j;
	double P[OSSZI], Q[OSSZI], temp_vec[OSSZI];
	double sum;

	for(j = 0; j < DIM; j++)
	{
		#pragma omp parallel for 
		for(i = 0; i < OSSZI; i++)
		{
			Q[i] = yin[(2 + i * 2)  * DIM  + j ];  //get P und Q for easier indexing
			P[i] = yin[(3 + i * 2)  * DIM  + j ];	
		}
		sum = 0.0; 
		#pragma omp parallel for
		for (i = 0 ; i < OSSZI; i++)
		{ // Sum
			if(KUPF_BINARY)
			{
				temp_vec[i] = -coupling[i] * (Q[i] - yin[j]);  
		
			}else
			{
				temp_vec[i] = coupling[i] * Q[i] - pow(coupling[i],2.0)/pow(ommega[i],2.0) * yin[j];  
			}
			
		}
		Sum_Method(OSSZI, temp_vec, &sum);
		ans[j] = yin[DIM+ j]/mass;  // dx/dx = p/m
		ans[DIM+j] = sum  ;   	//p-deriv
		#pragma omp parallel for
		for (i = 0 ; i < OSSZI; i++)
		{ 
			ans[(2 + i * 2)  * DIM  + j ] = P[i]/massq[i];      		// dq/dt = p_q/m_q
			if(KUPF_BINARY)
			{
				ans[(3 + i * 2)  * DIM  + j ] =  coupling[i] * (Q[i] - yin[j]); 		// dp_q/dt = -w^2 * q+ gamma * x
		
			}else
			{
				ans[(3 + i * 2)  * DIM  + j ] = - ommega[i] * ommega[i] * Q[i]  + coupling[i] * yin[j]; 		// dp_q/dt = -w^2 * q+ gamma * x
		
			}
			
		}
	}
}

void deriv_parallel(const double *yin, double *ans, double t, int N)		// deriv with i-j loops reversed for parralelization
{ 		// Ableitungsfunktion
	int i,j;
	double temp_vec[OSSZI];
	double sum;

	#pragma omp parallel for 		// Values for Bath particles parralized
	for (i = 0 ; i < OSSZI; i++)
	{ 
		for(j = 0; j < DIM; j++)
		{
			ans[(2 + i * 2)  * DIM  + j ] = yin[(3 + i * 2)  * DIM  + j ]/massq[i];      		// dq/dt = p_q/m_q
			if(KUPF_BINARY)
			{
				ans[(3 + i * 2)  * DIM  + j ] =  -coupling[i] *  (yin[(2 + i * 2)  * DIM  + j] - yin[j]); 		// dp_q/dt = -w^2 * q+ gamma * x
			}else
			{
				ans[(3 + i * 2)  * DIM  + j ] = - ommega[i] * ommega[i] * yin[(2 + i * 2)  * DIM  + j]\
										  + coupling[i] * yin[j]; 		// dp_q/dt = -w^2 * q+ gamma * x
			}
		}
	}

	for(j = 0; j < DIM; j++)		// Values for heavy particle
	{
		
		sum = 0.0; 
		#pragma omp parallel for 
		for (i = 0 ; i < OSSZI; i++)
		{ // berechne die Summe
			
			if(KUPF_BINARY)
			{
				temp_vec[i] =  coupling[i] *  (yin[(2 + i * 2)  * DIM  + j] - yin[j]); 		// dp_q/dt = -w^2 * q+ gamma * x
			}else
			{
				temp_vec[i] = coupling[i] * yin[(2 + i * 2)  * DIM  + j]\
						 - pow(coupling[i],2.0)/pow(ommega[i],2.0) * yin[j];  // for Kahan summation : vec gamma^2/ommega^2
			}
		}
		Sum_Method(OSSZI, temp_vec, &sum);
		ans[j] = yin[DIM+ j]/mass;  // dx/dx = p/m
		ans[DIM+j] = sum  ;   	//p-deriv
	}
}

