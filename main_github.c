#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include "global.h"
#include "utility.h"
#include "gnuplot.h"
#include "constants.h"
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
#include "steppmethods.h"

#define NUMBER_SAMPLES 1
#define NUMBER_VOLUMES 1	
#define SHIFT_LATE 200
#define BINS 100			
#define PDF_NR 12



// regular main but with multiple sampling runs over multiple volume fractions phi.


int
main (void){
	FILE *source, *target;
  	source = fopen("global.h", "r"); // get source global.h at program start so as to not copy edits during run

  	if( source == NULL )
  	{
  		printf("cannot open global.h\n");
  		exit(EXIT_FAILURE);
  	}
  	SETNUMTHREADS
  	int i,j,l,k,m ,n;
	Constants();									// draw coupling constants and other derived constants
	double *t = calloc(LENGTH_T, sizeof(double));
	time_t tme = time(NULL);
	char ordnerLABEL[300];
	char ueberordnerLABEL[300];
	char ordnerLABEL_original[300];
	strftime(ordnerLABEL, sizeof ordnerLABEL, "%A %c", localtime(&tme));
  	double squares_all_calcs[LENGTH_T][NUMBER_SAMPLES];			// for printing all squares in shared table
  	i = 0;
  	printf("\n%s\n", ordnerLABEL);
  	deleteSpaces(ordnerLABEL,ordnerLABEL);
  	strcpy (ordnerLABEL_original, ordnerLABEL);    	// keep original, change ordnerLABEL per run
  	double alpha_global = ALPHA;
  	double temp_global = TEMP;						// macros from global make trouble when used as print arguments, hence copy to memory
  	sprintf(ueberordnerLABEL,"N=%d_SIM=%d_T=%1.1E_ALPH=%1.1f", OSSZI, SIMULATIONEN, temp_global, alpha_global);
  	char tempLABEL[300];
  	strcpy (tempLABEL, ordnerLABEL);
  	strcat(tempLABEL, LABELPREFIX);
  	strcat(tempLABEL, ueberordnerLABEL);
  	strcpy (ueberordnerLABEL, tempLABEL);
  	struct stat st2 = {0};

  	double prob_density[BINS][PDF_NR];
  	for (i = 0; i < BINS ; i++)
  	{
  		for (j = 0; j < PDF_NR; j++)
  		{
  			prob_density[i][j] = 0.0;
  		}
  	}
  	int check_times[PDF_NR];
  	for (j = 0; j < PDF_NR; j++)
  	{
  		check_times[j] = (int) ( 10.0 * (j+1)) ;
  	}
  	check_times[PDF_NR-1] = LENGTH_T - 1;
  	for (j = 0; j < PDF_NR; j++)
  	{
  		printf("check times nr %d is %d \n", j, check_times[j]);
  	}
  	printf("\nprint results into ");
  	printf(ueberordnerLABEL);

  	if (stat(ueberordnerLABEL, &st2) == -1)
	{												//Teste ob Ordner existiert, erstelle Ordner
		mkdir(ueberordnerLABEL, 0700);	
	}
	if (chdir(ueberordnerLABEL))						// change into directory of simulation
	{
		printf("Error changing directory");
		return 1;
	}

  	// reserviere speicher fuer eine Vollständige trajektorie ------------------
	m = LENGTH_T; n = ORDER;
	double **z = (double **) malloc(m * sizeof(double *));
	z[0] = (double *) malloc(m* n * sizeof(double));
	for (i=1; i<m; i++) z[i] = z[0] + n *i; 


		double *squares = calloc(LENGTH_T, sizeof(double));
	double *short_correlation_x = calloc(LENGTH_T, sizeof(double));	//short term (1 timestep) correlations
	double *short_correlation_y = calloc(LENGTH_T, sizeof(double));
	double *total_momentum_x = calloc(LENGTH_T, sizeof(double));
	double *total_momentum_y = calloc(LENGTH_T, sizeof(double));
	double *long_correlation_x = calloc(LENGTH_T, sizeof(double)); //longterm correlation from first timestep
	double *long_correlation_y = calloc(LENGTH_T, sizeof(double));
	double *squares_increase = calloc(LENGTH_T, sizeof(double));	// Increase a from assumed form x^2 0 = a t
	double real_volumes[10];

	double *bathp = calloc(LENGTH_T, sizeof(double));
	double *bathq = calloc(LENGTH_T, sizeof(double));

	double *px_correlation = calloc(LENGTH_T, sizeof(double));	// save correlation of impulse_x to px(t=0)
	double *py_correlation = calloc(LENGTH_T, sizeof(double));	
	double *px_correlation_late =  calloc(LENGTH_T, sizeof(double));	// corelation after time SHIFT_LATE
	double *px_correlation_late2 = calloc(LENGTH_T, sizeof(double));	// corelation after time SHIFT_LATE
	

	int n_zplots = SIMULATIONEN;
	if (SIMULATIONEN > MAX_NR_PLOTS ) n_zplots = MAX_NR_PLOTS;    // speichere maximal 30*DIM trajectorien
  	// reserviere speicher fuer samplepfade massives Teilchen zum plotten ------------------	
	m = LENGTH_T; n = DIM*n_zplots;
	double **zplots = (double **) malloc(m * sizeof(double *));
	zplots[0] = (double *) malloc(m* n * sizeof(double));
	for (i=1; i<m; i++) zplots[i] = zplots[0] + n *i; 


	// reserviere speicher fuer samplepfade 1. bad  Teilchen zum plotten ------------------	
		m = LENGTH_T; n = DIM*n_zplots;
	double **qplots = (double **) malloc(m * sizeof(double *));
	qplots[0] = (double *) malloc(m* n * sizeof(double));
	for (i=1; i<m; i++) qplots[i] = qplots[0] + n *i; 


	// reserviere speicher fuer samplepfade letztes bad  Teilchen zum plotten ------------------	
		m = LENGTH_T; n = DIM*n_zplots;
	double **qlplots = (double **) malloc(m * sizeof(double *));
	qlplots[0] = (double *) malloc(m* n * sizeof(double));
	for (i=1; i<m; i++) qlplots[i] = qlplots[0] + n *i; 

	// reserviere speicher fuer samplepfade zusätzlicher Kac_Zwanzig Kraft Term ------------------	
		m = LENGTH_T; n = DIM*n_zplots;
	double **Zwanzig_Force = (double **) malloc(m * sizeof(double *));
	Zwanzig_Force[0] = (double *) malloc(m* n * sizeof(double));
	for (i=1; i<m; i++) Zwanzig_Force[i] = Zwanzig_Force[0] + n *i; 



	// reserviere speicher fuer probability Density ------------------	
		m = 80; n = 80;
	double **P_density = (double **) malloc(m * sizeof(double *));
	P_density[0] = (double *) malloc(m* n * sizeof(double));
	for (i=1; i<m; i++) P_density[i] = P_density[0] + n *i; 
  	// reserviere speicher fuer Gitterkoordinaten ------------------

		double DENSITY_POINTS = 0.0;
	
	for(m=0; m<80; m++)
	{
		for(n=0; n<80; n++)
		{
			P_density[m][n] = 0.0;
		}
	}	

	double *EKIN = calloc(LENGTH_T, sizeof(double));
	double *EBAD = calloc(LENGTH_T, sizeof(double));
	double *ETOT = calloc(LENGTH_T, sizeof(double));
	double *PTOT = calloc(LENGTH_T, sizeof(double));
	double *PTOTY = calloc(LENGTH_T, sizeof(double));
	double *LTOT = calloc(LENGTH_T, sizeof(double));


	// GSL random number Setup for Taus Generator
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_ranlxd2;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r, time(NULL)); // Seed with time
	// RNG setup End // 
	// GSL random number Setup for Taus Generator for global rand
	const gsl_rng_type * T2;
	T2 = gsl_rng_ranlxd2;
	RAND_GLOBAL = gsl_rng_alloc (T2);
	gsl_rng_set(RAND_GLOBAL, time(NULL)); // Seed with time
	// RNG setup End // 

	clock_t start, end; 

	printf("evenly spread t \n");	//set up time vec
  	for (i = 0; i < 10; i++) // first ten values in small steps
  	{
  		t[i] =  i * TIME_STEPS/10.0;
  	}
  	for (i = 10; i < LENGTH_T; i++)
  	{
  		t[i] = (i-9) * TIME_STEPS;
  	}
  	double squarespace = 0.0;

  	DIFF_COEFF = KBOLTZ * TEMP/GAMMA * sin(M_PI * ALPHA ) / (M_PI * ALPHA);
  	if(ALPHA > 0.95)
  	{
  		DIFF_COEFF = KBOLTZ * TEMP/GAMMA;
  	}
  	LATTICE_SPACING =  sqrt(2*DIM*DIFF_COEFF * pow( TIME_ARRIVAL, ALPHA) );
  	printf("LATTICE_SPACING = %3.3E\n", LATTICE_SPACING);
  	double latlength = LATTICE_SPACING;

		//------------------------------EIGENTLICHE SIMULATION ------------------------------------------------------------------
  	int sim_count = 0;
  	start = clock();
  	time_first_contact = 0.0;
  	int error_count = 0;
		// 
  	for (i = 0; i < SIMULATIONEN; i++)
  	{

		vec_zero_i(lattice_position,DIM);		// set stating cell to zero !
		TARGET_CONTACT = 0;
		FLAG_DUMP_VALUE = 0;					// set to 1 if error occures
		// -------------------------BAD
		Bath_Setup(y, LABEL, r);				// draw intial conditions
		//---------------------Bad aufgesetzt	
		VVerlet_parallel(ORDER, y, z, t, Poti_Handle);
		end = clock();
		//------------------------------auswertung Integrationsergebnisse
		if ( !(FLAG_DUMP_VALUE))			// only count result without err
		{
		  	double *kahan_c = calloc(LENGTH_T, sizeof(double));                //
		  	double *kahan_t = calloc(LENGTH_T, sizeof(double));
		  	double *kahan_y = calloc(LENGTH_T, sizeof(double));
		  	#pragma omp parallel for
		  	for(l = 0; l < LENGTH_T; l++)
		  	{
		  		if (l>10)
		  		{
		  			long_correlation_x[l] += z[1][0] *z[l][0]/((double) SIMULATIONEN);
		  			long_correlation_y[l] += z[1][1] *z[l][1]/((double) SIMULATIONEN);
		  		}
		  		kahan_y[l] = z[10][DIM] * z[l][DIM] / ((double) SIMULATIONEN) / (mass * KBOLTZ * TEMP)\
		  					- kahan_c[l];
		  		kahan_t[l] = px_correlation[l] + kahan_y[l];
		  		kahan_c[l] = (kahan_y[l] - px_correlation[l]) - kahan_y[l];
		  		px_correlation[l] = kahan_t[l];
		  		py_correlation[l] += z[10][DIM + 1] * z[l][DIM + 1] / ((double) SIMULATIONEN) / (mass * KBOLTZ * TEMP);
		  		px_correlation_late[l] += z[SHIFT_LATE][DIM] * z[l][DIM] / ((double) SIMULATIONEN) / (mass * KBOLTZ * TEMP);
		  		px_correlation_late2[l] += z[SHIFT_LATE * 2][DIM] * z[l][DIM] / ((double) SIMULATIONEN) / (mass * KBOLTZ * TEMP);
		  		for(j=0; j< DIM; j++)
		  		{
		  			if (i < n_zplots)
		  			{
		  				zplots[l][DIM * i + j] = z[l][j];
		  				qplots[l][DIM * i + j] = z[l][2*DIM + j];
		  				double Force_TEMP = 0.0;
		  				for(int os = 0; os < OSSZI; os++)
		  					{
		  						Force_TEMP += (pow(coupling[os],2.0)/pow(ommega[os],2.0) - coupling[os]) * z[l][(2 + 2 * os) * DIM + j];
		  					}
		  					Zwanzig_Force[l][DIM * i + j]=Force_TEMP;
		  					int i_last = OSSZI -1;
		  					qlplots[l][DIM * i + j] = z[l][2*OSSZI*DIM + j];
		  			} 
		  			squares[l] = squares[l] + pow(z[l][j],2.0)/((double) SIMULATIONEN);
		  			EKIN[l] = EKIN[l] + pow(z[l][DIM+j],2.0)/(mass * SIMULATIONEN * KBOLTZ *TEMP * DIM);
		  			for(k=0;k<OSSZI;k++)
		  			{
		  				bathq[l] = bathq[l] + pow(z[l][(2 + 2*k) *DIM +j],2.0)/( (double) (OSSZI*SIMULATIONEN)\
		  						*KBOLTZ*TEMP) * pow(ommega[k],2.0)/2.0;// mittlere potentielle E pro Badteilchen
		  																//reskaliert durch kT
		  				bathp[l] = bathp[l] + pow(z[l][(3 + 2*k) *DIM +j],2.0)/( (double) (OSSZI*SIMULATIONEN)\
		  						*KBOLTZ*TEMP) /2.0/massq[k];						// mittlere kinetische E pro Badteilchen
		  																//reskaliert durch kT
		  				EBAD[l] = EBAD[l] + pow(z[l][(3 + 2*k) *DIM +j], 2.0)/2.0\
		  						+ 0.5  * pow(ommega[k],2.0) \
		  						* pow(z[l][(2 + 2*k) *DIM +j] - coupling[k]/(pow(ommega[k],2.0) ) * z[l][j] , 2.0) /((double) SIMULATIONEN);
		  			}
		  		}
			  	for(k=0;k<OSSZI;k++)	// add p and angular momentum L for bath
			  	{ 
			  		PTOT[l] 	+= z[l][(3 + 2*k) *DIM + 0];
			  		PTOTY[l]	+= z[l][(3 + 2*k) *DIM + 1];
			  		LTOT[l] += 	(z[l][(2 + 2*k) *DIM + 0] * z[l][(3 + 2*k) *DIM + 1]\
			  					- z[l][(2 + 2*k) *DIM + 1] * z[l][(3 + 2*k) *DIM + 0]) /((double) SIMULATIONEN);
			  	}
			  	PTOT[l] += z[l][2];
			  	PTOT[l] += z[l][3];
			  	LTOT[l] += (z[l][0] * z[l][3] - z[l][1] * z[l][2]) / ((double) SIMULATIONEN);
			}
			if(DIM == 1)
			{
			  	// setup pdf checks
			  	double bins = BINS * 1.0;
			  	double dx = LATTICE_SPACING/bins;
			  	double sims = SIMULATIONEN;
			  	for (j = 0; j < PDF_NR; j++)
			  	{
			  		int t_check = check_times[j];
			  		for (int i_bin = 0; i_bin < BINS; i_bin++)
			  		{
			  			double lower = -LATTICE_SPACING / 2.0 + i_bin * dx;
			  			double upper =  -LATTICE_SPACING / 2.0 + (i_bin + 1) * dx;
			  			if 	( 	( z[t_check][0] <= upper) &&
			  					( z[t_check][0] > lower) 
			  				)
			  			{
			  				prob_density[i_bin][j] += 1.0/sims;
			  			}
			  		}
			  	}
			}
			sim_count += 1;
			//printf("%d  und t %4.2f \n", i ,((double) (end - start)));
			printf("\r%d von %d mit t/count = %4.2f s und average t_rest =%4.2f h ", sim_count , SIMULATIONEN, (double) (end - start) /  sim_count/ THREADS / CLOCKS_PER_SEC ,\
			  		((double) (end - start) / sim_count/ THREADS/ CLOCKS_PER_SEC * (SIMULATIONEN - sim_count)/3600));
			printf(" %d hits on Lat  and avrgtcntct = %4.2f",TARGET_CONTACT, time_first_contact/(i+1) );
			fflush(stdout);
			for(l = 0; l < LENGTH_T; l++)
			{
			  	for (int oss_i=0; oss_i < OSSZI; oss_i++)
			  	{
			  		total_momentum_x[l] += z[l][(3 + 2*oss_i) *DIM + 0];
			  		total_momentum_y[l] += z[l][(3 + 2*oss_i) *DIM + 1];
			  	}
			  	total_momentum_x[l] += z[l][(0 + 2*0) *DIM + 0];
			  	total_momentum_y[l] += z[l][(0 + 2*0) *DIM + 1];
			} 
			// -------------- end evaluation if----------------------
			free(kahan_y); free(kahan_t); free(kahan_c);
		}else
		{
			error_count++;
			printf("\n err occured at calc %d, calculation dumped, %d totatl errors\n",  i,error_count );
			i -= 1;		// do one more calculation
		}

	}
	tme = time(NULL);
	char end_label[90];
	strftime( end_label, sizeof end_label, "%A %c", localtime(&tme));
	printf("\n%s\n", end_label);
	for(l = 0; l < LENGTH_T; l++)
	{
		ETOT[l] = EKIN[l] + EBAD[l];
			  	PTOT[l] = sqrt(pow(PTOT[l]/((double) SIMULATIONEN),2.0) + pow(PTOTY[l]/((double) SIMULATIONEN),2.0)); 
	}
	for(l = 1; l < LENGTH_T; l++)
	{
		squares_increase[l] = (squares[l]- squares[l-1])/(t[l] - t[l-1]);
	}
			//------------------------------ENDE SIMULATION, speichere daten ------------------------------------------------------------------
	for(m=0; m<80; m++)
	{
		for(n=0; n<80; n++)
		{
			P_density[m][n] *= 1.0/DENSITY_POINTS;
		}
	}
	FILE *fp;

	struct stat st = {0}; 			

	if (stat(ordnerLABEL, &st) == -1){		//Teste ob Ordner existiert, erstelle Ordner
		mkdir(ordnerLABEL, 0700);	
	}

	fp = fopen ("shellscript.sh", "w");		//create Shellscript to start Gnuplot from c main()
	fprintf(fp, "cd %s\n", ordnerLABEL);
	fprintf(fp, "gnuplot gnuplot.txt\n");
	fprintf(fp, "cd ..");
	fclose (fp);


	 // copy global.h 	/


   	// copy to same name into directory, clould be anything different

   	if (chdir(ordnerLABEL))			// change into directory of simulation
   	{
   		printf("Error changing directory");
   		return 1;
   	}

   	mkdir("plots", 0700);
   	mkdir("trajec", 0700);

   	target = fopen("global.h", "w");

   	if( target == NULL )
   	{
   		fclose(source);
   		printf("Press any key to exit...\n");
   		exit(EXIT_FAILURE);
   	}
   	char ch;
   	while( ( ch = fgetc(source) ) != EOF )
   	{
   		fputc(ch, target);
   	}

   	printf("File copied successfully.\n");


   	fclose(target);

   	fp = fopen ("latlength.dat", "w");
   	fprintf(fp, "%lf  \n", latlength);
   	fclose (fp);


   	fp = fopen ("squares_rohdaten.dat", "w");
   	for(l = 0; l < LENGTH_T; l++){
   		fprintf(fp, "%1.3E  %1.3E\n", t[l], squares[l]);
   	}
   	fclose (fp);

   	fp = fopen ("ekin.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%lf  %lf\n", t[l]/TIME_ARRIVAL, EKIN[l]);
   	}
   	fclose (fp);

   	fp = fopen ("ekinbath.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%lf  %lf\n", t[l]/TIME_ARRIVAL, bathp[l]);
   	}
   	fclose (fp);

   	fp = fopen ("ommega.dat", "w");
   	for(l = 0; l < OSSZI; l++){
   		fprintf(fp, "%lf\n", ommega[l]);
   	}
   	fclose (fp);

   	fp = fopen ("PTOT.dat", "w");
   	for(l = 0; l < LENGTH_T; l++){
   		fprintf(fp, "%lf  %lf\n", t[l]/TIME_ARRIVAL, PTOT[l]);
   	}
   	fclose (fp);

   	fp = fopen ("P_X.dat", "w");
   	for(l = 0; l < LENGTH_T; l++){
   		fprintf(fp, "%lf  %1.3e\n", t[l]/TIME_ARRIVAL, total_momentum_x[l]);
   	}
   	fclose (fp);

   	fp = fopen ("P_Y.dat", "w");
   	for(l = 0; l < LENGTH_T; l++){
   		fprintf(fp, "%lf  %1.3e\n", t[l]/TIME_ARRIVAL, total_momentum_y[l]);
   	}
   	fclose (fp);

   	fp = fopen ("xshortcorellation.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%lf  %lf\n", t[l]/TIME_ARRIVAL, short_correlation_x[l]/latlength/latlength);
   	}
   	fclose (fp);

   	fp = fopen ("yshortcorellation.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%lf  %lf\n", t[l]/TIME_ARRIVAL, short_correlation_y[l]/latlength/latlength);
   	}
   	fclose (fp);

   	fp = fopen ("xlongcorellation.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%lf  %lf\n", t[l]/TIME_ARRIVAL, long_correlation_x[l]/latlength/latlength);
   	}
   	fclose (fp);

   	fp = fopen ("ylongcorellation.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%lf  %lf\n", t[l]/TIME_ARRIVAL, long_correlation_y[l]/latlength/latlength);
   	}
   	fclose (fp);

   	fp = fopen ("PX_corr.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%1.3E  %1.3E\n", t[l]/TIME_ARRIVAL, px_correlation[l]);
   	}
   	fclose (fp);

   	fp = fopen ("PY_corr.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%1.3E  %1.3E\n", t[l]/TIME_ARRIVAL, py_correlation[l]);
   	}
   	fclose (fp);

   	fp = fopen ("PX_corr_late.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%1.3E  %1.3E\n", t[l]/TIME_ARRIVAL, px_correlation_late[l]);
   	}
   	fclose (fp);

   	fp = fopen ("PX_corr_late2.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%1.3E  %1.3E\n", t[l]/TIME_ARRIVAL, px_correlation_late2[l]);
   	}
   	fclose (fp);

   	fp = fopen ("PX_corr_abs.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%1.3E  %1.3E\n", t[l]/TIME_ARRIVAL, fabs(px_correlation[l]));
   	}
   	fclose (fp);

   	fp = fopen ("PX_corr_late_abs.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%1.3E  %1.3E\n", t[l]/TIME_ARRIVAL, fabs(px_correlation_late[l]));
   	}
   	fclose (fp);

   	fp = fopen ("PX_corr_late2_abs.dat", "w");
   	for(l = 1; l < LENGTH_T; l++){
   		fprintf(fp, "%1.3E  %1.3E\n", t[l]/TIME_ARRIVAL, fabs(px_correlation_late2[l]));
   	}
   	fclose (fp);


   	if (DIM == 1)
   	{
   		fp = fopen ("PDF_1D.dat", "w");
   		for (int pos = 0; pos < BINS; pos++)
   		{
   			fprintf(fp, "%1.3E ", -LATTICE_SPACING/2.0 + pos*LATTICE_SPACING/((double) BINS));
   			for(j = 0; j < PDF_NR; j++)
   			{
   				fprintf(fp, "%1.3E ", prob_density[pos][j]);
   			}
   			fprintf(fp,"\n");
   		}
   		fclose (fp);
   	}


	// finde größtes s^2/t^alpha zum endzeitpunkt
   	int t_start = 10;
   	double temp_ende = 0.0;
   	for(int i_alpha = 1; i_alpha < 12; i_alpha ++){
   		double alpha_temp = 0.45 + 0.05 * i_alpha;
   		l = LENGTH_T-1;
   		if (temp_ende < (squares[l] / (pow(t[l], alpha_temp))) )
   		{
   			temp_ende = squares[l] / (pow(t[l], alpha_temp));
   		}
   	}
   	if (DIM>1) 
   	{

   		char zplotsLABEL[30];
   		for (i = 0; i < n_zplots; i++)
   		{
   			sprintf(zplotsLABEL, "trajec/zplots%d.dat",i);
   			fp = fopen (zplotsLABEL, "w");
   			for(l = 0; l < LENGTH_T; l++){
   				for (j = 0; j < DIM ; j++){
   					fprintf(fp, "%lf  ", ((zplots[l][j +i*DIM])/latlength) );
   				}	
   				fprintf(fp, "\n");
   			}
   			fclose (fp);
   		}
   		for (i = 0; i < n_zplots; i++)
   		{
   			sprintf(zplotsLABEL, "trajec/qplots%d.dat",i);
   			fp = fopen (zplotsLABEL, "w");
   			for(l = 0; l < LENGTH_T; l++){
   				for (j = 0; j < DIM ; j++){
   					fprintf(fp, "%lf  ", ((qplots[l][j +i*DIM])) );
   				}	
   				fprintf(fp, "\n");
   			}
   			fclose (fp);
   		}
   		for (i = 0; i < n_zplots; i++)
   		{
   			sprintf(zplotsLABEL, "trajec/Zwanzig_Force%d.dat",i);
   			fp = fopen (zplotsLABEL, "w");
   			for(l = 0; l < LENGTH_T; l++){
   				for (j = 0; j < DIM ; j++){
   					fprintf(fp, "%lf  ", ((Zwanzig_Force[l][j +i*DIM])) );
   				}	
   				fprintf(fp, "\n");
   			}
   			fclose (fp);
   		}
   		for (i = 0; i < n_zplots; i++)
   		{
   			sprintf(zplotsLABEL, "trajec/qlplots%d.dat",i);
   			fp = fopen (zplotsLABEL, "w");
   			for(l = 0; l < LENGTH_T; l++){
   				for (j = 0; j < DIM ; j++){
   					fprintf(fp, "%lf  ", ((qlplots[l][j +i*DIM])) );
   				}	
   				fprintf(fp, "\n");
   			}
   			fclose (fp);
   		}
   	}
   	if (DIM==1) 
   	{

   		char zplotsLABEL[30];
   		for (i = 0; i < n_zplots; i++)
   		{
   			sprintf(zplotsLABEL, "trajec/zplots%d.dat",i);
   			fp = fopen (zplotsLABEL, "w");
   			for(l = 0; l < LENGTH_T; l++){
   				fprintf(fp, " %lf  %lf \n" ,t[l]/TIME_ARRIVAL, zplots[l][i*DIM]/latlength);
   			}
   			fclose (fp);
   		}
   	}



   	fp = fopen ("ETOT.dat", "w");
   	for(l = 0; l < LENGTH_T; l++){
   		fprintf(fp, "%lf  %lf\n", t[l], ETOT[l]);
   	}
   	fclose (fp);

   	fp = fopen ("DENSITY.dat", "w");
   	for(m=0; m<80; m++)
   	{
   		for(n=0; n<80; n++)
   		{
   			fprintf(fp, "%1.2e ", P_density[m][n]);
   		}
   		fprintf(fp, "\n");
   	}
   	fclose (fp);


	gsl_rng_free (r); gsl_rng_free (RAND_GLOBAL);
	free(zplots[0]); free(zplots);
	free(qplots[0]); free(qplots);
	free(Zwanzig_Force[0]); free(Zwanzig_Force);
	free(qlplots[0]); free(qlplots);
	free (P_density[0]); free (P_density);
	free(t);
	free(bathp); free(bathq);
	free(squares_increase);
	free(short_correlation_x);
	free(total_momentum_x);
	free(total_momentum_y);
	free(short_correlation_y);
	free(long_correlation_x);
	free(px_correlation); free(py_correlation);
	free(px_correlation_late);
	free(px_correlation_late2);
	free(long_correlation_y);
	free(ETOT); free(EKIN); free(EBAD);
	free(PTOT); free(PTOTY);free(LTOT);
	free(squares);
	free(z[0]); free(z);
}
