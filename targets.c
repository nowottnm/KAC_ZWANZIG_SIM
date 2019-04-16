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

int Get_Lattice_Targets_Center_Circle( int *pos, double *targets, 	// leaves esc_cell radius free
						 double lattice_spacing,  int escape_cells)	// targets one in cell center
{
	int ret_value = 0; //leave starting cell empty, else particle starts in Target
	if 	(	(pos[0] == 0)&&
			(pos[1] == 0) 
		)

	{
		ret_value = 0; //catch mistake if zero cell was still supposed to have target
		return (ret_value);
	}
	if 	(pos[0] * pos[0] + pos[1] * pos[1] > escape_cells * escape_cells ) //If outside circle
	{
		ret_value = 1;    //single target center cell
		targets[0] = pos[0] * lattice_spacing; 
		targets[1] = pos[1] * lattice_spacing;
		return (ret_value);
	}
	return (ret_value);
}

int Get_Lattice_Targets_PORE_3( int *pos, double *targets, 	// leaves esc_cell radius free
						 double lattice_spacing,  int escape_cells)	// targets one in cell center
{	
	double corner_vec[2];
	int corner = 0;
	//corner one oben links
	for (double i_x = -1.0; i_x < 2.0; i_x+= 2.0)
	{
		for (double i_y = -1.0; i_y < 2.0; i_y+= 2.0)
		{ 
			int stride = 2 * corner * (NUMBER_PORE_OBSTACLES + 1);
			corner_vec[0] =  (pos[0] + i_x *0.5) * LATTICE_SPACING;
			corner_vec[1] =  (pos[1] + i_y *0.5) * LATTICE_SPACING;
			targets[0 + stride] = corner_vec[0];
			targets[1 + stride] = corner_vec[1];
			for (int i_pore = 0; i_pore < NUMBER_PORE_OBSTACLES; i_pore++)
			{
				targets[2 + stride + i_pore * 4] = corner_vec[0] - i_x * TARGET_LENGTH * (i_pore + 1);
				targets[3 + stride + i_pore * 4] = corner_vec[1];

				targets[4 + stride + i_pore * 4] = corner_vec[0];
				targets[5 + stride + i_pore * 4] = corner_vec[1] - i_y * TARGET_LENGTH * (i_pore + 1);
			}
			corner++;
		}
		corner++;
	}
	return (4 *(1 + 2* NUMBER_PORE_OBSTACLES));
}

int Get_Lattice_Targets_Center_Square( int *pos, double *targets,	// leaves esc_cell^2 square free
						 double lattice_spacing,  int escape_cells)	// targets one in cell center
{
	int ret_value = 0; //leave starting cell empty, else particle starts in Target
	if 	(	(pos[0] == 0)&&
			(pos[1] == 0) 
		)

	{
		ret_value = 0; //catch mistake if zero cell was still supposed to have target
		return (ret_value);
	}
	if 	( 	(abs(pos[0]) > escape_cells)||
			(abs(pos[1]) > escape_cells) 
		)
	{
		ret_value = 1;    //single target center cell
		targets[0] = pos[0] * lattice_spacing; 
		targets[1] = pos[1] * lattice_spacing;
		return (ret_value);
	}
	return (ret_value);
}

int Get_Lattice_Targets_Corner_Square( int *pos, double *targets,
						 double lattice_spacing,  int escape_cells)
// get target vectors for corner, middle or porous region. save in targest[4*DIM], return number targets
{
	int ret_value = 0;
	if(escape_cells == 0)		// catch case no free cells
	{
		ret_value = 4; // 4 targets per lattice cell
		int count = -1;
		for (int i = -1; i< 2; i+=2) //i,j in (-1,1)
		{
			for (int j = -1; j< 2; j+=2)
			{
				count++;
				targets[0 + count * DIM] = i * lattice_spacing/2.0 + pos[0] * lattice_spacing; 
				targets[1 + count * DIM] = j * lattice_spacing/2.0 + pos[1] * lattice_spacing;
			}
		}
		return (ret_value);
	}
	// check porous are
	if( (abs(pos[0]) > escape_cells)||
		(abs(pos[1]) > escape_cells) )
	{
		ret_value = 4; // 4 targets per lattice cell
		int count = -1;
		for (int i = -1; i< 2; i+=2) //i,j in (-1,1)
		{
			for (int j = -1; j< 2; j+=2)
			{
				count++;
				targets[0 + count * DIM] = i * lattice_spacing/2.0 + pos[0] * lattice_spacing; 
				targets[1 + count * DIM] = j * lattice_spacing/2.0 + pos[1] * lattice_spacing;
			}
		}
		return (ret_value);
	}
	//Check corner 
	for (int i = -1; i< 2; i+=2) //i,j in (-1,1)
	{
		for (int j = -1; j< 2; j+=2)
		{
			if( (pos[0] == i * escape_cells)&&
				(pos[1] == j * escape_cells))
			{	
				ret_value = 3; // corner and two sides
				targets[0] = i * lattice_spacing/2.0 + pos[0] * lattice_spacing; 		//corner
				targets[1] = j * lattice_spacing/2.0 + pos[1] * lattice_spacing;

				targets[2] = i * (-1) * lattice_spacing/2.0 + pos[0] * lattice_spacing;	//first side
				targets[3] = j * lattice_spacing/2.0 + pos[1] * lattice_spacing;
				targets[4] = i * lattice_spacing/2.0 + pos[0] * lattice_spacing;	//first side
				targets[5] = j * (-1) * lattice_spacing/2.0 + pos[1] * lattice_spacing;
				return (ret_value);
			}	// corner
		}
	}
	// check sides after sure you are not in corner
	for( int d = 0; d < DIM; d++)// loop over dimension
	{
		int e = 1-d;
		for (int i = -1; i< 2; i+=2)
		{
			if(pos[d] == i * escape_cells)
			{
				ret_value = 2; //two sides
				targets[d] 			= i * lattice_spacing/2.0 + pos[d] * lattice_spacing;
				targets[e] 			= (-1) * lattice_spacing/2.0 + pos[e] * lattice_spacing;
				targets[d +DIM] 	= i * lattice_spacing/2.0 + pos[d] * lattice_spacing;
				targets[e +DIM] 	= (1) * lattice_spacing/2.0 + pos[e] * lattice_spacing;
				return (ret_value);
			}
		}
	}
	return (ret_value);
}