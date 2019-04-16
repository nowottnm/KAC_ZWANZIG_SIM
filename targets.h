
// sammlung sunktionen zur Erzeugung der Targets

#ifndef HEADER_FILE_TARGET
#define HEADER_FILE_TARGET

int Get_Lattice_Targets_Center_Circle( int *pos, double *targets, 	// leaves esc_cell radius free, centers targets in lattice cell
						 double lattice_spacing,  int escape_cells);

int Get_Lattice_Targets_PORE_3( int *pos, double *targets, 	// leaves esc_cell radius free
						 double lattice_spacing,  int escape_cells);

int Get_Lattice_Targets_Center_Square( int *pos, double *targets,	// leaves esc_cell^2 square free, centers targets
						 double lattice_spacing,  int escape_cells);

int Get_Lattice_Targets_Corner_Square( int *pos, double *targets,	// corner targets, leave squae size escape cells empty
						 double lattice_spacing,  int escape_cells);

#endif