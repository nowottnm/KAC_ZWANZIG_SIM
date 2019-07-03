README for the KAC ZWANZIG sim, used in the Master thesis of M.Nowottnick. This program approximates generalized Langevin Equations of the type
d/dtv = K * v + random Force
K*v = int(K(t') v(t-t'))dt'
with an explicit deterministic bath providing both the friciton Kernel K and the random Force (on average) via interaction with a distinguished particle. We focus on the Kernel of power law type (leading to subdiffusion), so the convolution K*v becomes a fractal derivative and the whole equation is called the fractal Langevin equation.

Read the uploaded introductory chapters of the thesis for details of the algorithm as well as some subtle changes made on Kuperman's implementation of the bath of Kac-Zwanzig type for numerical implementation.

For this to compile, you need the gnu scientific library and openmp installed

please use the makefile:    >> make clean 
                            >> make
                            >> ./Badsym   

ideally, you should only change data in global.h. A short list of the parameters you can change are:

OSSZI: bath size, please use multiples of ten
SIMULATIONEN: number of runs for averaging
KUPF_BINARY: set to 1 for fully mobile bath, 0 for Zwanzig's original   
model, when in doubt leave 1
DIM : change dimension of the problem, currently only 1 or 2 supported
TIME_STEPS: stepsize at which to record data
TIME_END: number of such steps. NEEDS!! format (I + 9)
THREADNR: number threads
TEMP: Temperature in Kelvin
GAMMA: constant before bath coupling constants, K_0 in masterh thesis
ALPHA: desired scaling of the MSD ~t^ALPHA
mass: mass distinguished particle
TIME_ARRIVAL : geometries are specified in terms of time for diffusion over lattice cell, i.e. if one cell should be diffused over in 1000s set TIMEA_ARRIVAL to 1000.0
MAXSTEPSIZEMULT: integrator stepsize is ballstic time scale time this, at least .1 recommended
LABEL: "ZWANZIG" for P_0=0, "ZWANZIG_THERMALIZED" for equilibrated start
STEPPER METHOD: sets target. Pick from 'VVerlet_Step_deriv_Box' (1-D box), 'VVerlet_Step_deriv_Kupf' for free diffusion or 
'VVerlet_Step_hardsphere_reflect_Kupf' for square lattice of spheres
MAX_NR_PLOTS: Number of example trajectorier in you trajec folder of the simulation

everything that is not mentioned is best left alone, it might be legacy code needed for some obscure subroutines.

Now, at the end of the run you get a list of files into your results folder

These files are:

ekin.dat : paris of time and kinetic energy of distinguished particle
ekinbath.dat : same for total bath kinetic energy
ETOT: times and total bath energy
global.h: copy of used global.h
latlngth.dat: lattice size used
ommega.dat: list of frequencies of bath
P_X.dat : times and total X-Momentum
squares_rohdaten.dat: times and MSD

the folder /trajec also contains a number of sample trajectories (number specified in global.h)

in 1-D in the format t-x, else in x-y
the plot folder is just for convenience for eventual pictures.


Some rough notes on how to modify this:

If you just want to specify another geometry WITHIN square lattices, write your own method in target.c and set TARGET_METHOD to your version in the globalh. Note that the particle starts in the cell center, so leave it free or do a case separation for the first lattice cell.

If you want a completely new geometry, you will have to write a new Steppmethod in steppmethods.c.

Note that in constants.c new stepp methods have to be made aware to the routine setting the bath size.

If you want to add soft potentials, be aware that you have to give the Verlet method the correct first derivative on top of implementing a new step method.


