README for the KAC ZWANZIG sim

For this to compile, you need the gnu scientific library and openmp installed

please use the makefile:    >> make clean 
                            >> make
                            >> ./Badsym   

ideally, you should only change data in global.h. A short list of the parameters you can change are:

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

the folder trajec also contains a number of sample trajectories (number specified in global.h)

in 1-D in the format t-x, else in x-y
the plot folder is just for convenience for eventual pictures.


