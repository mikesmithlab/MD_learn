"""Common - Contains code and declarations of global variables / functions"""

import numpy as np
from sphere import Sphere


lx, ly #System size
x0, y0 # coord lower left corner
no_of_particles 
nstep
nprint
nenergy
time=0 # elapsed time
timestep

output_folder = '~/Documents/Programming/MD_learn/Sim_Data/'
fphase ='phase.dat'# output file - phase space
fenergy='energy.dat' # output file - system energy
G # vector for gravity

particle =[]# vector containing all objects of type Sphere

def step():
    pass

def make_forces():
    for i in range(particle.size()):
        for k in range(particle.size()):
            if (particle[i].ptype()==0) | (particle[k].ptype()==0):
                force(particle[i], particle[k],lx,ly)


def integrate():
    #Mobile particles have particle.ptype()==0
    for i in range(particle.size()):
        if particle[i].ptype()==0:
            particle[i].set_force_to_zero()
            particle[i].predict(timestep)
        else:
            particle[i].boundary_conditions(i, timestep, time)
    make_forces()
    for i in range(particle.size()):
        if particle[i].ptype()==0:
            particle[i].correct(timestep)
    for i in range(particle.size()):
        particle[i].periodic_bc(x_0, y_0, lx, ly)
    time+=timestep


def init_algorithm():
    pass

def phase_plot():
    pass

def init_system():
    pass

def total_kinetic_energy():
    pass


def main(init_file):
    init_system(init_file)
    init_algorithm()
    phase_plot()
    for i in range(nstep):
        step()
        if (i+1)%nprint == 0:
            phase_plot(fphase)
        if (i+1)%nenergy == 0:
            total_kinetic_energy(fenergy)





if __name__ == '__main__':
    main(init_file)