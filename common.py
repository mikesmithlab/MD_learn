"""Common - Contains code and declarations of global variables / functions"""

import numpy as np
from sphere import Sphere
import json


<<<<<<< HEAD

=======
lx, ly #System size
x0, y0 # coord lower left corner
no_of_particles
nstep
nprint
nenergy
time=0 # elapsed time
timestep
>>>>>>> a06d25bfafd2b26e1302978aa0fdc4da0d46815b

output_folder = '~/Documents/Programming/MD_learn/Sim_Data/'
fphase ='phase.dat'# output file - phase space
fenergy='energy.dat' # output file - system energy

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
<<<<<<< HEAD
    
=======
    """Dumps all particle coords and params to a file"""

>>>>>>> a06d25bfafd2b26e1302978aa0fdc4da0d46815b

def init_system(init_filename=None):
    with open(init_filename,"r") as f:
        params = json.load(f)
    return params

    particles = []

def dump_particle(f, index, x, y, 0, vx, vy, 0, radius, mass, ptype):
    "Open f handle using with open('test3.dat','a') as f:"
    np.savetxt(f, [[index, x, y, 0, vx, vy, 0, radius, mass, ptype]])






def total_kinetic_energy():
    """Calculates total kinetic energy of all free particles"""
    sum=0
    for particle in particles:
        if particle.ptype == 0:
            sum+=particle[i].kinetic_energy
    return sum


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