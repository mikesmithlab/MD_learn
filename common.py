"""Common - Contains code and declarations of global variables / functions"""

import numpy as np
from materials import material_a



output_folder = '~/Documents/Programming/MD_learn/Sim_Data/'
fphase ='phase.dat'# output file - phase space
fenergy='energy.dat' # output file - system energy

particle =[]# vector containing all objects of type Sphere

def step(dt,time,particles):
    integrate(dt,time,particles)
    return particles

def make_forces():
    for i in range(particle.size()):
        for k in range(particle.size()):
            if (particle[i].ptype()==0) | (particle[k].ptype()==0):
                force(particle[i], particle[k],lx,ly)


def integrate(dt,time,particles):
    #Mobile particles have particle.ptype()==0
    for particle in particles:
        if particle.ptype==0:
            particle.set_force_to_zero()
            particle.predict(dt)
        else:
            particle.boundary_conditions(i, dt, time)
    make_forces()
    for particle in particles:
        if particle.ptype()==0:
            particle.correct(dt)
    for particle in particles:
        particle.periodic_bc(x_0, y_0, lx, ly)
    time+=dt
    return particles


def init_algorithm():
    pass

def phase_plot():
    pass


def init_system(init_filename=None):
    """Generate initial state of simulation

    Reads coords from .dat file created in setup.py and stored in sim_input_data

    returns a list of particle instances
    """
    with open(init_filename, "r") as f:
        initial_conds =np.loadtxt(f)

    particles = []
    #Each row represents a new particle
    for condition in initial_conds:
        particles.append(Sphere(particle_id=condition[1],pos=condition[2:5],vel=condition[5:8],radius=condition[8], mass=condition[9], ptype=condition[10], material=material_a, rtd2=condition[11:14], rtd3=condition[14:17],rtd4=condition[17:20],force=condition[20:23]))
    return particles






def total_kinetic_energy():
    """Calculates total kinetic energy of all free particles"""
    sum=0
    for particle in particles:
        if particle.ptype == 0:
            sum+=particle[i].kinetic_energy
    return sum





