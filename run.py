from common import init_system, init_algorithm

import json


def main(init_file, sim_params):
    """This is the main function that intialises and runs the simulation loop"""

    #Setup a vector of particles from initialisation file.
    particles = init_system(init_file)
    init_algorithm()
    phase_plot()
    for i in range(nstep):
        step()
        if (i+1)%nprint == 0:
            phase_plot(fphase)
        if (i+1)%nenergy == 0:
            total_kinetic_energy(fenergy)



if __name__ == '__main__':
    
    


    main(sim_name)