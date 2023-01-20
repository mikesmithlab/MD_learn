from common import init_system, init_algorithm

import json


def main(particle_file, sim_file):
    """This is the main function that intialises and runs the simulation loop"""
    with open(sim_file,"r") as f_sim:
        sim_params = json.load(f_sim)
    
    #Setup a vector of particles from initialisation file.
    particles = init_system(particle_file)
    #init_algorithm()
    #phase_plot()
    for i in range(sim_params['nstep']):
        print
        #step()
        if (i+1)%sim_params['nprint'] == 0:
            print('print')
            #phase_plot(fphase)
        if (i+1)%sim_params['nenergy'] == 0:
            print('energy')
            #total_kinetic_energy(fenergy)
    


if __name__ == '__main__':
    
    


    main('Sim_Data/init.dat', 'sim_init.json')