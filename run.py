from common import init_system, init_algorithm

import json


def main(input_files, output_files):
    """This is the main function that intialises and runs the simulation loop"""
    with open(sim_file,"r") as f_sim_params:
        sim_params = json.load(f_sim_params)
        
    with open(sim_output_file, "a") as f_out:
    
    #Setup a vector of particles from initialisation file.
    particles = init_system(particle_file)
    #init_algorithm()
    #phase_plot()
    for i in range(sim_params['nstep']):
        #step()
        if (i+1)%sim_params['nprint'] == 0:
            print('print')
            #phase_plot(fphase)
        if (i+1)%sim_params['nenergy'] == 0:
            print('energy')
            #total_kinetic_energy(fenergy)
    


if __name__ == '__main__':
    
    init_path = 'sim_data/'    
    init_file = 'init20230123_105242'
    output_path = 'output_data/'

    main(input_files = init_path + init_file)