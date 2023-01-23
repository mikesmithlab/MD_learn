from common import init_system, init_algorithm, step
from input_output_file import record_state

import json


def main(input_path=None, input_file=None, output_path=None):
    """This is the main function that intialises and runs the simulation loop"""

    #Read simulation parameters and write copy to new file
    with open(input_path + input_file + '.json',"r") as f_sim_params:
        sim_params = json.load(f_sim_params)
        with open(output_path + input_file + '.json',"w") as f_sim_params_out:
            json.dump(sim_params, f_sim_params_out)

    dt = sim_params['timestep']
    nstep = sim_params['nstep']
    nprint = sim_params['nprint']
    time=0

    with open(output_path + input_file + '.dat', "a") as f_out:
        #Setup a vector of particles from initialisation file.
        particles = init_system(init_filename = input_path + input_file + '.dat')
        record_state(f_out, 0, particles)
        #init_algorithm()
        #phase_plot()
        for i in range(nstep):
            if (i+1)%nprint == 0:
                record_state(f_out, i, step(dt, particles))
                #if (i+1)%sim_params['nenergy'] == 0:
                #    print('energy')
                #    #total_kinetic_energy(fenergy)
            else:
                step(dt, particles)


if __name__ == '__main__':

    init_path = 'sim_input_data/'
    init_file = 'init20230123_170608'
    output_path = 'sim_output_data/'

    main(input_path=init_path, input_file=init_file, output_path=output_path)