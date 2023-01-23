import numpy as np
from pydates.pydates import now, format_datetime_to_str
import shutil
from sphere import Sphere
from materials import material_a
from input_output_file import record_state


def write_sim_initfile(sim_param_path='sim_param_templates/',sim_params='sim_init.json', pathname='sim_data/', filename='init.dat'):
    ext= filename.split('.')[1]
    filename = filename.split('.')[0] + format_datetime_to_str(now(), "%Y%m%d_%H%M%S") + '.'
    filename_params = filename + 'json'
    filename_init = filename + ext


    particles=[]
    vel=np.array([0, 0, 0]) #[x,y,phi]
    ptype=0
    radius=0.005
    mass=0
    rtd2=np.array([0,0,0])
    rtd3=np.array([0,0,0])
    rtd4=np.array([0,0,0])
    force = np.array([0,0,0])


    """Boundary Particles"""
    for i in range(11):
        pos=np.array([0.45+i*0.01, 0.19,0,])#[x,y,phi]
        particles.append(Sphere(pos=pos, vel=vel, radius=radius, mass=mass, ptype=ptype, material=material_a, rtd2=rtd2, rtd3=rtd3, rtd4=rtd4, force=force))
    for i in range(50):
        pos=np.array([0.55+(i+0.5)*0.005, i*0.01+0.2,0])
        particles.append(Sphere(pos=pos, vel=vel, radius=radius, mass=mass, ptype=ptype, material=material_a, rtd2=rtd2, rtd3=rtd3, rtd4=rtd4, force=force))
        pos=np.array([0.45-(i+0.5)*0.005, i*0.01+0.2,0])
        particles.append(Sphere(pos=pos, vel=vel, radius=radius, mass=mass, ptype=ptype, material=material_a, rtd2=rtd2, rtd3=rtd3, rtd4=rtd4, force=force))

    """free particles"""
    Rmax=0.006
    Rmin=0.004
    for i in range(67):
        for k in range(30):
            centerx=0.3115+0.013*k
            centery=0.6+0.013*i
            pos=np.array([centerx, centery, 0])
            variation=np.random.uniform()
            r=Rmin*Rmax/(Rmax-variation*(Rmax-Rmin))
            particles.append(Sphere(pos=pos, vel=vel, radius=r, mass=mass, ptype=ptype, material=material_a, rtd2=rtd2, rtd3=rtd3, rtd4=rtd4, force=force))

    with open(pathname + filename_init, "w") as f:
        record_state(f, 0, particles)

    print(sim_param_path)
    print(sim_params)
    print(pathname)
    print(filename_params)
    shutil.copy(sim_param_path + sim_params, pathname + filename_params)

if __name__=='__main__':
    sim_params='sim_init.json'
    write_sim_initfile(sim_params=sim_params)

