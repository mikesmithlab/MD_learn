import numpy as np


def dump_particle(f, tstep, p):
    np.savetxt(f, [[tstep, p.particle_id, p.pos[0], p.pos[1], p.pos[2], p.vel[0], p.vel[1], p.vel[2], p.radius, p.mass, p.ptype, p.rtd2[0], p.rtd2[1], p.rtd2[2], p.rtd3[0], p.rtd3[1], p.rtd3[2], p.rtd4[0], p.rtd4[1], p.rtd4[2], p.force[0], p.force[1], p.force[2]]])

def record_state(f, tstep, particles):
    for particle in particles:
        dump_particle(f, tstep, particle)