import numpy as np

def dump_particle(f,tstep,x,y,z,vx,vy,vz, rad, mass, ptype, rtd2x,rtd2y,rtd2z,rtd3x,rtd3y,rtd3z, rtd4x,rtd4y,rtd4z,fx,fy,fz):
    np.savetxt(f, [[tstep, x, y, z, vx, vy, vz, rad, mass, ptype,rtd2x,rtd2y,rtd2z,rtd3x,rtd3y,rtd3z, rtd4x,rtd4y,rtd4z,fx,fy,fz]])

def write_sim_initfile(pathname='Sim_Data/', filename='init.dat'):
    with open(pathname + filename, "w") as f:
        """Boundary Particles"""
        for i in range(11):
            dump_particle(f, 0,0.45+i*0.01, 0.19,0, 0, 0, 0, 0.005, 1, 1, 0,0,0,0,0,0,0,0,0,0,0,0)
        for i in range(50):
            dump_particle(f, 0,0.55+(i+0.5)*0.005, i*0.01+0.2,0, 0, 0, 0, 0.005, 1, 1, 0,0,0,0,0,0,0,0,0,0,0,0)
            dump_particle(f, 0,0.45-(i+0.5)*0.005, i*0.01+0.2, 0, 0, 0, 0, 0.005, 1, 1,0,0,0,0,0,0,0,0,0,0,0,0)

        """free particles"""
        Rmax=0.006
        Rmin=0.004
        for i in range(67):
            for k in range(30):
                centerx=0.3115+0.013*k
                centery=0.6+0.013*i
                variation=np.random.uniform()
                r=Rmin*Rmax/(Rmax-variation*(Rmax-Rmin))
                dump_particle(f, 0, centerx, centery,0, 0, 0, 0, r, r*r/(Rmax*Rmax), 0, 0,0,0,0,0,0,0,0,0,0,0,0)


if __name__=='__main__':
    write_sim_initfile()
