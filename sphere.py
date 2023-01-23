"""Sphere - code for particle class Sphere"""
import numpy as np


def normalize(dx, L):
    if dx<-L/2:
        dx+=L
    if dx>=L/2:
        dx-=L
    return dx


class Sphere:
    """Sphere particle class

    Instances are created by common.init_system which reads particle
    pos, vel etc from a .dat file. The material type is a dictionary
    imported from materials.py

    particles are updated by:
        1. self.predict yielding new values,
        2. calculating forces
        3. self.correct to create more accurate answer

    pos = [x, y, phi]
    """
    particle_id=0

    def __init__(self,
                pos=np.array([0,0,0]),
                vel=np.array([0,0,0]),
                radius : float =0.0,
                mass : float =0.0,
                ptype : int=0.0,
                material : dict=None,
                rtd2=np.array([0,0,0]),
                rtd3=np.array([0,0,0]),
                rtd4=np.array([0,0,0]),
                force = np.array([0,0,0])):

        Sphere.particle_id += 1
        self.particle_id = Sphere.particle_id
        self.radius = radius
        self.mass = mass
        self.inertia = (2/5)*self.mass * self.radius*self.radius # ???p45 says 2/5MR**2 but code shows 1/2MR**2
        self.ptype = ptype
        self._Y = material['Y_modulus']
        self._A = material['A_damping']
        self._mu = material['mu_coulomb_friction']
        self._gamma = material['gamma_tangential_damping']

        self.pos = pos
        self.vel = vel
        self.rtd2 = rtd2
        self.rtd3 = rtd3
        self.rtd4 = rtd4

        self.force = force

    def predict(self, dt):
        """predict is the first half of the Gear's integration scheme"""
        a1 = dt
        a2 = a1 * dt/2
        a3 = a2 * dt/2
        a4 = a3 * dt/4

        self.pos += a1 * self.vel + a2 * self.rtd2 +a3 * self.rtd3 + a4 * self.rtd4
        self.vel += a1 * self.rtd2 +a2 * self.rtd3 + a3 * self.rtd4
        self.rtd2 += a1 * self.rtd3 + a2 * self.rtd4
        self.rtd3 += a1 * self.rtd4

    def correct(self, dt, force, G):
        """correct is the second half of the Gear's integration scheme"""
        dtrez = 1/dt
        coeff0 = (19/90)*(dt * dt / 2)
        coeff1 = (3/4) * (dt / 2)
        coeff3 = (1/2) * 3 * dtrez
        coeff4 = (1/12) * 12 * dtrez * dtrez

        accel = self._accel(force, G)

        corr = accel - self.rtd2
        self.pos += coeff0 * corr
        self.vel += coeff1 * corr
        self.rtd2 += accel
        self.rtd3 += coeff3 * corr
        self.rtd4 += coeff4 * corr


    def _accel(self, force, G):
        return np.array([
            (1/self.mass)*force.x + G.x,
            (1/self.mass)*force.y + G.y,
            (1/self.inertia)*force.phi + G.phi
        ])

    def boundary_conditions(self, n, timestep, time):
        """Specify the behaviour of particles. Operates
        only on non-standard particles e.g walls etc."""
        if self.ptype == 0:
            #normal dynamic particle
            pass
        elif self.ptype == 1:
            # static particles
            pass
        elif self.ptype == 2:
            # particle oscillates horizontally
            self.pos[0]= 0.5-0.4*np.cos(10*time)
            self.pos[1] = 0.1
            self.vel[0] = 10*0.4*np.sin(time)
            self.vel[1] = 10


    def periodic_bc(self, x_0, y_0, lx, ly):
        if self.pos[0] < x_0:
            self.pos[0] += lx
        elif self.pos[0] > x_0+lx:
            self.pos[0] -= lx
        if self.pos[1] < y_0:
            self.pos[1] += ly
        elif self.pos[1] > y_0+ly:
            self.pos[1] -= ly

    def kinetic_energy(self):
        return self.mass *(self.vel[0] * self.vel[0] / 2 + self.vel[1] * self.vel[1] / 2) + self.inertia * (self.vel[2] * self.vel[2] / 2)





def calc_force(p1, p2, lx, ly):
    dx=normalize(p1.x-p2.x,lx)
    dy=normalize(p1.y-p2.y,ly)
    rr=np.sqrt(dx*dx + dy*dy)
    r1=p1.r
    r2=p2.r
    xi=r1+r2-rr

    if xi > 0:
        Y=p1.Y*p2.Y/(p1.Y + p2.Y)
        A=0.5*(p1.A+p2.A)
        if p1.mu < p2.mu:
            mu=p1.mu
        else:
            mu=p2.mu
        if p1.gamma < p2.gamma:
            gamma = p1.gamma
        else:
            gamma = p2.gamma
        reff = (r1*r2)/(r1+r2)
        dvx = p1.vx - p2.vx
        dvy = p1.vy - p2.vy
        rr_rez = 1 / rr
        ex = dx * rr_rez
        ey = dy * rr_rez
        xidot = -(ex*dvx + ey*dvy)
        vtrel=-dvx*ey + dvy*ex + p1.omega * p1.r - p2.omega * p2.p.r
        fn = np.sqrt(xi) * Y * np.sqrt(reff) * (xi + A*xidot)
        ft=-gamma*vtrel

        if fn < 0:
            fn = 0
        if ft < -mu * fn:
            ft = -mu*fn
        if ft > mu * fn:
            ft = mu*fn
        if p1.type == 0:
            p1.add_force(np.array([fn*ex-ft*ey,fn*ey+ft*ex, r1*ft]))
        if p2.type == 0:
            p2.add_force(np.array([-fn*ex-ft*ey,-fn*ey+ft*ex, -r2*ft]))


