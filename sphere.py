"""Sphere - code for particle class Sphere"""
import numpy as np


def normalize(dx, L):
    while dx<-L/2:
        dx+=L
    while dx>=L/2:
        dx-=L
    return dx


class Sphere:
    def __init__(self, _r, _J, _m, _ptype, material):
        self._r = _r
        self._J = _J
        self._m = _m
        self._ptype = _ptype
        self._Y = material['Y']
        self._A = material['A']
        self._mu = material['mu']
        self._gamma = material['gamma']

        self.rtd0
        self.rtd1
        self.rtd2
        self.rtd3
        self.rtd4

    def predict(self, dt):
        """predict is the first half of the Gear's integration scheme"""
        a1 = dt
        a2 = a1 * dt/2
        a3 = a2 * dt/2
        a4 = a3 * dt/4

        self.rtd0 += a1 * self.rtd1 + a2 * self.rtd2 +a3 * self.rtd3 + a4 * self.rtd4
        self.rtd1 += a1 * self.rtd2 +a2 * self.rtd3 + a3 * self.rtd4
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
        self.rtd0 += coeff0 * corr
        self.rtd1 += coeff1 * corr
        self.rtd2 += accel
        self.rtd3 += coeff3 * corr
        self.rtd4 += coeff4 * corr


    def _accel(self, force, G):
        return np.array([
            (1/self._m)*force.x + G.x, 
            (1/self._m)*force.y + G.y,
            (1/self._J)*force.phi + G.phi
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
            self.x = 0.5-0.4*np.cos(10*time)
            self.y = 0.1
            self.vx = 10*0.4*np.sin(time)
            self.vy = 10

    def periodic_bc(self, x_0, y_0, lx, ly):
        if self.rtd0.x < x_0:
            self.rtd0.x += lx
        elif self.rtd0.x > x_0+lx:
            self.rtd0.x -= lx
        if self.rtd0.y < y_0:
            self.rtd0.y += ly
        elif self.rtd0.y > y_0+ly:
            self.rtd0.y -= ly
    
    def kinetic_energy(self):
        return self._m *(self.rtd1.x * self.rtd1.x / 2 + self.rtd1.y * self.rtd1.y / 2) + self._J * (self.rtd1.phi self.rtd1.phi / 2)
    
    def output(self, f):
        data = np.array([self.rtd0, self.rtd1, self._r, self._m, self._ptype, self._Y, self._A, self._mu, self._gamma, self._force, self.rtd3, self._rtd4])
        np.savetxt("temp.txt",data, "w")
        
        



def force(p1, p2, lx, ly):
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
           

