"""Sphere - code for particle class Sphere"""

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

    

