#Copyright © 2024, UChicago Argonne, LLC
"""Utility functions to transform physical quantities to the right units"""

from .constants import Nav

class Conversion:

    def __init__(self, c_add=None, c_mul=None, inv=False):
        self.c_add = c_add
        self.c_mul = c_mul
        self.inv = False

    def __call__(self, value):

        c = 1.0/value if self.inv else value
        c = value if self.c_mul is None else self.c_mul * value
        if self.c_add is None:
            return c
        else:
            return c + self.c_add
        
    def inverse(self):
        if self.c_add is None:
            if self.c_mul is None:
                return Conversion()
            else:
                return Conversion(c_mul=1.0/self.c_mul, inv=self.inv)
        else:
            if self.inv:
                raise ValueError("Cannot compute inverse when inv is True")
            else:
                return Conversion(c_mul=1.0/self.c_mul,
                                c_add=-self.c_add/self.c_mul)


_t_conv = {
    ('C','K') : Conversion(c_add=273.15),
    ('K','C') : Conversion(c_add=-273.15)
}

_area_conv = {
    ('nm2','m2') : Conversion(c_mul=1e-18),
    ('A2', 'm2') : Conversion(c_mul=1e-20),
    ('cm2', 'm2') : Conversion(c_mul=1e-4),
    ('mol/m2', 'm2') : Conversion(c_mul = 1.0/Nav, inv=True)
}

def temperature(value, units='C', to_units='K'):
    if units == to_units:
        return value
    else:
        return _t_conv[(units, to_units)](value)

def site_area(value, units='nm2', to_units='m2'):
    if units == to_units:
        return value
    else:
        return _area_conv[(units, to_units)](value)


