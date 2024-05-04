"""Utility functions to transform physical quantities to the right units"""

class Conversion:

    def __init__(self, c_add=None, c_mul=None):
        self.c_add = c_add
        self.c_mul = c_mul

    def __call__(self, value):

        c = value if self.c_mul is None else self.c_mul * value
        if self.c_add is None:
            return c
        else:
            return c + self.c_add


_t_conv = {
    ('C','K') : Conversion(c_add=273.15),
    ('K','C') : Conversion(c_add=-273.15)
}

def temperature(value, units='C', to_units='K'):
    if units == to_units:
        return value
    else:
        return _t_conv[(units, to_units)](value)



