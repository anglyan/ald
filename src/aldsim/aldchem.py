#Copyright © 2024, UChicago Argonne, LLC

from .aldutils import calc_vth
from .constants import kb
from .core.ideal import ZeroD

import numpy as np

class Precursor:
    """
    Implement a precursor molecule
    """

    def __init__(self, name='None', mass=100, ligands=None):
        self.name = name
        self.mass = mass
        self.ligands = ligands

    def vth(self, T):
        """Calculate the mean thermal velocity at temperature T (in K)"""
        return calc_vth(self.mass, T)


class SurfaceKinetics:

    def __init__(self, prec, nsites, f=1):
        self.prec = prec
        self.f = f
        self.nsites = nsites

    @property
    def site_area(self):
        return self._s0
    
    @site_area.setter
    def site_area(self, value):
        self._s0 = value
        self._nsites = self.f/self._s0
    
    @property
    def nsites(self):
        return self._nsites
    
    @nsites.setter
    def nsites(self, value):
        self._nsites = value
        self._s0 = self.f/self._nsites

    def beta(self, *args):
        pass

    def beta_av(self, *args):
        pass

    def vth(self, T):
        return self.prec.vth(T)

        

class ALDideal(SurfaceKinetics):
    """Ideal first-order irreversible Langmuir kinetics"""

    name = 'ideal'

    def __init__(self, prec, nsites, beta0, f=1, dm=1):
        self.beta0 = beta0
        self.dm = dm
        super().__init__(prec, nsites, f)

    def beta(self, cov=0):
        return self.f*self.beta0*(1-cov)

    def beta_av(self, av):
        return self.f*self.beta0*av
    
    def model(self, model_name, **kwargs):
        if model_name.lower() == 'zerod':
            return ZeroD(self.prec, self.nsites, self.beta0, self.f, **kwargs)
        


class SoftSaturating(SurfaceKinetics):
    """First-order irreversible Langmuir kinetics with two reaction pathways"""

    name = 'softsat'

    def __init__(self, prec, nsites, beta1, beta2, f1, f2=None):
        self.beta1 = beta1
        self.beta2 = beta2
        self.f1 = f1
        if f2 is None:
            self.f2 = 1-self.f1
        else:
            self.f2 = f2
        super().__init__(prec, nsites, self.f1+self.f2)

    def beta(self, cov1=0, cov2=0):
        return self.f1*self.beta1*(1-cov1) + self.f2*self.beta2*(1-cov2)

    def beta_av(self, av1, av2):
        return self.f1*self.beta1*av1 + self.f2*self.beta2*av2


class ALDProcess:

    def __init__(self, chem1, chem2, n1, n2, T=None):

        self.chem1 = chem1
        self.n1 = n1
        self.n2 = n2
        self.chem2 = chem2
        self.T = T

    @property
    def T(self):
        return self._T
    
    @T.setter
    def T(self, value):
        self._T = value
        self.chem1.T = value
        self.chem2.T = value

    
    


    





