#Copyright © 2024-2025, UChicago Argonne, LLC

from .aldutils import calc_vth
from .constants import kb, Rgas, Nav

import numpy as np

_precursor_mass = {
    'TMA' : 144.17,
    'H2O' : 18.01
}

class Precursor:
    """
    Defines a precursor molecule
    """

    def __init__(self, name='None', mass=None, ligands=None):
        self.name = name
        if mass is None:
            self.mass = _precursor_mass[self.name]
        else:
            self.mass = mass
        self.ligands = ligands

    def vth(self, T):
        """Calculate the mean thermal velocity at temperature T (in K)"""
        return calc_vth(self.mass, T)
    
    def Jwall(self, T, p, in_mols=False):
        """Calculate the flux per unit area for a given temperature (in K) and pressure (in Pa)"""
        if in_mols:
            return 0.25*self.vth(T)*p/(Rgas*T)
        else:
            return 0.25*self.vth(T)*p/(kb*T)


class SurfaceKinetics:
    """
    Base class for self-limited kinetics

    It assumes that a fraction of the surface is reactive and comprised of reaction
    sites of equal surface area. The default is that all the surface is reactive.

    """

    def __init__(self, prec, nsites, f=1):
        self.prec = prec
        self._f = f
        self.nsites = nsites

    @property
    def site_area(self):
        """Area of a single reaction site"""
        return self._s0
    
    @site_area.setter
    def site_area(self, value):
        self._s0 = value
        self._nsites = self._f/self._s0
    
    @property
    def nsites(self):
        """Number of reactive sites per surface area"""
        return self._nsites

    @nsites.setter
    def nsites(self, value):
        self._nsites = value
        self._s0 = self._f/self._nsites

    @property
    def nsites_mol(self):
        """Number of reactive sites per surface area in mols"""
        return self.nsites/Nav
    
    @property
    def f(self):
        """Fraction of reactive sites"""
        self._f 

    @f.setter
    def f(self, value):
        self._f = value
        self._nsites = self._f/self._s0

    def beta(self, *args):
        pass

    def beta_av(self, *args):
        pass

    def vth(self, T):
        return self.prec.vth(T)
    
    def Jwall(self, T, p):
        return self.prec.Jwall(T, p)

        

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
    
    def t0(self, T, p):
        """Characteristic time for saturation"""
        return 1.0/(self.site_area*self.Jwall(T, p)*self.beta0)
    
    def saturation_curve(self, T, p):
        """Return the saturation curve as a (time, coverage) tuple """
        t0 = self.t0(T,p)
        tscale = 5*t0
        logtscale = np.log10(tscale)
        scale = int(logtscale)
        if logtscale < 0:
            scale -= 1
        factor = int(10**(logtscale-scale))+1
        tmax= factor*10**scale
        print(t0, tmax)
        dt = tmax/100
        x = np.arange(0, tmax, dt)
        y = 1-np.exp(-x/t0)
        return x, y


class ALDsoft(SurfaceKinetics):
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

    def t0(self, T, p):
        """Characteristic time for saturation"""
        t1 = 1.0/(self.site_area*self.Jwall(T, p)*self.beta1)
        t2 = 1.0/(self.site_area*self.Jwall(T, p)*self.beta2)
        return t1, t2
    
    def saturation_curve(self, T, p):
        """Return the saturation curve as a (time, coverage) tuple """
        t1, t2 = self.t0(T,p)
        t0 = max(t1, t2)
        tscale = 5*t0
        logtscale = np.log10(tscale)
        scale = int(logtscale)
        if logtscale < 0:
            scale -= 1
        factor = int(10**(logtscale-scale))+1
        tmax= factor*10**scale
        print(t0, tmax)
        dt = tmax/100
        x = np.arange(0, tmax, dt)
        y = (self.f1*(1-np.exp(-x/t1)) + self.f2*(1-np.exp(-x/t2)))/(self.f1+self.f2)
        return x, y



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

    
    


    





