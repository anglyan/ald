#Copyright © 2024, UChicago Argonne, LLC

class Precursor:

    def __init__(self, name='None', M=100, ligands=None):
        self.name = name
        self.M = M
        self.ligands=ligands


class ALDKinetics:

    name = 'ideal'

    def __init__(self, nsites, beta0, f=1):
        self.f = f
        self.nsites = nsites
        self.beta = beta0

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
        


class SoftSaturating:

    name = 'softsat'

    def __init__(self, nsites, beta1, beta2, f1, f2=None):
        self.nsites = nsites
        self.beta1 = beta1
        self.beta2 = beta2
        self.f1 = f1
        if self.f2 is None:
            self.f2 = 1-self.f1
        self.f = self.f1 + self.f2
        self._s0 = self.f/self.nsites


class ALDChem:

    def __init__(self, prec, kinetics, dm=None):
        self.prec = prec
        self.kinetics = kinetics
        if dm is None:
            self.default_dm = True
            self.dm = 1
        else:
            self.dm = dm




