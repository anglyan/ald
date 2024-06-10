
class IdealDoseModel:

    def __init__(self, chem, p, T):
        self.chem = chem
        self.p = p
        self.T = T

    @property
    def T(self):
        return self._T
 
    @T.setter
    def T(self, value):
        self._vth = self.chem.vth(value)
        self._T = value

    @property
    def vth(self):
        return self._vth

    @property
    def site_area(self):
        return self.chem.site_area
    
    @site_area.setter
    def site_area(self, value):
        self.chem.site_area = value
    
    @property
    def mass(self):
        return self.prec.mass
   
    @property
    def p(self):
        return self._p
 
    @p.setter
    def p(self, value):
        self._p = value