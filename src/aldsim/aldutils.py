import numpy as np

kb = 1.38e-23
amu = 1.660e-27

def calc_vth(M,T):
    return np.sqrt(8*kb*T/(np.pi*amu*M))


def calc_sitearea_fromgpc(M, density, gpc, nmol=1):
    """Average area of a surface site from growth per cycle

    Calculate the average area of a surface site

    Parameters
    ----------

    M : float
        Molecular mass in atomic mass units
    density : float
        Density of the film, in g/cm3
    gpc : float
        Growth per cycle, in Angstroms
    nmol : int, optional (default 1)
        Number of precursor molecules per unit formula of the solid

    Returns
    -------
    float
        Average area of a surface site in sq. meters

    """

    masscm2 = density*gpc*1e-8
    molcm2 = masscm2/M*6.022e23
    return 1e-4/(nmol*molcm2)


def calc_sitearea_fromqcm(M, mpc, nmol=1):
    """Average area of a surface site from QCM mass

    Calculate the average area of a surface site from qcm data

    Parameters
    ----------

    M : float
        Molecular mass in atomic mass units
    mpc : float
        Mass per cycle in  ng/cm2
    nmol : int, optional (default 1)
        Number of precursor molecules per unit formula of the solid

    Returns
    -------
    float
        Average area of a surface site in sq. meters

    """

    return M/(mpc*1e-5*6.022e23*nmol)

def calc_sitearea_fromrbs(atoms_area, atoms_permol=1.0):
    """Average area of a surface site from RBS data

    Calculate the average area of a surface site from RBS data

    Parameters
    ----------

    atoms_area : float
        Atoms per unit area for one ALD cycle (atoms per sq. meter)
    atoms_permol : int, optional (default 1)
        Number of atoms per precursor molecule

    Returns
    -------
    float
        Average area of a surface site in sq. meters

    """

    return atoms_permol/atoms_area



class ALDchem:

    """
    Ideal self-limited process

    Codifies all the parameters required to work with self-limited processes

    All parameters are in SI units.


    """

    def __init__(self, p, M, beta0, s0, T):
        self._p = p
        self._M = M
        self._beta0 = beta0
        self._s0 = s0
        self._T = T
        self._update()


    def _update(self):
        self.vth = calc_vth(self.M, self.T)
        self.alpha = 0.25*self.vth*self.beta0
        self.n0 = self.p/(kb*self.T)
        self.nu0 = self.alpha*self.s0*self.n0

  
    @property
    def s0(self):
        return self._s0

    @s0.setter
    def s0(self, s0):
        self._s0 = s0
        self._update()

    @property
    def M(self):
        return self._M

    @M.setter
    def M(self, M):
        self._M = M
        self._update()

    @property
    def beta0(self):
        return self._beta0

    @beta0.setter
    def beta0(self, beta0):
        self._beta0 = beta0
        self._update()
    
    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, T):
        self._T = T
        self._update()

    @property
    def p(self):
        return self._p

    @p.setter
    def p(self, p):
        self._p = p
        self._update()

    def from_qcm(self, mpc, nmol=1):
        self.mpc = mpc
        self.s0 = calc_sitearea_fromqcm(self.M, mpc, nmol)


