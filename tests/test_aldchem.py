from aldsim.aldchem import Precursor, ALDKinetics, SoftSaturating, ALDChem
import pytest

class TestPrecursor:

    def test_precursor(self):
        c = Precursor()
        assert c.name == 'None'
        assert c.mass == 100


    def test_customprecursor(self):
        c = Precursor(name='TMA', mass=101.5)
        assert c.name == 'TMA'
        

class TestALDKinetics:

    def test_nsites_s0(self):
        k = ALDKinetics(1e19, 1e-2)
        assert k.site_area == pytest.approx(1e-19)
        k.site_area = 1e-18
        assert k.nsites == pytest.approx(1e18)


class TestSoftSaturating:

    def test_nsites_s0(self):
        k = SoftSaturating(1e19, 1e-2, 1e-3, 0.8, 0.2)
        assert k.site_area == pytest.approx(1e-19)
        k.site_area = 1e-18
        assert k.nsites == pytest.approx(1e18)


class TestALDChem:

    def test_setT(self):
        k = ALDKinetics(1e19, 1e-2)
        c = Precursor(name='TMA', mass=101.5)
        ac = ALDChem(c, k)
        ac.T = 300
        assert ac.T == 300

    def test_sitearea(self):
        k = ALDKinetics(1e19, 1e-2)
        c = Precursor(name='TMA', mass=101.5)
        ac = ALDChem(c, k)
        assert ac.site_area == pytest.approx(1e-19)
        ac.site_area = 1e-18
        assert ac.kinetics.site_area == pytest.approx(1e-18)
        

