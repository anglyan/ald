from aldsim.chem import Precursor, SurfaceKinetics, ALDsoft, ALDideal
import pytest


def test_precursor():
    c = Precursor(mass=100)
    assert c.name == 'None'
    assert c.mass == 100
    c = Precursor(name='TMA', mass=101.5)
    assert c.name == 'TMA'
        

def test_surfacekinetics():
    p = Precursor(mass=100)
    k = SurfaceKinetics(p, 1e19, 1)
    assert k.site_area == pytest.approx(1e-19)

def test_surfacekinupdate():
    p = Precursor(mass=100)
    k = SurfaceKinetics(p, 1e19, 1)
    k.site_area = 1e-18
    assert k.nsites == pytest.approx(1e18)

def test_aldideal():
    p = Precursor(mass=100)
    ald = ALDideal(p, 1e19, 0.001)


class TestSoftSaturating:

    def test_nsites_s0(self):
        p = Precursor(mass=100)
        k = ALDsoft(p, 1e19, 1e-2, 1e-3, 0.8, 0.2)
        assert k.site_area == pytest.approx(1e-19)
        k.site_area = 1e-18
        assert k.nsites == pytest.approx(1e18)
