from aldsim.models import WellStirred, PlugFlowMixed
import pytest

class TestPlugFlowMixed:

    def test_saturationcurve(self):
        pfm = PlugFlowMixed(10)
        x, y = pfm.saturation_curve()
        assert x.shape == y.shape


