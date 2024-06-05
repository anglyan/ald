from aldsim.models import WellStirred, PlugFlowMixed

class TestPlugFlowMixed:

    def test_saturationcurve(self):
        pfm = PlugFlowMixed(10)
        x, y = pfm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = PlugFlowMixed(10)
        x,y,z = pfm.run()
        assert x.shape == y.shape


class TestWellStirred:

    def test_saturationcurve(self):
        wsm = WellStirred(10)
        x, y = wsm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = WellStirred(10)
        x,y,z = pfm.run()
        assert x.shape == y.shape


