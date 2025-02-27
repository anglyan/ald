from aldsim.core.ideal.particle.continuous import PlugFlowSpatial, WellMixedSpatial

class TestPlugFlowSpatial:

    def test_saturationcurve(self):
        pfm = PlugFlowSpatial(10)
        x, y = pfm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = PlugFlowSpatial(10)
        x,y,z = pfm.run()
        assert x.shape == y.shape


class TestWellMixedSpatial:

    def test_saturationcurve(self):
        wsm = WellMixedSpatial(10)
        x, y = wsm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = WellMixedSpatial(10)
        x,y,z = pfm.run()
        assert x.shape == y.shape


