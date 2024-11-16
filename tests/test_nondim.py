from aldsim.nondim.ideal.particle import MixedFlow, MixedPlugFlow

class TestPlugFlowMixed:

    def test_saturationcurve(self):
        pfm = MixedPlugFlow(10)
        x, y = pfm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = MixedPlugFlow(10)
        x,y,z = pfm.run()
        assert x.shape == y.shape


class TestMixedFlow:

    def test_saturationcurve(self):
        wsm = MixedFlow(10)
        x, y = wsm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = MixedFlow(10)
        x,y,z = pfm.run()
        assert x.shape == y.shape


