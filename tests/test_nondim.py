from aldsim.core.ideal.batch import MixedFlow, PlugFlowMixed

class TestPlugFlowMixed:

    def test_saturationcurve(self):
        pfm = PlugFlowMixed(10)
        x, y = pfm.saturation_curve()
        assert x.shape == y.shape

    def test_run(self):
        pfm = PlugFlowMixed(10)
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


