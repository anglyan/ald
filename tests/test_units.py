from aldsim.units import Conversion, temperature

class TestConversion:

    def test_conv(self):
        c = Conversion()
        assert c(1) == 1
        c = Conversion(c_mul=0.5)
        assert c(1) == 0.5
        c = Conversion(c_add=-1)
        assert c(1) == 0

class TestTemperature:

    def test_CtoK(self):
        assert temperature(100) == 373.15

    def test_KtoC(self):
        assert temperature(373.15, units='K', to_units='C') == 100

    def test_KtoK(self):
        assert temperature(300, units='K') == 300





