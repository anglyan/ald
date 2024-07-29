import numpy as np

class PlugFlowMixed:

    def __init__(self, Da):
        self.Da = Da

    def calc_coverage(self, t):
        Da = self.Da
        return 1 - 1/Da*np.log(1+(np.exp(Da)-1)*np.exp(-Da*t))
    
    def saturation_curve(self, tmax=5, dt= 0.01):
        t = np.arange(0, tmax, dt)
        Da = self.Da
        c = 1 - 1/Da*np.log(1+(np.exp(Da)-1)*np.exp(-Da*t))
        return t, c
    
    def run(self, tmax=5, dt=0.01):
        t = np.arange(0, tmax, dt)
        Da = self.Da
        y = 1/Da*np.log(1+(np.exp(Da)-1)*np.exp(-Da*t))
        c = 1-y
        x = np.exp(-Da*y)
        return t, c, x


