#Copyright © 2024, UChicago Argonne, LLC

"""Batch models"""

import numpy as np
from scipy.integrate import solve_ivp

class WellStirred:

    def __init__(self, Da):
        self.Da = Da

    def calc_coverage(self, tau):
        ep = 1
        t = 0.5
        Da = self.Da
        damp = 0.1
        while ep > 1e-6:
            f_t = self._f_t(t, Da, tau)
            fp_t = self._fp_t(t, Da, tau)
            damp = 0.1
            tn = 1.1
            while tn > 1 or tn < 0:
                damp = 0.5*damp
                tn = t - damp*f_t/fp_t
            ep = abs(t-tn)/t
            t = tn
        return t
    
    def _f_t(self, theta, Da, tau):
        return theta - np.log(1-theta)/Da - tau

    def _fp_t(self, theta, Da, tau):
        return 1 + 1/(Da*(1-theta))

    def _f(self, t, y):
        return -self.Da*y/(1+self.Da*y)

    def run(self, tmax=5, dt=0.01):
        out = solve_ivp(self._f, [0,tmax], [1], t_eval=np.arange(0,tmax,dt))
        cov = 1-out.y[0,:]
        x = 1/(1+self.Da*out.y[0,:])
        return out.t, cov, x

    def saturation_curve(self, tmax=5, dt=0.01):
        t, cov, _ = self.run(tmax, dt)
        return t, cov
    
    def saturation_curve_implicit(self, theta_max=0.9999):
        Da = self.Da
        theta = np.arange(0,theta_max,0.0001)
        tau = theta - np.log(1-theta)/Da
        return tau, theta

    def fraction_out(self, theta_max=0.999):
        Da = self.Da
        theta = np.arange(0,theta_max,0.0001)
        tau = theta - np.log(1-theta)/Da
        return tau, 1/(1+Da*(1-theta))


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



class PlugFlowMixedSoft:

    def __init__(self, D1, D2, f1, f2):
        self.D1 = D1
        self.D2 = D2
        self.f1 = f1
        self.f2 = f2
        self.a = D2/D1

    def _f(self, t, y):
        """Provides the gradient for the fraction of available sites
        of the first reaction pathway"""

        dec = self.f1*self.D1*y+self.f2*self.D2*np.float_power(y, self.a)
        return - self.D1*y*(1-np.exp(-dec))/dec

    def saturation_curve(self, tmax=5):

        out = solve_ivp(self._f, [0,tmax], [1], t_eval=np.arange(0,tmax,0.01))
        y1 = out.y[0,:]
        y2 = np.float_power(y1, self.a)
        cov = self.f1*(1-y1) + self.f2*(1-y2)
        return out.t, cov

    def run(self, tmax, dt):
        raise NotImplementedError
    
    def calc_coverage(self, t):
        raise NotImplementedError
    
