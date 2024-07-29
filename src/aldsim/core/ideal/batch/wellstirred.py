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
        return -y/(1/self.Da+y)

    def run(self, tmax=5, dt=0.01):
        out = solve_ivp(self._f, [0,tmax], [1], method='LSODA',
                t_eval=np.arange(0,tmax,dt))
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


