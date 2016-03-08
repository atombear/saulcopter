from __future__ import division
import scipy
import scipy.integrate

class flyer(object):
    
    g = 9.8 # gravity m/s^2
    
    def __init__(self, L, R, Omega, ramp = 'none', rt = 0):
        self.L = L # length of rope m
        self.R = R # radius of circular trajectory m
        self.Omega = Omega # angular frequency of circular trajectory m
        self.ramp = ramp # How the radius should change to its ss value
        self.rt = rt # How long a linear ramp should take
        
        self.OmegaR = (R*Omega**2/L)**.5
        self.Omega0 = (self.g/L)**.5
        
    def OmegaRamp(self,t0):
        if self.ramp == 'linear':
            return self.OmegaR*(t0>self.rt) + self.OmegaR/self.rt*t0*(t0<self.rt)
        
    def deriv(self,y,t0):
        dydt = scipy.zeros(4)
        
        if self.ramp == 'none':
            dydt[0] = y[1]
            dydt[1] = self.OmegaR**2*scipy.cos(y[0])*scipy.cos(self.Omega*t0-y[2]) - self.Omega0**2*scipy.sin(y[0]) + scipy.cos(y[0])*scipy.sin(y[0])*y[3]**2
            dydt[2] = y[3]
            dydt[3] = (self.OmegaR**2*scipy.sin(self.Omega*t0-y[2]) - 2*scipy.cos(y[0])*y[1]*y[3])/scipy.sin(y[0])
        if self.ramp == 'linear':
            OmegaRt0 = self.OmegaRamp(t0)
            dydt[0] = y[1]
            dydt[1] = OmegaRt0**2*scipy.cos(y[0])*scipy.cos(self.Omega*t0-y[2]) - self.Omega0**2*scipy.sin(y[0]) + scipy.cos(y[0])*scipy.sin(y[0])*y[3]**2
            dydt[2] = y[3]
            dydt[3] = (OmegaRt0**2*scipy.sin(self.Omega*t0-y[2]) - 2*scipy.cos(y[0])*y[1]*y[3])/scipy.sin(y[0])
        
        return dydt
        
    def FlightPath(self,y0,t):
        
        if self.ramp == 'none':
            self.varsolved = scipy.integrate.odeint(self.deriv,y0,t)
            self.xpath = self.R*scipy.cos(self.Omega*t)+self.L*scipy.cos(self.varsolved[:,2])*scipy.sin(self.varsolved[:,0])
            self.ypath = self.R*scipy.sin(self.Omega*t)+self.L*scipy.sin(self.varsolved[:,2])*scipy.sin(self.varsolved[:,0])
            self.zpath = -self.L*scipy.cos(self.varsolved[:,0])
            
        if self.ramp == 'linear':
            self.varsolved = scipy.integrate.odeint(self.deriv,y0,t)
            self.xpath = (self.L*self.OmegaRamp(t)**2/self.Omega**2)*scipy.cos(self.Omega*t)+self.L*scipy.cos(self.varsolved[:,2])*scipy.sin(self.varsolved[:,0])
            self.ypath = (self.L*self.OmegaRamp(t)**2/self.Omega**2)*scipy.sin(self.Omega*t)+self.L*scipy.sin(self.varsolved[:,2])*scipy.sin(self.varsolved[:,0])
            self.zpath = -self.L*scipy.cos(self.varsolved[:,0])