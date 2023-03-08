__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross, Andrea Codd"

from esys.escript import *
import numpy as np
import cmath
from esys.escript.pdetools import Locator, MaskFromTag
from esys.escript.linearPDEs import LinearSinglePDE, SolverOptions

"""
MT 2D
"""

class MT2Dmodel(object):
    def __init__(self, domain, sigma_background, sensors, periods, mu=4*np.pi*1e-7, dd=(1.,1.), origin=(0.,0.),
        fixBottom = False, airLayer = 0.0):
        self.domain = domain
        self.sigmaBG = sigma_background
        self.fixBottom = fixBottom
        self.modelTE = MT2DTE(self.domain, sensors, mu=mu, fixBottom=self.fixBottom)
        self.modelTM = MT2DTM(self.domain, sensors, mu=mu, fixBottom=self.fixBottom, airLayer = airLayer)
        self.dd = dd
        self.periods = periods
        self.loc = Locator(ReducedFunction(self.domain),[s for s in sensors])
        self.numS = len(sensors)
        self.numP = len(periods)
        self.origin=origin

    def setSigmaRhoForward(self, sigma, rho, Ground=1.):
        self.modelTE.setConductivity(sigma, Ground=Ground)
        self.modelTM.setResistivity(rho, Ground=Ground)  
        return self

    def getSigmaRho(self):
        sigma = self.modelTE.getConductivity()
        rho = self.modelTM.getResistivity()
        return sigma, rho

    def getExHxZs(self):         
        Exs=[]
        Hxs=[]
        pZxys = []
        pZyxs = []
        for p in self.periods:
            om = 2.*np.pi/p
            Ex, Zxy = self.modelTE.getExZxy(om)
            Exs.append(Ex)
            pZxys.append(self.loc(Zxy))
            Hx, Zyx = self.modelTM.getHxZyx(om)
            Hxs.append(Hx)
            pZyxs.append(self.loc(Zyx))
        return Exs, Hxs, pZxys, pZyxs
    
    def getZs(self):
        ZComp = []
        for p in self.periods:
            om = 2.*np.pi/p 
            Z = self.modelTE.getZxy(om=om)
            ZComp.append(np.array(self.loc(Z)))
        for p in self.periods:
            om = 2.*np.pi/p 
            Z = self.modelTM.getZyx(om=om)
            ZComp.append(np.array(self.loc(Z)))
        return np.array(ZComp)  

########################   
########################     

 
class MT2DTE(object):
    """
    Solve ONE MT TE problem for frequency  om = 2 pi T,
    where T is the period in seconds.
    """
    def __init__(self, domain, sensors, mu = 4*np.pi*1e-7, fixBottom=False):
        self.domain = domain
        self.mu = mu
        self.sigma = None
        self.sigma_boundary = None
        self.fixBottom = fixBottom
        self.pde = self.__createPDE()
        self.loc = Locator(ReducedFunction(self.domain),[s for s in sensors])

    def __createPDE(self):
        pde = LinearSinglePDE(self.domain, isComplex=True,)
        options = pde.getSolverOptions()
        options.setSolverMethod(SolverOptions.DIRECT)
        pde.setSymmetryOn()
        pde.setValue(A = kronecker(self.domain.getDim()))
        z = self.domain.getX()[self.domain.getDim()-1]
        t = sup(z)
        if self.fixBottom:
            b = inf(z)
            pde.setValue(q = whereZero(z-t) + whereZero(z-b), r = (z-b)/(t-b))
        else:
            pde.setValue(q = whereZero(z-t), r = 1.)
        return pde
    
    def setConductivity(self, sigma, sigma_boundary=None, Ground=1.):
        self.sigma = Ground*interpolate(sigma, Function(self.domain))
        if not self.fixBottom:
            if sigma_boundary:
                self.sigma_boundaryinterpolate(sigma_boundary, FunctionOnBoundary(self.domain))
            else:
                self.sigma_boundary=interpolate(sigma, FunctionOnBoundary(self.domain))
        return self
    
    def getConductivity(self):
        return self.sigma 

    def getEx(self, om=1.): 
        """
        return electric field, Ex (solution space) for om = 2pi/T, T is period in s. 
        :param om: frequency in [Hz]
        :type om: `float`
        :returns: Ex in solution space (continuous function)
        """
        self.pde.setValue(D=1j*om*self.mu*self.sigma)
        if not self.fixBottom:
            z = FunctionOnBoundary(self.domain).getX()[self.domain.getDim()-1]
            b = inf(z)
            k = (1+1j)*sqrt(om*self.mu*self.sigma_boundary/2)
            self.pde.setValue(d = k*whereZero(z-b))      
        return self.pde.getSolution()

    def getZxy(self, om=1.):
        """
        return impedance Zxy (function space) for om = 2pi/T, T is period in s. 
        :param om: frequency in [Hz]
        :type om: `float`
        :returns: Zxy in function space
        """
        Ex = self.getEx(om)
        g=grad(Ex, Function(self.domain))
        Hy=-1./(1j*om*self.mu)*g[self.domain.getDim()-1]
        Exi=interpolate(Ex, Hy.getFunctionSpace())
        Zxy=Exi/Hy
        return Zxy

    def getpZxy(self,om=1):
        Zxy=self.getZxy(om)
        return np.array(self.loc(Zxy)) 

    def getApparentResitivity(self, om, Zxy):
        """
        return the apparent resistivity from a given frequency `f` and impedance `Zxy`

        :param f: frequency in [Hz]
        :type f: `float`
        :param Zxy: impedance
        :type Zxy: `Data` or `np.array`
        """
        return abs(Zxy)**2/(self.mu*om)

    def getPhase(self, Zxy):
        """
        return the phase in [deg] from a given frequency `f` and impedance `Zxy`
        
        :param f: frequency in [Hz]
        :type f: `float`
        :param Zxy: impedance
        :type Zxy: `Data` or `np.array`
        """
        return atan2(np.imag(Zxy),np.real(Zxy))/np.pi*180.

    def getExZxy(self,om=1.):
        """
        return Exi and Zxy (function spaces) for om = 2pi/T, T is period in s. 
        :param om: frequency in [Hz]
        :type om: `float`
        :returns: Exi and Zxy in function space
        """
        Ex = self.getEx(om)        
        g = grad(Ex, Function(self.domain))
        Hy = -1./(1j*om*self.mu)*g[self.domain.getDim()-1]
        Exi = interpolate(Ex, Hy.getFunctionSpace())
        Zxy = Exi/Hy
        return Exi, Zxy



####################################################################
####################################################################
####################################################################

class MT2DTM(object):
    """
    Solve ONE MT TM problem for frequency  om = 2 pi T,
    where T is the period in seconds.
    """
    def __init__(self, domain, sensors, mu=4*np.pi*1e-7, fixBottom=False, airLayer=None):
        """
        """
        self.domain=domain
        self.mu=mu
        self.rho=None
        self.rho_boundary=None
        self.fixBottom=fixBottom
        self.pde=self.__createPDE(airLayer)
        self.loc = Locator(ReducedFunction(self.domain),[s for s in sensors])

    def __createPDE(self, airLayer=None):
        pde=LinearSinglePDE(self.domain, isComplex=True,)
        optionsG=pde.getSolverOptions()
        optionsG.setSolverMethod(SolverOptions.DIRECT)
        pde.setSymmetryOn()
        z=self.domain.getX()[self.domain.getDim()-1]
        b=inf(z)        
        if airLayer is None:
            self.airLayer=whereNonNegative(z-sup(z))
        elif isinstance(airLayer, float) or  isinstance(airLayer, int):
            self.airLayer=whereNonNegative(z-airLayer)
        else:
            self.airLayer=wherePositive(interpolate(airLayer, Solution(self.domain)))
        if self.fixBottom:
            pde.setValue(q=self.airLayer+whereZero(z-b), r=self.airLayer)
        else:
            pde.setValue(q=self.airLayer, r=self.airLayer)
        return pde
    
    def setResistivity(self, rho, rho_boundary=None,Ground=1.):
        self.rho=Ground*interpolate(rho, Function(self.domain))+(1-Ground)*1.e10
        if not self.fixBottom:
            if rho_boundary:
                self.rho_boundary=interpolate(rho_boundary, FunctionOnBoundary(self.domain))
            else:
                self.rho_boundary=interpolate(rho, FunctionOnBoundary(self.domain))

        self.pde.setValue(A=self.rho*kronecker(self.domain.getDim())) 
        return self
    
    def getResistivity(self):
        return self.rho 
    
    def getHx(self, om=1.):
        self.pde.setValue(D=1j*om*self.mu)
        if not self.fixBottom:
            z=FunctionOnBoundary(self.domain).getX()[self.domain.getDim()-1]
            b=inf(z)
            k=(1+1j)*sqrt(om*self.mu*self.rho_boundary/2)
            self.pde.setValue(d=k*whereZero(z-b))
        return self.pde.getSolution()
 
    def getZyx(self, om=1):
        Hx=self.getHx(om) 
        g=grad(Hx, Function(self.domain))
        Ey=self.rho*g[self.domain.getDim()-1]
        Hxi=interpolate(Hx, Ey.getFunctionSpace())
        Zyx=Ey/Hxi
        return Zyx

    def getpZyx(self,om=1):
        Zyx=self.getZyx(om)
        return np.array(self.loc(Zyx)) 

    def getPhase(self, Zyx):
        """
        return the phase in [deg] from a given frequency f [Hz] and impedance Zyx
        :param Zyx: impedance
        :type Zyx: `Data` or `np.array`
        """
        return atan2(np.imag(Zyx),np.real(Zyx))/np.pi*180


    def getHxZyx(self, om=1):
        Hx=self.getHx(om) 
        g=grad(Hx, Function(self.domain))
        Ey=self.rho*g[self.domain.getDim()-1]
        Hxi=interpolate(Hx, Ey.getFunctionSpace())
        Zyx=Ey/Hxi
        return Hxi, Zyx        
    
    def getApparentResitivity(self, om, Zyx):
        """
        return the apparent resistivity from a given frequency `f` and impedance `Zxy`
        :param f: frequency in [Hz]
        :type f: `float`
        :param Zxy: impedance
        :type Zxy: `Data` or `np.array`
        """
        return abs(Zyx)**2/(self.mu*om)
        


