__copyright__ = "Copyright (c) 2020 by University of Queensland http://www.uq.edu.au"
__license__   = "Licensed under the Apache License, version 2.0 http://www.apache.org/licenses/LICENSE-2.0"
__credits__   = "Lutz Gross, Andrea Codd"

from esys.escript import *
import numpy as np
from esys.escript.pdetools import Locator, MaskFromTag
from MTApps.MTForApp import MT2DTE, MT2DTM
from esys.weipa import saveSilo

"""
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


