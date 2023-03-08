from esys.escript import *
import numpy as np
from esys.escript.pdetools import Locator, MaskFromTag
from MTApps.MTForApp import MT2DTE, MT2DTM
from esys.weipa import saveSilo

#from esys.downunder.apps import MT2DTEModel, MT2DTMModel
#from MTApps.MTInvApp import MT2DTE, MT2DTM

#from math import floor
#from scipy.interpolate import RegularGridInterpolator

"""
"""

class MT2Dmodel(object):
    def __init__(self, domain, sigma_background, sensors, periods, mu=4*np.pi*1e-7, dd=(1.,1.), origin=(0.,0.),
        fixBottom = False, airLayer = 0.0):
        print("model2D init")
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
        #xyz=ReducedFunction(self.domain).getX()
        #self.GrF = interpolate(whereNegative(xyz[1]-airLayer), Function(self.domain))        
        #self.GrF = Scalar(0.,Function(domain))
        #self.GrF.setTaggedValue("core",1.)
        #self.GrF.setTaggedValue("buffer",1.)
        #self.GrF.expand()


    def setSigmaRho(self,sigmaf):
        Rcore = Scalar(0,Function(self.domain))
        Rcore.setTaggedValue("core",1)
        Rcore.expand()
        Rbuffer = Scalar(0,Function(self.domain))
        Rbuffer.setTaggedValue("buffer",1)
        Rbuffer.expand()
        Rair = Scalar(0,Function(self.domain)) 
        Rair.setTaggedValue("air",1)
        Rair.expand()
        X=Function(self.domain).getX()
        sigmaI = Rcore*interpolateTable(sigmaf, X[:2], self.origin, self.dd)+1.e6*(1-Rcore)
        sigma = sigmaI*Rcore + Rbuffer*self.sigmaBG + Rair*0.0
        rho = (1./sigmaI)*Rcore + Rbuffer/self.sigmaBG+ Rair*1.e6
        self.modelTE.setConductivity(sigma)
        self.modelTM.setResistivity(rho) 
        return self 

    def setSigmaRhoForward(self, sigma, rho, Ground=1.):
        self.modelTE.setConductivity(sigma, Ground=Ground)
        self.modelTM.setResistivity(rho, Ground=Ground)  
        return self

    def getSigmaRho(self):
        sigma = self.modelTE.getConductivity()
        rho = self.modelTM.getResistivity()
        return sigma, rho

    def getExHxZs(self):           # used in getZDataBlob.py  not in inversion
        #print("getExHxZs")
        # returns electric and magnetic fields and Z values at sensors for each frequency
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

    def getZsMPI(self, isrc, periods):
        if isrc == 0:
            ZComp = np.empty((self.numP, self.numS), dtype=complex)
            for pi, p in enumerate(periods):
                om = 2.*np.pi/p 
                Z = self.modelTE.getZxy(om=om)
                ZComp[pi] = np.array(self.loc(Z))
        elif isrc == 1:
            ZComp = np.empty((self.numP, self.numS), dtype=complex)
            for pi, p in enumerate(periods):
                om = 2.*np.pi/p 
                Z = self.modelTM.getZyx(om=om)
                ZComp[pi] = np.array(self.loc(Z))
        else: # no MPI case
            ZComp = np.empty((2, self.numP, self.numS), dtype=complex) # 2 because TM and TE
            i = 0
            for pi, p in enumerate(periods):
                om = 2.*np.pi/p 
                Z = self.modelTE.getZxy(om=om)
                ZComp[i, pi] = np.array(self.loc(Z))
            i = 1
            for pi, p in enumerate(periods):
                om = 2.*np.pi/p 
                Z = self.modelTM.getZyx(om=om)
                ZComp[i, pi] = np.array(self.loc(Z))
        return ZComp

