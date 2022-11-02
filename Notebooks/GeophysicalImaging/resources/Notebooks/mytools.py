from esys.escript import *
from esys.escript.linearPDEs import LinearSinglePDE, LinearPDE
from esys.escript.minimizer import CostFunction
from esys.escript.pdetools import ArithmeticTuple
from esys.escript.linearPDEs import SolverOptions
import numpy as np
import logging, time

logger = logging.getLogger('esys')
logger.setLevel(logging.DEBUG)

G=6.67e-11  # m^3/kg/sec^2 gravity constant
mGal=100000.
Mu0=4*np.pi*1e-7

def getGravity(domain, rho, reduced=False):
        model=LinearSinglePDE(domain)
        z=domain.getX()[1]
        model.setValue(A=np.eye(2), q=whereZero(z-sup(z)), Y=-4*np.pi*G*interpolate(rho, Function(domain)))
        u=model.getSolution()
        if reduced:
             gz=-grad(u, ReducedFunction(domain))[domain.getDim()-1]
        else:
            gz=-grad(u)[domain.getDim()-1]
        return gz

class Inversion2DGravity(CostFunction):
    """
    this is simple 2D gravity inversion
    """
    def __init__(self, domain, data, Wd, rho_ref=0, rhodash=1, z_surface=0., mu=1.):
        """
        """
        super(Inversion2DGravity, self).__init__()
        self.domain=domain
        self.data=interpolate(data, Function(domain))
        self.Wd=interpolate(Wd, Function(domain))
        self.rho_ref=rho_ref
        self.rhodash=rhodash
        self.mu=mu
        # PDE for forward model:
        self.model=LinearSinglePDE(self.domain)
        z=domain.getX()[self.domain.getDim()-1]
        self.model.setValue(A=np.eye(self.domain.getDim()), q=whereZero(z-sup(z)))
        # PDE for Hessian inverse:
        self.H=LinearSinglePDE(self.domain)  
        self.H.setValue(A=mu*np.eye(self.domain.getDim()), q=whereNonNegative(z-z_surface))
    def getDualProduct(self, m, r):
        """
        to get <r|m>
        """
        return integrate(r[0]*m + inner(r[1], grad(m)))
    def getNorm(self, m):
        # used to terminate iteration
        return Lsup(m)
    def getArguments(self,m):
        """
        return values that are needed for the cost function evaluation and the gardient (to save compute time)
        """
        rho=self.rhodash*interpolate(m, Function(self.domain))+self.rho_ref

        self.model.setValue(Y=-4*np.pi*G*rho, X=Data())
        u=self.model.getSolution()
        gz=-grad(u)[self.domain.getDim()-1]
        return gz,
    def getValue(self, m, *args):
        """
        returns the value of the cost function
        """
        # misfit:
        if len(args)==0:
            args=self.getArguments(m)
        gz=args[0]
        phi_d=0.5*integrate(self.Wd*(gz-self.data)**2)
        # regularization
        phi_r=0.5*self.mu*integrate(length(grad(m))**2)
        logger.info("phi_d, phi_r = %e, %e"%(phi_d,phi_r))
        return phi_d+phi_r
    
    def getGradient(self, m, *args):
        """
        returns the gradient in the form (X,Y)
        """
        # adjoint problem:
        if len(args)==0:
            args=self.getArguments(m)

        gz=args[0]
        numDim=self.domain.getDim()
        self.model.setValue(X=self.Wd*(self.data-gz)*np.eye(numDim)[numDim-1], Y=Data())
        ustar=self.model.getSolution()
        X=self.mu*grad(m)
        Y=-4*np.pi*G*self.rhodash*interpolate(ustar, X.getFunctionSpace())
        return ArithmeticTuple(Y, X)
    
    def getInverseHessianApproximation(self, m, r, *args, initializeHessian=False):
        self.H.setValue( Y=r[0], X=r[1])
        return self.H.getSolution()

def getMTTE(domain, f, sigma, reduced=False):
        omega=2*np.pi*f
        model=LinearSinglePDE(domain, isComplex=True)

        model.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
        z=domain.getX()[1]
        model.setValue(A=np.eye(2), D=1j*Mu0*omega*interpolate(sigma, Function(domain)), q=whereZero(z-sup(z)), r=1)
        u=model.getSolution()
        if reduced:
             g1=grad(u, ReducedFunction(domain))[1]
        else:
             g1=grad(u)[1]
        Z=1j*omega*Mu0*u/g1
        return Z

class InversionMTTE(CostFunction):
    """
    this is simple 2D gravity inversion
    """
    def __init__(self, domain, f, data, Wd, sigma_ref=1, z_surface=0., mu=1.):
        """
        """
        super(InversionMTTE, self).__init__()
        self.domain=domain
        self.omega=2*np.pi*f
        self.data=interpolate(data, Function(domain))
        self.Wd=interpolate(Wd, Function(domain))
        self.sigma_ref=interpolate(sigma_ref,  Function(domain))
        self.mu=mu
        # PDE for forward model:
        self.model=LinearSinglePDE(self.domain, isComplex=True)
        self.model.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
        z=domain.getX()[self.domain.getDim()-1]
        self.model.setValue(A=np.eye(self.domain.getDim()), q=whereZero(z-sup(z)))
        # PDE for Hessian inverse:
        self.H=LinearSinglePDE(self.domain)  
        self.H.setValue(A=mu*np.eye(self.domain.getDim()), q=whereNonNegative(z-z_surface))
    def getDualProduct(self, m, r):
        """
        to get <r|m>
        """
        return integrate(r[0]*m + inner(r[1], grad(m)))
    def getNorm(self, m):
        # used to terminate iteration
        return Lsup(m)
    def getArguments(self,m):
        """
        return values that are needed for the cost function evaluation and the gradient (to save compute time)
        """
        im=interpolate(m, Function(self.domain))
        sigma=self.sigma_ref*exp(clip(im, minval=-15, maxval=15))
        self.model.setValue(D=1j*Mu0*self.omega*sigma, r=1, X=Data(), Y=Data())
        u=self.model.getSolution()
        g1=grad(u)[1]
        Z=1j*self.omega*Mu0*u/g1
        return sigma, Z, u
    def getValue(self, m, *args):
        """
        returns the value of the cost function
        """
        # misfit:
        if len(args)==0:
            args=self.getArguments(m)
        sigma, Z, u=args
        phi_d=0.5*integrate(self.Wd*abs(Z-self.data)**2)
        # regularization
        phi_r=0.5*self.mu*integrate(length(grad(m))**2)
        logger.info("phi_d, phi_r, sigma= %e, %e, %s"%(phi_d,phi_r, str(sigma)))
        return phi_d+phi_r
    
    def getGradient(self, m, *args):
        """
        returns the gradient in the form (X,Y)
        """
        # adjoint problem:
        if len(args)==0:
            args=self.getArguments(m)

        sigma, Z, u=args
        dsigma=sigma
        g1=grad(u)[1]
        ui=interpolate(u, Function(self.domain))
        self.model.setValue(D=1j*Mu0*self.omega*sigma, r=Data())
        A=1j*self.omega*Mu0*(Z-self.data).conjugate()*self.Wd
        self.model.setValue(X=-A*ui/g1**2*np.eye(2)[1], Y=A/g1)
        ustar=self.model.getSolution()
        X=self.mu*grad(m)
        Y=-(1j*Mu0*self.omega*ui*ustar*dsigma).real()
        return ArithmeticTuple(Y, X)
    
    def getInverseHessianApproximation(self, m, r, *args, initializeHessian=False):
        self.H.setValue( Y=r[0], X=r[1])
        return self.H.getSolution()

class InversionJoint2DGravMTTE(CostFunction):
    """
    this is simple joint inversion for 2D gravity inversion anr MT TE (one frequence)
    """
    def __init__(self, domain, z_surface=0., theta=1.,
                 data_g=0., Wd_g=1, rho_ref=0, rhodash=1, mu_g=1.,
                 f=1, data_mt=0j, Wd_mt=1, sigma_ref=1, mu_mt=1. ):
        """
        """
        super(InversionJoint2DGravMTTE, self).__init__()
        self.EPS=1e-13
        self.domain=domain
        V=vol(self.domain)
        self.theta=theta/V
        z=domain.getX()[self.domain.getDim()-1]
        
        # ===============  Gravity part ==============================
        self.data_g=interpolate(data_g, Function(domain))
        self.Wd_g=interpolate(Wd_g, Function(domain))
        self.rho_ref=rho_ref
        self.rhodash=rhodash
        self.mu_g=mu_g/V
        #      PDE for forward model:
        self.modelG=LinearSinglePDE(self.domain)
        self.modelG.setValue(A=np.eye(self.domain.getDim()), q=whereZero(z-sup(z)))       
        # =============== MT part ===============
        self.omega=2*np.pi*f
        self.data_mt=interpolate(data_mt, Function(domain))
        self.Wd_mt=interpolate(Wd_mt, Function(domain))
        self.sigma_ref=interpolate(sigma_ref,  Function(domain))
        self.mu_mt=mu_mt/V
        #    PDE for forward model:
        self.modelMT=LinearSinglePDE(self.domain, isComplex=True)
        self.modelMT.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
        self.modelMT.setValue(A=np.eye(self.domain.getDim()), q=whereZero(z-sup(z)))
        # =========== PDE for Hessian inverse ======
        self.Hg=LinearSinglePDE(self.domain)  
        self.Hg.setValue(A=self.mu_g*np.eye(self.domain.getDim()), q=whereNonNegative(z-z_surface))
        self.Hmt=LinearSinglePDE(self.domain)  
        self.Hmt.setValue(A=self.mu_mt*np.eye(self.domain.getDim()), q=whereNonNegative(z-z_surface))
    def getDualProduct(self, m, r):
        """
        to get <r|m>
        """
        return integrate(inner(r[0],m) + inner(r[1], grad(m)))
    def getNorm(self, m):
        # used to terminate iteration
        return Lsup(m)
    def getArguments(self,m):
        """
        return values that are needed for the cost function evaluation and the gradient (to save compute time)
        """
        im=interpolate(m, Function(self.domain))
        # ===============  Gravity part ==============================
        rho=self.rhodash*im[0]+self.rho_ref
        self.modelG.setValue(Y=-4*np.pi*G*rho, X=Data())
        u=self.modelG.getSolution()
        gz=-grad(u)[self.domain.getDim()-1]
        # =============== MT part ===============
        sigma=self.sigma_ref*exp(clip(im[1], minval=-15, maxval=15))
        self.modelMT.setValue(D=1j*Mu0*self.omega*sigma, r=1, X=Data(), Y=Data())
        u_mt=self.modelMT.getSolution()
        Z=1j*self.omega*Mu0*u_mt/grad(u_mt)[1]
        logger.info("sigma= %s"%(str(sigma),))  
        logger.info("rho= %s"%(str(rho),))  
        
        return gz, sigma, Z, u_mt
    def getValue(self, m, *args):
        """
        returns the value of the cost function
        """
        # misfit:
        if len(args)==0:
            args=self.getArguments(m)
        gz, sigma, Z, u_mt=args
        # ===============  Gravity part ==============================
        phi_d_g=0.5*integrate(self.Wd_g*(gz-self.data_g)**2)
        print("gz",str(gz*wherePositive(self.Wd_g)), str(self.data_g*wherePositive(self.Wd_g)) )
        # regularization
        phi_r_g=0.5*self.mu_g*integrate(length(grad(m[0]))**2)
        logger.info("phi_d_g , phi_r_g  = %e, %e"%(phi_d_g,phi_r_g))

        # =============== MT part ===============
        phi_d_mt=0.5*integrate(self.Wd_mt*abs(Z-self.data_mt)**2)
        # regularization
        phi_r_mt=0.5*self.mu_mt*integrate(length(grad(m[1]))**2)
        logger.info("phi_d_mt, phi_r_mt = %e, %e"%(phi_d_mt,phi_r_mt))

        
        # === Xgrad term ========================
        gm=grad(m)
        gm0=gm[0]
        gm1=gm[1]
        
        lgm0=length(gm0)
        lgm1=length(gm1)
        d01=inner(gm0,gm1)
        lgm0+=whereZero(lgm0, tol=0)*self.EPS
        lgm1+=whereZero(lgm1, tol=0)*self.EPS
        C=((lgm0*lgm1)**2-d01**2)*(1/lgm0**2+1/lgm1**2)
        phi_c=0.5*self.theta*integrate(C)         
        logger.info("phi_c            = %e"%(phi_c))
        
        return phi_d_mt+phi_r_mt+phi_d_g+phi_r_g
    
    def getGradient(self, m, *args):
        """
        returns the gradient in the form (X,Y)
        """
        # adjoint problem:
        if len(args)==0:
            args=self.getArguments(m)
        gz, sigma, Z, u_mt=args

        Y=Data(0., (2,) , Function(self.domain))
        X=Data(0., (2,self.domain.getDim()), Function(self.domain))
        X[0]=self.mu_g*grad(m[0])
        X[1]=self.mu_mt*grad(m[1])

        # ===============  Gravity part ==============================
        numDim=self.domain.getDim()
        self.modelG.setValue(X=self.Wd_g*(self.data_g-gz)*np.eye(numDim)[numDim-1], Y=Data())
        ustar=self.modelG.getSolution()
        Y[0]+=-4*np.pi*G*self.rhodash*interpolate(ustar, X.getFunctionSpace())
        # ===============  MT part ==============================
        dsigma=sigma
        g1=grad(u_mt)[1]
        ui=interpolate(u_mt, Function(self.domain))
        self.modelMT.setValue(D=1j*Mu0*self.omega*sigma, r=Data())
        A=1j*self.omega*Mu0*(Z-self.data_mt).conjugate()*self.Wd_mt
        self.modelMT.setValue(X=-A*ui/g1**2*np.eye(2)[1], Y=A/g1)
        ustar=self.modelMT.getSolution()
        Y[1]+=-(1j*Mu0*self.omega*ui*ustar*dsigma).real()

        # === Xgrad term ========================       
        gm=grad(m)
        gm0=gm[0]
        gm1=gm[1]
        lgm0=length(gm0)
        lgm1=length(gm1)
        d01=inner(gm0,gm1)
            
        lgm0+=whereZero(lgm0, tol=0)*self.EPS
        lgm1+=whereZero(lgm1, tol=0)*self.EPS
        f=d01*(1/lgm0**2+1/lgm1**2)
        X[0,:]+=self.theta*( (1 + d01**2/lgm0**4)*gm0-f*gm1)
        X[1,:]+=self.theta*( (1 + d01**2/lgm1**4)*gm1-f*gm0)
        return ArithmeticTuple(Y, X)
    
    def getInverseHessianApproximation(self, m, r, *args, initializeHessian=False):
        p=Data(0., m.getShape(), m.getFunctionSpace())
        self.Hg.setValue( Y=r[0][0], X=r[1][0])
        p[0]=self.Hg.getSolution()
        self.Hmt.setValue( Y=r[0][1], X=r[1][1])
        p[1]=self.Hmt.getSolution()
        return p
