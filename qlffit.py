#!/usr/bin/env python

import copy
import itertools
import numpy as np
from numpy.ma import mrecords
from scipy.integrate import quad,dblquad,romberg,simps
from scipy.ndimage.filters import convolve
from scipy import optimize
from astropy.table import Table
from simqso.lumfun import interp_dVdzdO
from simqso import lumfun

def arr_between(a,b):
	return np.logical_and(a>=b[0],a<b[1])

class QuasarSurvey(object):
	def __init__(self,m,z,m_lim,skyArea,m2M):
		'''__init__(self,m,z,m_lim,skyArea,m2M):
		    Inputs are survey description.
		    m:       apparent magnitudes of objects
		    z:       redshifts
		    m_lim:   limiting apparent magnitude of survey
		    skyArea: area of survey in deg^2
		    m2M:     function taking apparent mag and redshift as arguments,
	                 along with a keyword ``inverse'', and returning the
		             conversion from apparent mag to absolute mag, or
		             the reverse if inverse=True.
		             Must include both cosmological and k-corrections.
		             i.e., M = m - m2M(m,z) = m - DM(z) - K(m,z)
		              and  m = M + m2M(M,z,inverse=True)
		            Allows for luminosity-dependent k-corrections.
		'''
		self.m = m
		self.z = z
		self.m_lim = m_lim
		self.skyArea = skyArea
		self.skyFraction = skyArea/41252.96
		self.area_srad =  self.skyFraction * 4*np.pi
		self.N = len(m)
		# convert apparent to absolute magnitudes
		self.m2M = m2M
		self.m2M_val = m2M(m,z)
		self.M = m - self.m2M_val
	def set_selection_function(self,selfun):
		self.selfun = selfun
		self.p_Mz = lambda M,z: self.selfun(M,z,absMag=True)
		self.weights = np.clip(self.selfun(self.m,self.z),1e-20,1.0)**-1
	def Nofz(self,zedges):
		N = np.empty(zedges.shape[0]-1)
		Ncorr = np.empty_like(N)
		for i,z1,z2 in zip(itertools.count(),zedges[:-1],zedges[1:]):
			ii = np.where(arr_between(self.z,(z1,z2)))[0]
			N[i] = len(ii)
			Ncorr[i] = np.sum(self.weights[ii])
		return N,Ncorr
	def take(self,ii):
		rv = copy.copy(self)
		for k in ['m','z','M','weights']:
			rv.__dict__[k] = rv.__dict__[k][ii]
		return rv
	def __getitem__(self,index):
		if type(index) is np.ndarray:
			return self.take(index)
		else:
			return (self.m[index],self.z[index],self.M[index])
	@staticmethod
	def init_lf_table(Mbins,zbins):
		lfShape = Mbins.shape + zbins.shape
		lfTab = Table(masked=True)
		lfTab['absMag'] = Mbins.astype(np.float32)
		lfTab['counts'] = np.zeros(lfShape,dtype=np.float32)
		lfTab['rawCounts'] = np.zeros(lfShape,dtype=np.int32)
		lfTab['countUnc'] = np.zeros(lfShape,dtype=np.float32)
		lfTab['filled'] = np.zeros(lfShape,dtype=np.bool)
		lfTab['phi'] = np.zeros(lfShape,dtype=np.float64)
		lfTab['rawPhi'] = np.zeros(lfShape,dtype=np.float64)
		lfTab['sigPhi'] = np.zeros(lfShape,dtype=np.float64)
		return lfTab
	def getinbounds(self,Medges,zedges):
		# identify which bins are within the flux limit by converting the
		# the luminosity bins to fluxes
		Mbounds,zbounds = np.meshgrid(Medges,zedges,indexing='ij')
		mbounds = Mbounds + self.m2M(Mbounds,zbounds,inverse=True)
		inbounds = mbounds < self.m_lim
		# this sums the bin edges 2x2: 
		#   4=full covg, 0=no covg, otherwise partial
		inbounds = convolve(inbounds.astype(int),np.ones((2,2)))[:-1,:-1]
		return inbounds
	def calcBinnedLF(self,Medges,zedges,**kwargs):
		'''calcBinnedLF(self,Medges,zedges,**kwargs)
		    calculate binned luminosity function from the stored survey data.
		    Medges: array defining bin edges in absolute mag
		    zedges: array defining bin edges in redshift
		'''
		# kind of hacky to access cosmo through m2M...
		dVdzdO = interp_dVdzdO(zedges,self.m2M.cosmo)
		#
		Mbins = Medges[:-1] + np.diff(Medges)/2
		zbins = zedges[:-1] + np.diff(zedges)/2
		lfShape = Mbins.shape + zbins.shape
		# assign data points to bins and trim out-of-bounds objects
		Mi = np.digitize(self.M,Medges) - 1
		zi = np.digitize(self.z,zedges) - 1
		ii = np.where( (Mi>=0) & (Mi<len(Mbins)) &
		               (zi>=0) & (zi<len(zbins)) )[0]
		# do the counting in bins
		lf = self.init_lf_table(Mbins,zbins)
		np.add.at( lf['rawCounts'], (Mi[ii],zi[ii]),                1   )
		np.add.at(    lf['counts'], (Mi[ii],zi[ii]), self.weights[ii]   )
		np.add.at(  lf['countUnc'], (Mi[ii],zi[ii]), self.weights[ii]**2)
		#
		inbounds = self.getinbounds(Medges,zedges)
		lf['filled'][:] = (inbounds==4)
		# calculate bin volumes by integrating dVdM = (dV/dz)dzdM
		#   ... note if there were many redshift bins, could save time
		#       by only calculating dV once for each filled bin within
		#       each redshift slice
		binVol = np.zeros(lfShape)
		for i,j in zip(*np.where(inbounds > 0)):
			Mlim = lambda z: np.clip(self.m_lim-self.m2M(self.m_lim,z),
			                         Medges[i],Medges[i+1])
			binVol[i,j],_ = dblquad(lambda M,z: dVdzdO(z),
			                        zedges[j],zedges[j+1],
			                        lambda z: Medges[i],Mlim)
		# calculate luminosity function from ~ counts/volume
		mask = (lf['rawCounts']==0) | (binVol == 0)
		binVol = np.ma.array(binVol * self.area_srad, mask=mask)
		lf['phi'] = np.ma.divide(lf['counts'],binVol)
		lf['rawPhi'] = np.ma.divide(lf['rawCounts'],binVol)
		lf['sigPhi'] = np.ma.divide(np.ma.sqrt(lf['countUnc']),binVol)
		return lf

class QLFIntegrator(object):
	def __init__(self,Mrange,zrange,dVdzdO):
		self.Mrange = Mrange
		self.zrange = zrange
		self.dVdzdO = dVdzdO
		self.int_kwargs = {}

class FullQLFIntegrator(QLFIntegrator):
	def __init__(self,Mrange,zrange,dVdzdO,**kwargs):
		super(FullQLFIntegrator,self).__init__(Mrange,zrange,dVdzdO)
		self.nM = kwargs.pop('nM',20)
		self.nz = kwargs.pop('nz',10)
		self.int_kwargs.setdefault('epsabs',kwargs.pop('epsabs',1e-3))
		self.int_kwargs.setdefault('epsrel',kwargs.pop('epsrel',1e-3))
		self.zz = np.linspace(self.zrange[0],self.zrange[1],self.nz)
		self.MM = np.linspace(self.Mrange[0],self.Mrange[1],self.nM)
	def __call__(self,Phi_Mz,p_Mz,par):
		#
		integrand = lambda M,z: Phi_Mz(M,z,par) * p_Mz(M,z) * self.dVdzdO(z)
		lfsum = 0
		for z1,z2 in zip(self.zz[:-1],self.zz[1:]):
			for M1,M2 in zip(self.MM[:-1],self.MM[1:]):
				intp,err = dblquad(integrand, z1, z2,
				                   lambda z: M1,lambda z: M2,
				                   **self.int_kwargs)
				lfsum += intp
		return lfsum

class FastQLFIntegrator(QLFIntegrator):
	def __init__(self,Mrange,zrange,dVdzdO,**kwargs):
		super(FastQLFIntegrator,self).__init__(Mrange,zrange,dVdzdO)
		self.int_kwargs.setdefault('divmax',kwargs.pop('divmax',20))
		self.int_kwargs.setdefault('tol',kwargs.pop('epsabs',1e-3))
		self.int_kwargs.setdefault('rtol',kwargs.pop('epsrel',1e-3))
	def __call__(self,Phi_Mz,p_Mz,par):
		#
		integrand = lambda M,z: Phi_Mz(M,z,par) * p_Mz(M,z) * self.dVdzdO(z)
		inner = lambda z: romberg(integrand,*self.Mrange,args=(z,),
		                          **self.int_kwargs)
		outer = romberg(inner,*self.zrange,**self.int_kwargs)
		return outer

class FasterQLFIntegrator(QLFIntegrator):
	def __init__(self,Mrange,zrange,dVdzdO,**kwargs):
		super(FasterQLFIntegrator,self).__init__(Mrange,zrange,dVdzdO)
		self.minProb = kwargs.pop('minProb',1e-3)
		in_MBinW = kwargs.pop('MBinWidth',0.1)
		in_zBinW = kwargs.pop('zBinWidth',0.05)
		self.nM = int( np.diff(self.Mrange) / in_MBinW ) + 1 
		self.nz = int( np.diff(self.zrange) / in_zBinW ) + 1
		#
		self.Medges = np.linspace(self.Mrange[0],self.Mrange[1],self.nM)
		self.zedges = np.linspace(self.zrange[0],self.zrange[1],self.nz)
		self.MBinW = np.diff(self.Medges)[0]
		self.zBinW = np.diff(self.zedges)[0]
		#
		self.dV = self.dVdzdO(self.zedges)
		self.Mi,self.zi = np.meshgrid(self.Medges,self.zedges,indexing='ij')
		#
		self.p_Mz_cache = {}
		self.lowProbMask = {}
	def _get_p_Mz_grid(self,p_Mz):
		p_Mz_grid = self.p_Mz_cache.get(p_Mz)
		if p_Mz_grid is None:
			p_Mz_grid = p_Mz(self.Mi,self.zi)
			self.p_Mz_cache[p_Mz] = p_Mz_grid
			self.lowProbMask[p_Mz] = p_Mz_grid > self.minProb
		return p_Mz_grid,self.lowProbMask[p_Mz]
	def __call__(self,Phi_Mz,p_Mz,par):
		#
		p_Mz_grid,mask = self._get_p_Mz_grid(p_Mz)
		Phi_Mz_grid = Phi_Mz(self.Mi,self.zi,par)
		#
		lfsum_z = simps(Phi_Mz_grid * p_Mz_grid * self.dV, dx=self.zBinW)
		lfsum = simps(lfsum_z, dx=self.MBinW)
		return lfsum

def joint_qlf_likelihood_fun(par,surveys,lfintegrator,Phi_Mz,verbose):
	min_prob = 1e-3
	first_term,second_term = 0.0,0.0
	for s in surveys:
		# first term: sum over each observed quasar
		p_Mizi = s.weights**-1
		ii = np.where(p_Mizi > min_prob)[0]
		prod = p_Mizi[ii] * Phi_Mz(s.M[ii],s.z[ii],par) 
		first_term += -2*np.sum(np.log(prod))
		# second term: integral of LF over available volume
		lfsum = lfintegrator(Phi_Mz,s.p_Mz,par)
		second_term += 2 * s.area_srad * lfsum
	if verbose:
		print 'testing ',par,first_term,second_term
	return first_term + second_term

class FitMethod(object):
	def __init__(self):
		pass
	def __call__(self,*args,**kwargs):
		return self.routine(*args,**kwargs)
	def set_bounds(self,exclude_list=[]):
		pass

class NelderMeadFit(FitMethod):
	def __init__(self,verbose=False):
		self.routine = optimize.fmin
		self.args = ()
		self.kwargs = {'full_output':True,'xtol':1e-3,'ftol':1e-3,
		               'disp':verbose}

class JointQLFFitter(object):
	def __init__(self,Mrange,zrange,cosmo,qlfModel,**kwargs):
		self.likefun = joint_qlf_likelihood_fun
		self.Mrange = Mrange
		self.zrange = zrange
		self.dVdzdO = lumfun.interp_dVdzdO(zrange,cosmo)
		self.qlfModel = qlfModel
		self.qlfModel.set_scale('linear')
		self.fitMethod = kwargs.get('fit_method',NelderMeadFit())
		self.set_integrate_mode(kwargs.get('integrate_mode','fast'),
		                        kwargs.get('integrate_kwargs',{}))
		self.verbose = kwargs.get('verbose',False)
	def set_integrate_mode(self,mode,integrate_kwargs={}):
		self.integrate_mode = mode
		self.integrate_kwargs = integrate_kwargs
		if self.integrate_mode == 'full':
			self.lfintegrator = FullQLFIntegrator(self.Mrange,self.zrange,
			                                      self.dVdzdO,
			                                      **self.integrate_kwargs)
		elif self.integrate_mode == 'fast':
			self.lfintegrator = FastQLFIntegrator(self.Mrange,self.zrange,
			                                      self.dVdzdO,
			                                      **self.integrate_kwargs)
		elif self.integrate_mode == 'reallyfast':
			self.lfintegrator = FasterQLFIntegrator(self.Mrange,self.zrange,
			                                        self.dVdzdO,
			                                        **self.integrate_kwargs)
		else:
			raise ValueError
	def fit(self,surveys,qlfModel=None,initVals=None):
		if qlfModel is None:
			qlfModel = self.qlfModel
		if initVals is None:
			initVals = list(qlfModel.getpar())
		likefunArgs = (surveys,self.lfintegrator,qlfModel,self.verbose)
		res = self.fitMethod(self.likefun,initVals,*self.fitMethod.args,
		                     args=likefunArgs,**self.fitMethod.kwargs)
		self.lastFit = res
		return res
	def getModel(self):
		rv = self.qlfModel.copy()
		rv.setpar(self.lastFit[0])
		return rv
	def getS(self,surveys,qlfModel=None,par=None):
		if qlfModel is None:
			qlfModel = self.qlfModel
		if par is None:
			par = qlfModel.getpar()
		likefunArgs = (surveys,self.lfintegrator,qlfModel,self.verbose)
		return self.likefun(par,*likefunArgs)
	def varyFitParam(self,paramName,surveys,ntry=None,logRange=None):
		if ntry is None:
			ntry = 50
		# XXX all of this is not right if these params have more than one
		#     free value
		if logRange is None:
			logrange = {
			  'logPhiStar':(-1.5,0.0), 'MStar':(-1.5,0.0),
			  'alpha':(-2.0,0.3), 'beta':(-2.0,0.3),
			}[paramName]
		logbins = logrange + (ntry,)
		#
		S0 = self.getS(surveys)
		print 'S0 is ',S0,' at ',self.qlfModel.params[paramName].get()
		rv = {}
		#
		for i,pval0 in self.qlfModel.params[paramName].iterfree():
			fitvals = [(pval0,S0)]
			qlfModel = self.qlfModel.copy()
			print 'trying %s[#%d]' % (paramName,i)
			for sgn in [-1,1]:
				delv = sgn*np.logspace(*logbins)
				for dv in delv:
					qlfModel.params[paramName].set(pval0+dv,i=i)
					qlfModel.params[paramName].fix(i)
					S = self.fit(surveys,qlfModel=qlfModel)[1]
					qlfModel.params[paramName].free(i)
					if sgn < 0:
						fitvals.insert(0, (pval0+dv, S) )
					else:
						fitvals.append(   (pval0+dv, S) )
					print ' '.join(['%.3f']*6) % (pval0,S0,pval0+dv,S,dv,S-S0)
					if S-S0 > 10:
						# this is more than 3 sigma
						break
			rv[i] = np.array(fitvals)
		return rv
	def sampleModels(self,sigParam,surveys,n=100):
		S0 = self.getS(surveys)
		par0 = self.qlfModel.getpar()
		qlfModel = self.qlfModel.copy()
		S = np.zeros(n)
		allpar = np.zeros((n,len(par0)))
		for i in range(n):
			par = par0 + sigParam*np.random.normal(size=len(sigParam))
			S[i] = self.getS(surveys,qlfModel,par)
			allpar[i] = par
		return Table(dict(par=allpar,dS=(S-S0)))

