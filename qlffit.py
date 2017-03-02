#!/usr/bin/env python

import copy
import itertools
import numpy as np
from numpy.ma import mrecords
from scipy.integrate import dblquad,romberg
from scipy.ndimage.filters import convolve
from scipy import optimize
from astropy.table import Table
from simqso.lumfun import interp_dVdzdO

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
		dVdzdO = interp_dVdzdO((zedges[0]-0.1,zedges[1]+0.1),self.m2M.cosmo)
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

def _integrate_fast(integrand,zrange,Mrange,integrate_kwargs):
	integrate_kwargs.setdefault('tol',integrate_kwargs.get('epsabs',1e-3))
	integrate_kwargs.setdefault('rtol',integrate_kwargs.get('epsrel',1e-3))
	inner = lambda z: romberg(integrand,*Mrange,args=(z,),**integrate_kwargs)
	outer = romberg(inner,*zrange,**integrate_kwargs)
	return outer

def _integrate_full(integrand,zrange,Mrange,nM,nz,integrate_kwargs):
	integrate_kwargs.setdefault('epsabs',1e-3)
	integrate_kwargs.setdefault('epsrel',1e-3)
	# do a full double integration, but break it up into chunks
	zz = np.linspace(zrange[0],zrange[1],nz)
	MM = np.linspace(Mrange[0],Mrange[1],nM)
	_sum = 0
	for z1,z2 in zip(zz[:-1],zz[1:]):
		for M1,M2 in zip(MM[:-1],MM[1:]):
			intp,err = dblquad(integrand, z1, z2,
			                   lambda z: M1,lambda z: M2,
			                   **integrate_kwargs)
			_sum += intp
	return _sum

def joint_qlf_likelihood_fun(par,surveys,Phi_Mz,dV_dzdO,zrange,Mrange,
                             nM,nz,fast,integrate_kwargs):
	first_term,second_term = 0.0,0.0
	for s in surveys:
		# first term: sum over each observed quasar
		# this used to work...
		#first_term += -2*np.sum(np.log(Phi_Mz(s.M,s.z,par)))
		prod = [ p_Mz*Phi_Mz(M,z,par) 
		            for M,z,p_Mz in zip(s.M,s.z,s.weights**-1) ]
		prod = np.clip(prod,1e-20,np.inf)
		first_term += -2*np.sum(np.log(prod))
		# second term: integral of LF over available volume
		integrand = lambda M,z: Phi_Mz(M,z,par) * s.p_Mz(M,z) * dV_dzdO(z)
		if fast:
			_sum = _integrate_fast(integrand,zrange,Mrange,integrate_kwargs)
		else:
			_sum = _integrate_full(integrand,zrange,Mrange,
			                       nM,nz,integrate_kwargs)
		second_term += 2 * s.area_srad * _sum
	#
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
	def __init__(self):
		self.routine = optimize.fmin
		self.args = ()
		self.kwargs = {'full_output':True,'xtol':1e-3,'ftol':1e-3}

def joint_qlf_mle(surveys,Phi_Mz,ival,dVdzdO,zrange,Mrange,
                  minimizer=None,nM=20,nz=10,fast=True,integrate_kwargs={}):
	'''
	Maximum likelihood estimation of QLF parameters, using data points
	from multiple surveys.
	-
	surveys: list of QuasarSurvey objects, contains data points and
	         survey parameters for each survey
	 Phi_Mz: function Phi(M,z) that returns comoving number density at M,z
       ival: vector of initial values for parameters in Phi
	 dVdzdO: dV/dz*dOmega 
	'''
	if minimizer is None:
		minimizer = NelderMeadFit()
	scale = Phi_Mz.scale
	Phi_Mz.set_scale('linear')
	fit = minimizer(joint_qlf_likelihood_fun,ival,*minimizer.args,
	                args=(surveys,Phi_Mz,dVdzdO,zrange,Mrange,
	                      nM,nz,fast,integrate_kwargs),
	                **minimizer.kwargs)
	Phi_Mz.set_scale(scale)
	return fit

