#!/usr/bin/env python

import copy
import itertools
import numpy as np
from numpy.ma import mrecords
from scipy.integrate import dblquad,romberg
from scipy import optimize
from astropy.table import Table

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
	def calcBinnedLF(self,Medges,zedges,volumefun,**kwargs):
		'''calcBinnedLF(self,Medges,zedges,volumefun,**kwargs)
		    calculate binned luminosity function from the stored survey data.
		    Medges: array defining bin edges in absolute mag
		    zedges: array defining bin edges in redshift
		    volumefun: function which returns volume of bin in Mpc^3 * srad
		               in redshift range z1,z2; called as volumefun(z1,z2) 
		'''
		Nsubz = kwargs.get('Nsubz',20)
		NsubM = kwargs.get('NsubM',20)
		ii = np.where(arr_between(self.M,(Medges[0],Medges[-1])) &
		              arr_between(self.z,(zedges[0],zedges[-1])))[0]
		Mbins = Medges[:-1] + np.diff(Medges)/2
		zbins = zedges[:-1] + np.diff(zedges)/2
		Mi = np.digitize(self.M,Medges)
		zi = np.digitize(self.z,zedges)
		# create a masked structured array to hold the LF bin data
		lfShape = Mbins.shape + zbins.shape
		lf = Table(masked=True)
		lf['absMag'] = Mbins.astype(np.float32)
		lf['counts'] = np.zeros(lfShape,dtype=np.float32)
		lf['rawCounts'] = np.zeros(lfShape,dtype=np.int32)
		lf['countUnc'] = np.zeros(lfShape,dtype=np.float32)
		lf['filled'] = np.zeros(lfShape,dtype=np.bool)
		lf['phi'] = np.zeros(lfShape,dtype=np.float64)
		lf['rawPhi'] = np.zeros(lfShape,dtype=np.float64)
		lf['sigPhi'] = np.zeros(lfShape,dtype=np.float64)
		# do the counting in bins
		np.add.at(lf['rawCounts'],(Mi[ii]-1,zi[ii]-1),1)
		np.add.at(lf['counts'],(Mi[ii]-1,zi[ii]-1),self.weights[ii])
		np.add.at(lf['countUnc'],(Mi[ii]-1,zi[ii]-1),self.weights[ii]**2)
		# ???
		print '''need to be careful -- bins may also not be filled in dM.
		         need to change the last M bin to have dM = M1-m2M(mlim)'''
		# calculate volume densities within the bins, accounting for 
		# unfilled bins
		dVdM = np.zeros(lfShape)
		dM = np.diff(Medges)
		m_min = np.empty(lfShape,dtype=float)
		m_max = np.empty(lfShape,dtype=float)
		for j in range(zbins.shape[0]):
			m_min[:,j] = Medges[:-1] + \
			    self.m2M(Medges[:-1],zedges[j],inverse=True)
			m_max[:,j] = Medges[1:] + \
			    self.m2M(Medges[1:],zedges[j+1],inverse=True)
			filled_bin = np.where((m_min[:,j] < self.m_lim) & 
			                      (m_max[:,j] < self.m_lim))[0]
			partial_bin = np.where((m_min[:,j] < self.m_lim) & 
			                      (m_max[:,j] > self.m_lim))[0]
			lf['filled'][filled_bin,j] = True
			# filled bins
			dV = volumefun(zedges[j],zedges[j+1])
			dVdM[filled_bin,j] = dV * dM[j]
			# partial bins
			for i in partial_bin:
				e_z = np.linspace(zedges[j],zedges[j+1],Nsubz)
				e_M = np.linspace(Medges[i],Medges[i+1],NsubM)
				e_dVdM = 0.0
				for ki in range(NsubM-1):
					for kj in range(Nsubz-1):
						m_min_ = e_M[ki] +  \
						         self.m2M(e_M[ki],e_z[kj],inverse=True) 
						m_max_ = e_M[ki+1] + \
						         self.m2M(e_M[ki+1],e_z[kj+1],inverse=True)
						if (m_min_ < self.m_lim) and (m_max_ < self.m_lim):
							# only do filled bins in subgrid
							e_dVdM += volumefun(e_z[kj],e_z[kj+1]) * \
							           (e_M[ki+1] - e_M[ki])
				dVdM[i,j] = e_dVdM
		# calculate luminosity function from ~ counts/volume
		mask = (lf['rawCounts']==0) | (dVdM == 0)
		vol = np.ma.array(dVdM * self.area_srad, mask=mask)
		lf['phi'] = np.ma.divide(lf['counts'],vol)
		lf['rawPhi'] = np.ma.divide(lf['rawCounts'],vol)
		lf['sigPhi'] = np.ma.divide(np.ma.sqrt(lf['countUnc']),vol)
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

