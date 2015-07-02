#!/usr/bin/env python

import copy
import itertools
import numpy as np
from numpy.ma import mrecords

def arr_between(a,b):
	return np.logical_and(a>=b[0],a<=b[1])

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
		lf = np.ma.array(np.ma.zeros(lfShape),
		                 dtype=[('absMag','f4'),
		                        ('counts','f4'),('rawCounts','i4'),
		                        ('countUnc','f4'),('filled','i2'),
		                        ('phi','f8'),('rawPhi','f8'),('sigPhi','f8')],
		                 mask=np.zeros(lfShape,dtype=bool))
		lf = lf.view(mrecords.mrecarray)
		lf['absMag'] = np.repeat(Mbins,len(zbins)).reshape(lfShape)
		# do the counting in bins
		for i in ii:
			lf['rawCounts'][Mi[i]-1,zi[i]-1] += 1
			lf['counts'][Mi[i]-1,zi[i]-1] += self.weights[i]
			lf['countUnc'][Mi[i]-1,zi[i]-1] += self.weights[i]**2
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
			lf['filled'][filled_bin,j] = 1
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
		# construct mask from empty bins
		lf.mask[:] = (lf['rawCounts']==0) | (dVdM == 0)
		# calculate luminosity function from ~ counts/volume
		lf['phi'] = lf['counts'] / (dVdM * self.area_srad)
		lf['rawPhi'] = lf['rawCounts'] / (dVdM * self.area_srad)
		lf['sigPhi'] = np.sqrt(lf['countUnc']) / (dVdM * self.area_srad)
		return lf

