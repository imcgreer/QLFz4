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
	def set_color_selfun(self,selfun):
		self.color_selfun = selfun
		self.color_selval = selfun.value(self.M,self.z)
	def set_photo_complete(self,photo_complete):
		# have to translate from absolute to apparent mag 
		self.photo_complete = \
		     lambda M,z: photo_complete(M + self.m2M(M,z,inverse=True))
		# for the internal weights, use the known observed mag
		self.photo_complete_val = photo_complete(self.m)
	def set_spec_complete(self,spec_complete):
		self.spec_complete = \
		     lambda M,z: spec_complete(M + self.m2M(M,z,inverse=True),z)
		self.spec_complete_val = spec_complete(self.m,self.z)
	def calc_selection_function(self):
		self.p_Mz = lambda M,z: np.clip(self.color_selfun(M,z) * 
		                                self.photo_complete(M,z) * 
		                                self.spec_complete(M,z), 0, 1)
		self.weights = 1 / np.clip(self.color_selval *
		                           self.photo_complete_val *
		                           self.spec_complete_val, 0, 1)
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
		for k in ['m','z','M','Kcorrval',
		          'color_selval','photo_complete_val','spec_complete_val',
		          'weights']:
			rv.__dict__[k] = rv.__dict__[k][ii]
		return rv
	def __getitem__(self,index):
		if type(index) is np.ndarray:
			return self.take(index)
		else:
			return (self.m[index],self.z[index],self.M[index])
	def calcLF(self,Medges,zedges,dVdzdO,**kwargs):
		Nsubz = kwargs.get('Nsubz',20)
		NsubM = kwargs.get('NsubM',20)
		ii = np.where(arr_between(self.M,(Medges[0],Medges[-1])) &
		              arr_between(self.z,(zedges[0],zedges[-1])))[0]
		Mbins = Medges[:-1] + np.diff(Medges)/2
		zbins = zedges[:-1] + np.diff(zedges)/2
		Mi = np.digitize(self.M,Medges)
		zi = np.digitize(self.z,zedges)
		# create a structured array to hold the LF bin data
		lfShape = Mbins.shape + zbins.shape
		lf = np.ma.array(np.ma.zeros(lfShape),
		                 dtype=[('counts','f4'),('rawCounts','i4'),
		                        ('countUnc','f4'),('filled','i2'),
		                        ('phi','f8'),('rawPhi','f8'),('sigPhi','f8')],
		                 mask=np.zeros(lfShape,dtype=bool))
		lf = lf.view(mrecords.mrecarray)
		#
		for i in ii:
			lf['rawCounts'][Mi[i]-1,zi[i]-1] += 1
			lf['counts'][Mi[i]-1,zi[i]-1] += self.weights[i]
			lf['countUnc'][Mi[i]-1,zi[i]-1] += self.weights[i]**2
		#
		print '''need to be careful -- bins may also not be filled in dM.
		         need to change the last M bin to have dM = M1-m2M(mlim)'''
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
			dV = dVdzdO(zedges[j],zedges[j+1])
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
							e_dVdM += dVdzdO(e_z[kj],e_z[kj+1]) * \
							           (e_M[ki+1] - e_M[ki])
				dVdM[i,j] = e_dVdM
		# construct mask from empty bins
		lf.mask[:] = (lf['rawCounts']==0) | (dVdM == 0)
		# calculate luminosity function from ~ counts/volume
		lf['phi'] = lf['counts'] / (dVdM * self.area_srad)
		lf['rawPhi'] = lf['rawCounts'] / (dVdM * self.area_srad)
		lf['sigPhi'] = np.sqrt(lf['countUnc']) / (dVdM * self.area_srad)
		return lf

