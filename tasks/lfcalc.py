#!/usr/bin/env python

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from simqso import sqanalysis,lumfun
import qlffit

class mag2Lum(object):
	def __init__(self,simName,simDir,cosmo):
		self.lcorr = sqanalysis.calcKCorrFromGrid(simName+'_Lum',
		                                          simDir,bandNum=3)
		self.fcorr = sqanalysis.calcKCorrFromGrid(simName+'_Flux',
		                                          simDir,bandNum=3)
		self.cosmo = cosmo
	def __call__(self,m,z,inverse=False):
		if inverse:
			return self.lcorr(m,z) + self.cosmo.distmod(z).value
		else:
			return self.fcorr(m,z) + self.cosmo.distmod(z).value

def qlfz4_selection_function(mag,magErr,flux,fluxErr):
	# no selection criteria - just return True for all
	return (mag[:,3] > 0)

import pandas as pd

def loadSelectionGrid(m_arr, z_arr, fraction, binsM, binsZ):
	"""Take luminosity and redshift arrays and return selection
	   function array
	"""
	df = pd.DataFrame()
	df['M'] = m_arr
	df['z'] = z_arr
	df['Mbin'] = pd.cut(df.M, binsM)
	df['zbin'] = pd.cut(df.z, binsZ)
	fractionCutDF = fraction #.reset_index()
	fractionCutDF.columns = ['M', 'z', 'fraction']
	joined = pd.merge(df, fractionCutDF, how='left', right_on=('M', 'z'),
	                  left_on=('Mbin', 'zbin'))['fraction'].values
	return joined.filledna(0).values

class SelectionFunction(object):
	def __init__(self,simName,simDir,photo_complete,spec_complete,m2M,
	             sfx='',**kwargs):
		self.m2M = m2M
		if True:
			self.selcrit = qlfz4_selection_function
			self.M_selfun = sqanalysis.calcSelectionFunctionFromGrid(
			                                            simName+'_Lum',
			                                            self.selcrit,simDir,
			                                            **kwargs)
			self.m_selfun = sqanalysis.calcSelectionFunctionFromGrid(
			                                            simName+'_Flux',
			                                            self.selcrit,simDir,
			                                            **kwargs)
		else:
			self.binsM = np.arange(-27.5, -23.5, 0.2)
			self.binsZ = np.arange(3.25, 4.75, 0.1)
			self.fractionDF = pd.read_csv('../fractionNov30.csv')
			self.selfun = lambda M,z: loadSelectionGrid(M,z,
			                            self.fractionDF,self.binsM,self.binsZ)
		self.photo_complete = photo_complete
		self.spec_complete = spec_complete
	def _photo_complete_call(self,m,z,absMag=False):
		if type(self.photo_complete) is float:
			return self.photo_complete
		elif absMag:
			return self.photo_complete(m-self.m2M(m,z,inverse=True),z)
		else:
			return self.photo_complete(m,z)
	def _spec_complete_call(self,m,z,absMag=False):
		if type(self.spec_complete) is float:
			return self.spec_complete
		elif absMag:
			return self.spec_complete(m-self.m2M(m,z,inverse=True),z)
		else:
			return self.spec_complete(m,z)
	def _color_complete_call(self,m,z,absMag=False):
		if absMag:
			#return self.selfun(m,z)
			return self.M_selfun(m,z)
		else:
			#return self.selfun(m-self.m2M(m,z,inverse=True),z)
			return self.m_selfun(m,z)
	def __call__(self,m,z,absMag=False):
		return self._color_complete_call(m,z,absMag) * \
		       self._photo_complete_call(m,z,absMag) * \
		       self._spec_complete_call(m,z,absMag)

from scipy.integrate import quad

def get_dVdzdO(cosmo):
	_zfun = lumfun.interp_dVdzdO((3.0,5.0),cosmo)
	dVdzdO = lambda z1,z2: quad(_zfun,z1,z2)[0]
	return dVdzdO

def placeholder_survey():
	i_lim = 22.0
	surveyArea = 60.6
	simName = 'SDSS'
	simDir = '../simulations/'
	cosmo = FlatLambdaCDM(H0=70,Om0=0.3)
	m2M = mag2Lum(simName,simDir,cosmo)
	dat = Table.read('../placeholder_north_survey.csv')
	survey = qlffit.QuasarSurvey(dat['median_mag_i']-2*2.086*dat['ebv'],
	                             dat['z'],i_lim,surveyArea,m2M)
	photo_complete = 1.0
	spec_complete = 1.0
	selfun = SelectionFunction(simName,simDir,photo_complete,spec_complete,m2M)
	survey.set_selection_function(selfun)
	Medges = np.arange(-29.,-21.,0.5)
	zedges = np.arange(3.25,4.51,0.55)
	lf = survey.calcBinnedLF(Medges,zedges,get_dVdzdO(cosmo))
	return lf

