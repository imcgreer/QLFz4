#!/usr/bin/env python

import os
import numpy as np
from scipy.interpolate import LSQUnivariateSpline
from astropy.stats import sigma_clip

mjd_bins = [51500,52000,52700,53200,53500,53850,54200,54600]

def make_lcstruct(mjds,fluxes,ivars,mask=None):
	'''make a lightcurve from a single time series'''
	if mask is None:
		mask = np.zeros(len(mjds),dtype=bool)
	dtype = dict(names=('mjd','flux','ivar','mask'),
	             formats=('f8','f8','f8','b1'))
	lc = np.array(zip(mjds,fluxes,ivars,mask),dtype=dtype)
	return lc

def bin_lc(lc,weights='ivar',thresh=3.0,niter=2):
	'''also does outlier masking'''
	bins = np.digitize(lc['mjd'],mjd_bins)
	nbins = len(mjd_bins)-1
	blc = np.zeros(nbins,dtype=lc.dtype)
	for i in range(nbins):
		binpts = np.where((bins==i+1) & ~lc['mask'])[0]
		if len(binpts) == 0:
			blc['mask'][i] = True
			continue
		f = sigma_clip(lc['flux'][binpts],sig=thresh,iters=niter)
		lc['mask'][binpts] |= f.mask
		blc['mjd'][i] = np.average(lc['mjd'][binpts[~f.mask]])
		flux,ivar = np.ma.average(f,weights=lc['ivar'][binpts],returned=True)
		blc['flux'][i] = flux
		blc['ivar'][i] = ivar
	return lc,blc

def combine_lcstructs(lcs,medianfits):
	'''combine several lightcurves into one (e.g., multiple bands)
	   and sort by MJD'''
	lcs_rel = []
	for lc,mfit in zip(lcs,medianfits):
		lc_rel = lc.copy()
		lc_rel['flux'] -= mfit['median']
		lcs_rel.append(lc_rel)
	lc_rel = np.concatenate(lcs_rel)
	ii = lc_rel['mjd'].argsort()
	return lc_rel[ii]

def fit_median(lc):
	ii = np.where(~lc['mask'])[0]
	median_flux = np.median(lc['flux'][ii])
	chi2 = np.sum((lc['flux'][ii] - median_flux)**2 * lc['ivar'][ii])
	return dict(median=median_flux,chi2=chi2,ndof=len(ii))

def fit_spline(lc,niter=1,reject_thresh=5.0):
	#knots = np.array([52220.,52560.,52940.,53300.,53670.,54030,54400.])
	# knots in the middle of the annual samplings
	#knots = np.array([52240.,52940.,53670.,54420.])
	# knots between the annual samplings
	#knots = np.array([52000.,52740.,53160.,53830.,54220.])
	knots = np.array([52740.,53160.,53830.,54220.])
	iternum = 1
	ii = np.where(~lc['mask'])[0]
	while True:
		knots = knots[(knots>lc['mjd'][ii[0]]) & (knots<lc['mjd'][ii[-1]])]
		spfit = LSQUnivariateSpline(lc['mjd'][ii],lc['flux'][ii],knots,
		                            np.sqrt(lc['ivar'][ii]),k=2)
		chi2v = (lc['flux'][ii]-spfit(lc['mjd'][ii])) * np.sqrt(lc['ivar'][ii])
		if iternum < niter:
			ii = ii[np.abs(chi2v) < reject_thresh]
			iternum += 1
		else:
			break
	chi2 = np.sum(chi2v**2)
	ndof = len(ii)
	return dict(spfit=spfit,fitvals=spfit(lc['mjd']),
	            chi2=chi2,ndof=ndof,good=ii)

