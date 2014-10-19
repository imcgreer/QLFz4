#!/usr/bin/env python

import os
import numpy as np
from scipy.interpolate import LSQUnivariateSpline
from astropy.stats import sigma_clip

mjd_bins = [51500,52000,52400,52750,53200,53500,53850,54200,54600]

def make_lcstruct(mjds,fluxes,ivars,mask=None):
	'''make a lightcurve from a single time series'''
	if mask is None:
		mask = np.zeros(len(mjds),dtype=bool)
	dtype = dict(names=('mjd','flux','ivar','mask'),
	             formats=('f8','f8','f8','b1'))
	lc = np.array(zip(mjds,fluxes,ivars,mask),dtype=dtype)
	return lc

def bin_lc(lc,weights='equal',invvar=True,thresh=2.2,niter=2):
	'''generate annually-averaged lightcurves
	   weights: weighting scheme to apply when averaging points
	       'equal' [default] points are given equal weight
	       'midpt' points are weighted in inverse proportion to their
	               distance from the midpoint of the observing season.
	               i.e., Aug/Dec have lowest weights, Oct highest
	   ivar: include inverse variance in weights
	   also does outlier masking using sigma_clip
	'''
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
		if weights == 'equal':
			w = 1.0
		elif weights == 'midpt':
			delt = lc['mjd'][binpts] - blc['mjd'][i]
			tw = lc['mjd'][binpts[-1]] - lc['mjd'][binpts[0]]
			w = 1 - delt/tw # =1 at midpoint, ~0.5 at edges
		if invvar:
			w *= lc['ivar'][binpts]
		flux,ivar = np.ma.average(f,weights=w,returned=True)
		blc['flux'][i] = flux
		blc['ivar'][i] = ivar
	return lc,blc

def combine_lcstructs(lcs,medianfits):
	'''combine several lightcurves into one (e.g., multiple bands)
	   and sort by MJD'''
	lcs_rel = []
	for lc,mfit in zip(lcs,medianfits):
		lc_rel = lc.copy()
		lc_rel['flux'] /= mfit['median']
		lc_rel['ivar'] *= mfit['median']**2
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
	if lc['mjd'][0] < 52000:
		knots = np.array([52500.,53500.,54000.])
	else:
		knots = np.array([52700.,53500.,54000.])
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

def load_lsst_stripe82_lightcurves(lcdir):
	procf = os.path.join(lcdir,'process_i.dat')
	objs = np.genfromtxt(procf,
	                     usecols=(1,2,3,5),names='ra,dec,z,id',
	                     dtype=('f8','f8','f4','i8'),
	                     converters={5: lambda s: long(s[3:-5])})
	lightcurves,annualcurves = {},{}
	mfits = {}
	for objid in objs['id']:
		lc,blc = {},{}
		mfit = {}
		for b in 'gri':
			lcdatfile = os.path.join(lcdir,b,'%d.dat'%objid)
			lcdat = np.loadtxt(lcdatfile,unpack=True)
			if True:
				# fluxes are ~1e-29, scale them
				lcdat[1:] *= 1e29
			lc[b] = make_lcstruct(*lcdat)
			lc[b],blc[b] = bin_lc(lc[b])
			mfit[b] = fit_median(lc[b])
		lc['gri'] = combine_lcstructs([lc[b] for b in 'gri'],
		                              [mfit[b] for b in 'gri'])
		lightcurves[objid] = lc
		annualcurves[objid] = blc
		mfits[objid] = mfit
	return dict(objs=objs,lcs=lightcurves,annual=annualcurves,medianfit=mfits)

