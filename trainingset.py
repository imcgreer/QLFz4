#!/usr/bin/env python

import os
import numpy as np
from astropy.io import fits
from astrotools.idmstuff import matchlists

import lcfit
import lcplot

qlfz4dir = os.path.join(os.environ['QLFZ4DATA'], '2014October')

def load_training_catalog():
	trainingdataf = os.path.join(qlfz4dir,'TRAININGwgriReduced.fits')
	return fits.getdata(trainingdataf,1)

def load_all():
	tset = {}
	tset['catalog'] = load_training_catalog()
	for c in ['QSO','GALAXY','STAR']:
		lcdir = os.path.join(qlfz4dir,'trainingDataFiles_flux',c)
		tset[c] = lcfit.load_lsst_stripe82_lightcurves(lcdir)
		m1,m2 = matchlists(tset[c]['objs']['id'],
		                   tset['catalog']['deepSourceId_1'])
		assert np.all(m1==np.arange(len(tset[c]['objs']['id'])))
		tset[c+'cat'] = tset['catalog'][m2]
	return tset

def add_spline_fit(tset):
	for c in ['QSO','GALAXY','STAR']:
		sfits = {}
		objids = tset[c+'cat']['deepSourceId_1']
		tset[c]['dchi2spline'] = np.empty(len(objids),dtype=np.float32)
		for i,objid in enumerate(objids):
			lc = tset[c]['lcs'][objid]
			sfits[objid] = lcfit.fit_spline(lc['gri'])
			mchi2 = np.sum([tset[c]['medianfit'][objid][b]['chi2'] 
			                  for b in 'gri'])
			tset[c]['dchi2spline'][i] = mchi2 - sfits[objid]['chi2']
		tset[c]['splinefit'] = sfits
	return tset

def plotone(tset,objclass,objid,**kwargs):
	if objid < 1e4:
		i = objid
		objid = tset[objclass+'cat']['deepSourceId_1'][i]
	else:
		i = np.where(tset[objclass+'cat']['deepSourceId_1']==objid)[0]
	lc = tset[objclass]['lcs'][objid]
	alc = tset[objclass]['annual'][objid]
	mfit = tset[objclass]['medianfit'][objid]
	fig = lcplot.lcplot(lc,alc,mfit,
	                    (objid,tset[objclass+'cat']['iPsfMag'][i],
               	         tset[objclass+'cat']['z'][i]),**kwargs)
	try:
		lcplot.add_spline_curve(fig,lc,mfit,
		                        tset[objclass]['splinefit'][objid],
		                        tset[objclass]['dchi2spline'][i])
	except KeyError:
		pass
	return fig

