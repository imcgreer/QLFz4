#!/usr/bin/env python

import os
import numpy as np
from astropy.io import fits
from astrotools.idmstuff import matchlists

import lcfit
import lcplot


#tsetdir = os.environ['HOME'] + \
#      '/data/projects/LSST/Stripe82/2014October/trainingDataFiles_flux/'

tsetdir = os.environ['QLFZ4DATA']+'2014October/trainingDataFiles_flux/'

def load_training_catalog():
	return fits.getdata(tsetdir+'../TRAININGwgriReduced.fits',1)

def load_class(objclass):
	objs = np.genfromtxt(tsetdir+objclass+'/process_i.dat',
	                     usecols=(1,2,3,5),names='ra,dec,z,id',
	                     dtype=('f8','f8','f4','i8'),
	                     converters={5: lambda s: long(s[3:-5])})
	lightcurves,annualcurves = {},{}
	mfits = {}
	for objid in objs['id']:
		lc,blc = {},{}
		mfit = {}
		for b in 'gri':
			lcdatfile = tsetdir+objclass+'/'+b+'/%d.dat' % objid
			lcdat = np.loadtxt(lcdatfile,unpack=True)
			if True:
				# fluxes are ~1e-29, scale them
				lcdat[1:] *= 1e29
			lc[b] = lcfit.make_lcstruct(*lcdat)
			lc[b],blc[b] = lcfit.bin_lc(lc[b])
			mfit[b] = lcfit.fit_median(lc[b])
		lc['gri'] = lcfit.combine_lcstructs([lc[b] for b in 'gri'],
		                                    [mfit[b] for b in 'gri'])
		lightcurves[objid] = lc
		annualcurves[objid] = blc
		mfits[objid] = mfit
	return dict(objs=objs,lcs=lightcurves,annual=annualcurves,medianfit=mfits)

def load_all():
	tset = {}
	tset['catalog'] = load_training_catalog()
	for c in ['QSO','GALAXY','STAR']:
		tset[c] = load_class(c)
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
	mfit = tset[objclass]['medianfit'][objid]
	fig = lcplot.lcplot(lc,mfit,
	                    (objid,tset[objclass+'cat']['iPsfMag'][i],
               	         tset[objclass+'cat']['z'][i]),**kwargs)
	try:
		lcplot.add_spline_curve(fig,lc,mfit,
		                        tset[objclass]['splinefit'][objid],
		                        tset[objclass]['dchi2spline'][i])
	except KeyError:
		pass
	return fig

def matchljcoadd():
	import deepcat
	train = load_training_catalog()
	m = deepcat.match(train['ra'],train['decl'])
	return m

