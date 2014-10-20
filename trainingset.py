#!/usr/bin/env python

import os
import numpy as np
from astropy.io import fits

import lcfit
import lcplot

# Directory containing lightcurve data
# Should have this structure:
#  2014October/
#    TRAININGwgriReduced.fits
#    trainingDataFiles_flux/
#      GALAXY/
#        process_[gri].dat
#        [gri]/
#          <DEEPSOURCEID>.dat
#      QSO/...
#      STAR/...
#
qlfz4dir = os.path.join(os.environ['QLFZ4DATA'], '2014October')

def load_training_catalog():
	trainingdataf = os.path.join(qlfz4dir,'TRAININGwgriReduced.fits')
	return fits.getdata(trainingdataf,1)

def load_all():
	'''Load the lightcurves and catalog data for the full training set,
	   divided by object class.
	'''
	tset = {}
	tset['catalog'] = load_training_catalog()
	for c in ['QSO','GALAXY','STAR']:
		lcdir = os.path.join(qlfz4dir,'trainingDataFiles_flux',c)
		tset[c] = lcfit.load_lsst_stripe82_lightcurves(lcdir)
		ii = [np.where(tset['catalog']['deepSourceId_1']==objid)[0][0]
		          for objid in tset[c]['objs']['id']]
		tset[c+'cat'] = tset['catalog'][ii]
	return tset

def add_spline_fit(tset):
	'''Add the spline fits to the training set structure.'''
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
	'''Plot an object from the training set. E.g.,
	     plotone(tset,'QSO',0) # plots the first QSO in the list
	     ii = where(tset['QSOcat']['z']>4)[0]; plotone(tset,'QSO',ii[0])
	       # plot a high-z quasar 
	   keyword arguments are passed to lcplot()
	   If the spline fit exists it will be plotted.
	'''
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

def compare_deltachi2(tset):
	import matplotlib.pyplot as plt
	plt.figure()
	for c,clr,sym in zip(['GALAXY','STAR','QSO'],'grb','^so'):
		plt.scatter(tset[c+'cat']['median_r'],tset[c]['dchi2spline'],
		            facecolor=clr,edgecolor='none',alpha=0.8,s=30,marker=sym)
	plt.yscale('log')
	plt.ylim(1,5e4)

