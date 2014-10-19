#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import scoreatpercentile

import lcfit

def lcplot(lc,alc,medianfit,objdat,fig=None):
	clr = {'g':'g','r':'r','i':'indigo'}
	objid,imag,z = objdat
	if fig is None:
		fig = plt.figure(figsize=(12,5.0))
		plt.subplots_adjust(0.04,0.04,0.99,0.95,0.20,0.24)
	else:
		fig.clear()
	ax = None
	for pnum,b in enumerate('gri',start=1):
		ax = plt.subplot(3,1,pnum,sharex=ax)
		ii = np.where(~lc[b]['mask'])[0]
		plt.errorbar(lc[b]['mjd'][ii],lc[b]['flux'][ii],
		             1/np.sqrt(lc[b]['ivar'][ii]),
		             fmt='s',color=clr[b])
		for j in np.where(~alc[b]['mask'])[0]:
			err = 1/np.sqrt(alc[b]['ivar'][j])
			span = Rectangle((alc[b]['mjd'][j]-90,alc[b]['flux'][j]-err),
			                 180,2*err,facecolor=clr[b],edgecolor='none',
			                 alpha=0.5)
			ax.add_patch(span)
		plt.axhline(medianfit[b]['median'],color=clr[b])
		plt.text(0.05,0.05,r'$\chi^2=%.1f/%d$' % 
		                   (medianfit[b]['chi2'],medianfit[b]['ndof']),
		         transform=ax.transAxes)
		if b=='g':
			ax.set_title('$%d, m_i=%.2f, z=%.2f$' % (objid,imag,z),
			             size=14)
	return fig

def add_spline_curve(fig,lc,medianfit,splinefit,deltaChi2=None):
	xx = np.arange(51500,55000,10)
	for ax,b in zip(fig.axes,'gri'):
		ylim = ax.get_ylim()
		yy = splinefit['spfit'](xx) * medianfit[b]['median']
		ax.plot(xx,yy,c='0.2')
		if b == 'g':
			ax.text(0.18,0.05,r'$\chi^2(spline)=%.1f/%d$' % 
			                   (splinefit['chi2'],splinefit['ndof']),
			         transform=ax.transAxes)
			if deltaChi2 is not None:
				ax.text(0.34,0.05,r'$(%.1f)$' % deltaChi2,
			         transform=ax.transAxes)
		elif b == 'i':
			mchi2 = np.sum([medianfit[b]['chi2'] for b in 'gri'])
			mndof = np.sum([medianfit[b]['ndof'] for b in 'gri'])
			ax.text(0.18,0.05,r'$\chi^2(gri)=%.1f/%d$' % (mchi2,mndof),
			         transform=ax.transAxes)
		ax.set_xlim(51500,54600)
		ax.set_ylim(*ylim)


