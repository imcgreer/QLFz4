#!/usr/bin/env python

import os
import numpy as np
from astropy.table import table
from astropy import units as u
from astropy.coordinates import SkyCoord,match_coordinates_sky
from astropy.table import Table,hstack

from simqso import hiforest,sqbase,sqgrids,sqmodels,sqphoto,sqrun,sqanalysis
from simqso.lumfun import QuasarSurvey

def read_s82_catalog(paperdir=None):
	if paperdir is None:
		paperdir = os.path.join(os.environ['HOME'],'research','paper_qlfz4')
	rows = []
	reading = False
	with open(os.path.join(paperdir,'catalog.tex')) as catf:
		for line in catf:
			if 'enddata' in line:
				break
			elif reading:
				dat = line.strip().split('&')
				coo = dat[0] + ' ' + dat[1]
				name = dat[-1][:-2]
				coo = SkyCoord(coo, unit=(u.hourangle, u.deg))
				i1 = dat[4].find('$') + 1
				i2 = dat[4].find("\\")
				try:
					iAB = float(dat[2])
				except:
					iAB = 0.0
				try:
					M1450 = float(dat[3])
				except:
					M1450 = 0.0
				z = float(dat[4][i1:i2])
				rows.append( (coo.ra.value,coo.dec.value,iAB,
				              M1450,z,dat[5],name) )
			elif 'startdata' in line:
				reading = True
	# non-quasars being added by hand
	rows.append( (333.545426,0.807974, 22.3731, 0, -1, '???', 'w08') )
	rows.append( (334.664088,0.968030, 22.3456, 0, -1, '???', 'w05') )
	tab = Table(rows=rows,
	            names=('ra','dec','iAB','M1450','z','obsdate','name'))
	return tab

def match_target_tab(paperdir=None,obsdir=None):
	zcat = read_s82_catalog()
	if obsdir is None:
		obsdir = os.path.join(os.environ['HOME'],'research','LSST',
		                      'Stripe82','2014October')
	cans = Table.read(os.path.join(obsdir,'cfhtw4_candidates_v3.fits'))
	cans = cans[ (cans['dec']<1.25) & (cans['mags'][:,3]>21.5) & 
	             (cans['mags'][:,3]<22.5) ]
	c1 = SkyCoord(cans['ra'],cans['dec'],unit=(u.deg,u.deg))
	c2 = SkyCoord(zcat['ra'],zcat['dec'],unit=(u.deg,u.deg))
	idx,d2d,d3c = match_coordinates_sky(c1,c2)
	ii = np.where(d2d.arcsec < 3.0)[0]
	assert np.all(ii == np.arange(len(cans)) )
	cfhtw4 = hstack([cans,zcat['M1450','z','obsdate','name'][idx[ii]]])
	cfhtw4.write('cfhtw4qsos.fits',overwrite=True)

def run_w4_sim(cosmo):
	wave = sqbase.fixed_R_dispersion(3000,7e4,500)
	#
	if False:
		m = sqgrids.AbsMagVar(sqgrids.UniformSampler(-25,-23),1450)
		units = 'luminosity'
	else:
		m = sqgrids.AppMagVar(sqgrids.UniformSampler(21.2,22.8),'CFHT-i')
		units = 'flux'
	z = sqgrids.RedshiftVar(sqgrids.UniformSampler(3.5,4.5))
	#
	qsoGrid = sqgrids.QsoSimGrid([m,z],(20,20),100,cosmo=cosmo,units=units)
	#
	tmp_m2M = lambda z: sqbase.mag2lum('SDSS-i',1450,z,cosmo)
	absMagVals = sqgrids.FixedSampler(qsoGrid.appMag-tmp_m2M(qsoGrid.z))
	qsoGrid.addVar(sqgrids.AbsMagVar(absMagVals,1450))
	#
	bossDr9model = sqmodels.get_BossDr9_model_vars(qsoGrid,wave,10)
	qsoGrid.addVars(bossDr9model)
	#
	qsoGrid.loadPhotoMap([('CFHT','CFHTLS_Wide'),('UKIRT','UKIDSS_DXS')])
	#
	_ = sqrun.buildSpectraBySightLine(wave,qsoGrid,verbose=10)
	#
	photoData = sqphoto.calcObsPhot(qsoGrid.synFlux,qsoGrid.photoMap)
	qsoGrid.addData(photoData)
	#
	qsoGrid.write('cfhtlsw4_z4sims.fits')
	return qsoGrid

class W4ColorSel(sqanalysis.SelectionFunction):
	photo_complete = 0.9
	spec_complete = 1.0
	@staticmethod
	def w4_cfhtls_dxs_sel(mags,errs,fluxes,fluxerrs,bands):
		bj = lambda b: bands.index(b)
		u,g,r,i,z,J,H,K = [ mags[:,bj(b)]
		                         for b in [ 'CFHT-CFHTLS_Wide-'+_b
		                                       for _b in 'ugriz' ] +
		                                  [ 'UKIRT-UKIDSS_DXS-'+_b
		                                       for _b in 'JHK' ] ]
		u_snr = ( fluxes[:,bj('CFHT-CFHTLS_Wide-u')] / 
		          fluxerrs[:,bj('CFHT-CFHTLS_Wide-u')] ) 
		#s = ( (i>21.5) & (i<22.5) )
		s = i > 0
		return s
	def __init__(self,*args,**kwargs):
		super(W4ColorSel,self).__init__(*args,**kwargs)
		self.selector = self.w4_cfhtls_dxs_sel

class m2M(object):
	def __init__(self,qsoGrid,kcorr):
		self.qsoGrid = qsoGrid
		self.kcorr = kcorr
		self.cosmo = self.qsoGrid.cosmo
	def __call__(self,m,z,inverse=False):
		return self.qsoGrid.distMod(z) + self.kcorr(m,z,inverse)

def CFHTW4_ColorSample_QLF(cosmo):
	qsos = Table.read('cfhtw4qsos.fits')
	simqsos = sqgrids.QsoSimGrid(cosmo=cosmo)
	simqsos.read('cfhtlsw4_z4sims.fits')
	# XXX
	simqsos.loadPhotoMap([('CFHT','CFHTLS_Wide'),('UKIRT','UKIDSS_DXS')])
	kcorr = sqanalysis.SimKCorr(simqsos)
	kcorr.calc_grid('CFHT-i')
	w4_m2M = m2M(simqsos,kcorr)
	selfun = W4ColorSel(simqsos,w4_m2M)
	selfun.calc_grid()
	# see w4sel.py
	w4survey_area = ( (335.714 - 332.20) * (1.25 - -0.15)  +
	                  (335.714 - 333.05) * (-0.15 - -1.022) )
	w4survey = QuasarSurvey(qsos['mags'][:,3],qsos['z'],22.5,
	                        w4survey_area,w4_m2M)
	w4survey.set_selection_function(selfun)
	lf = w4survey.calcBinnedLF([-24.5,-23.9,-23.3],[3.6,4.4])
	print lf

if __name__=='__main__':
	from astropy.cosmology import WMAP9
	run_w4_sim(WMAP9)
	CFHTW4_ColorSample_QLF(WMAP9)

