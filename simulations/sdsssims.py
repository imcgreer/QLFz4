#!/usr/bin/env python

from astropy.cosmology import WMAP9
from simqso import qsoSimulation, lumfun

lumGrid = {
  'GridType':'LuminosityRedshiftGrid',
#  'GridFileName':simName+'_grid',
  'mRange':(-27.0,-22.5,0.25),
  'zRange':(2.5,4.6,0.1),
  'nPerBin':5,
  'LumUnits':'M1450',
}

fluxGrid = {
  'GridType':'FluxRedshiftGrid',
  'mRange':(18.5, 23.5, 0.5),
  'zRange':(2.5,4.5,0.1),
  'nPerBin':25,
  'ObsBand':'SDSS-i',
  'RestBand':1450.,
}

czFluxGrid = {
  'GridType':'FluxRedshiftGrid',
  'LFSampledGrid':True,
  'QLFmodel':lumfun.DoublePowerLawLF(-6.0,-25.0,-1.5,-3.0),
  'mRange':(16.5,22.6,0.5),
  'zRange':(2.0,4.1,0.2),
  'nPerBin':5,
  'ObsBand':'SDSS-i',
  'RestBand':1450.,
}

lfFluxGrid = {
  'GridType':'LuminosityFunction',
  'QLFmodel':lumfun.DoublePowerLawLF(-6.0,-25.0,-1.5,-3.0),
  'QLFargs':{'skyArea':10.},
  'mRange':(16.5,22.6),
  'zRange':(2.0,4.1),
  'ObsBand':'SDSS-i',
  'RestBand':1450.,
}

#lum = True
lum = False
lf = False
cz = False

if lf:
	name,grid = ('testnew_qlflum',None) if lum else ('testnew_qlfflux',lfFluxGrid)
elif cz:
	name,grid = ('testnew_czlum',None) if lum else ('testnew_czflux',czFluxGrid)
else:
	name, grid = ('testnew_lum',lumGrid) if lum else ('testnew_flux',fluxGrid)

simpars = {
  'FileName':name,
  'waveRange':(3000.,10000.),
  'SpecResolution':500,
  'DispersionScale':'logarithmic',
  'GridParams':grid,
  'Cosmology':WMAP9,
  'ForestParams':{
    'FileName':name+'_forest',
    'ForestModel':'Worseck&Prochaska2011',
    'ForestType':'Sightlines', # 'Grid' 'OneToOne'
#    'GridzBins':(0.1,4.51,0.25),
    'zRange':(0.0,4.5),
    'NumLinesOfSight':25,
    'Rmin':30000.,
  },
  'QuasarModelParams':{
    'ContinuumParams':{
      'ContinuumModel':'GaussianPLawDistribution',
      'PowerLawSlopes':[(-0.5,0.2),1100,(-0.3,0.2),
                    5700,(-0.78,0.3),10850,(-1.81,0.3),22300,(-1.03,0.3)],
    },
    'EmissionLineParams':{
	  #'EmissionLineModel':'FixedVdBCompositeLines',
	  'EmissionLineModel':'VariedEmissionLineGrid',
      'fixLineProfiles':False,
      'minEW':0.0,
      'EmLineIndependentScatter':False,
      'scaleEWs':{'LyAb':1.6,'LyAn':1.6,
                  'CIVb':0.75,'CIVn':0.75,
                  'CIII]b':0.8,'CIII]n':0.8,
                  'MgIIb':0.8,'MgIIn':0.8},
    },
    'IronEmissionParams':{
	  'FeScalings':[(0,1540,0.5),(1540,1680,2.0),(1680,1868,1.6),
                    (1868,2140,1.0),(2140,3500,1.0)],
    },
    'DustExtinctionParams':{
      #'DustExtinctionModel':'Fixed E(B-V)',
      'DustExtinctionModel':'Exponential E(B-V) Distribution',
      'DustModelName':'SMC',
      'E(B-V)':0.033,
      #'DustLOSfraction':1.0,
    },
    'HostGalaxyParams':{
    },
  },
  'PhotoMapParams':{
    'PhotoSystems':[
      ('SDSS','Legacy'),
      ('UKIRT','UKIDSS_LAS'),
    ]
  },
#  'maxFeatureIter':1,
}

def comparefluxtolumgrids():
	from astropy.io import fits
	import matplotlib.pyplot as plt
	f = fits.getdata('testnew_flux.fits')
	l = fits.getdata('testnew_lum.fits')
	fM = f.M.reshape(12,10,5)
	fm = f.synMag.reshape(12,10,5,5)
	lM = l.M.reshape(12,10,5)
	lm = l.synMag.reshape(12,10,5,5)
	plt.figure()
	plt.scatter(lM[:,0,:],lm[:,0,:,3],c='b')
	plt.scatter(fM[:,0,:],fm[:,0,:,3],c='r')
	plt.figure()
	plt.scatter(l.z,lM,c='b')
	plt.scatter(f.z,fM,c='r')
	plt.figure()
	plt.scatter(l.z,lm[...,3],c='b')
	plt.scatter(f.z,fm[...,3],c='r')

if __name__=='__main__':
	qsoSimulation(simpars,noPhotoMap=True,saveSpectra=True,writeFeatures=True)#,forestOnly=True)


