#!/usr/bin/env python

import os,sys
from copy import copy
import pprint
from numpy import uint32,where
from astropy.cosmology import WMAP9
from simqso import qsoSimulation,lumfun,sqanalysis
from simqso.sqgrids import VariedEmissionLineGrid

simdir = './'

grids = {
  'CFHTLS_Wide_LumGrid':{
    'GridType':'LuminosityRedshiftGrid',
    'mRange':(-27.5,-21.59,0.1),
    'zRange':(3.4,4.61,0.05),
    'nPerBin':100,
    'LumUnits':'M1450',
  },
  'CFHTLS_Wide_FluxGrid':{
    'GridType':'FluxRedshiftGrid',
    'mRange':(18.6,24.51,0.1),
    'zRange':(3.4,4.61,0.05),
    'nPerBin':100,
    'ObsBand':'CFHT-i',
    'RestBand':1450.,
  },
  'CFHTLS_Wide_McG13QLFpoints':{
    'GridType':'LuminosityFunction',
    'QLFmodel':lumfun.DoublePowerLawLF(lambda z: -8.94 - 0.47*(z-6.),
                                       -27.21,-2.03,-4.00),
    'QLFargs':{'skyArea':150.},
    'mRange':(17.0,24.0),
    'zRange':(3.4,4.6),
    'ObsBand':'CFHT-i',
    'RestBand':1450.,
  },
}

fiducial_linetweak_model = {
  'ContinuumParams':{
    'ContinuumModel':'GaussianPLawDistribution',
    'PowerLawSlopes':[(-1.5,0.3),1100,(-0.5,0.3),
                  5700,(-0.37,0.3),9730,(-1.7,0.3),22300,(-1.03,0.3)],
  },
  'EmissionLineParams':{
    'EmissionLineModel':'VariedEmissionLineGrid',
    'fixLineProfiles':False,
    'minEW':0.0,
    'EmLineIndependentScatter':False,
    'scaleEWs':{'LyAb':1.1,'LyAn':1.1,
                'CIVb':0.75,'CIVn':0.75,
                'CIII]b':0.8,'CIII]n':0.8,
                'MgIIb':0.8,'MgIIn':0.8},
  },
  'IronEmissionParams':{
    'FeScalings':[(0,1540,0.5),(1540,1680,2.0),(1680,1868,1.6),
                  (1868,2140,1.0),(2140,3500,1.0)],
  },
}

def PaperI_forest(gridname):
	return {
	  'FileName':gridname+'_forest',
	  'ForestModel':'McGreer+2013',
	  'ForestType':'Sightlines', 
	  'zRange':(0.0,6.1),
	  'NumLinesOfSight':1000,
	  'Rmin':30000.,
	}

cfhtls_photo = {
  'PhotoSystems':[
    ('CFHT','CFHTLS_Wide'),
    ('UKIRT','UKIDSS_DXS','JK'),
    ('SDSS','Stripe82'),
  ]
}

base_params = {
  'FileName':None,
  'waveRange':(3000.,3e4),
  'SpecDispersion':500,
  'DispersionScale':'logarithmic',
  'GridParams':None,
  'GridFileName':None,
  'Cosmology':WMAP9,
  'ForestParams':None,
  'QuasarModelParams':None,
  'maxFeatureIter':4,
  'PhotoMapParams':None,
  'RandomSeed':1,
}

def gridName(survey,gridType):
	_gridName = {'flux':survey+'_FluxGrid',
	             'luminosity':survey+'_LumGrid',
	             'QLF':survey+'_McG13QLFpoints',
	}[gridType]
	return _gridName

def simName(survey,gridType):
	_simName = {'flux':survey+'_Flux',
	            'luminosity':survey+'_Lum',
	            'QLF':survey+'_McG13QLF',
	}[gridType]
	return _simName

def get_params(survey,gridtype,quasarmodel):
	p = copy(base_params)
	gridname = gridName(survey,gridtype)
	forestname = survey + '_Grid' if gridtype!='QLF' else gridname
	p['GridFileName'] = gridname
	p['GridParams'] = grids[gridname]
	p['FileName'] = simName(survey,gridtype)
	#p['ForestParams'] = PaperI_forest(forestname)
	p['ForestParams'] = {
	  'FileName':gridname+'_forest',
	  'ForestModel':'Worseck&Prochaska2011',
	  'zRange':(0.0,4.6),
	  'ForestType':'Sightlines',
	  #'NumLinesOfSight':2000,
	  'NumLinesOfSight':100, # however, 100 will run a lot faster for testing
	  'Rmin':30000.,
	}
	p['QuasarModelParams'] = copy(quasarmodel)
	p['PhotoMapParams'] = cfhtls_photo
	# this is a bit of a hack, but because each survey shares a random seed,
	# as long as the shape of the luminosity and flux grids are the same,
	# the random sampling of the redshifts should be identical, and thus
	# the forest model can be shared and only generated once for each survey
	p['RandomSeed'] = uint32(hash(survey))
	pprint.pprint(p)
	return p

if __name__=='__main__':
	#survey,gridtype = sys.argv[1:3]
	#try:
	#	extraargs = sys.argv[3:]
	#except:
	#	extraargs = ()
	survey = 'CFHTLS_Wide'
	gridtype = 'flux'
	quasarmodel = fiducial_linetweak_model
	forestscale = None
	p = get_params(survey,gridtype,quasarmodel)#,forestscale,extraargs)
	if survey=='test':
		simdir = './'
	qsoSimulation(p,
#	              noPhotoMap=True,
#	              saveSpectra=True,
	              writeFeatures=True,
#	              forestOnly=True,
	              outputDir=simdir)

