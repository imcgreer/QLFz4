#!/usr/bin/env python

from numpy import uint32
from astropy.cosmology import WMAP9
from simqso import qsoSimulation, lumfun

grids = {
  'SDSS_Lum':{
    'GridType':'LuminosityRedshiftGrid',
    'mRange':(-27.0,-22.09,0.1),
    'zRange':(3.25,4.5,0.05),
    'nPerBin':10,
    'LumUnits':'M1450',
  },
  'SDSS_Flux':{
    'GridType':'FluxRedshiftGrid',
    'mRange':(18.5, 23.41, 0.1),
    'zRange':(3.25,4.5,0.05),
    'nPerBin':10,
    'ObsBand':'SDSS-i',
    'RestBand':1450.,
  },
  'SDSS_LF':{
    'GridType':'LuminosityFunction',
    'QLFmodel':lumfun.DoublePowerLawLF(-6.0,-25.0,-1.5,-3.0),
    'QLFargs':{'skyArea':10.},
    'mRange':(16.5,22.6),
    'zRange':(2.0,4.1),
    'ObsBand':'SDSS-i',
    'RestBand':1450.,
  },
}

simpars = {
#  'FileName':name,
  'waveRange':(3000.,10000.),
  'SpecDispersion':500,
  'DispersionScale':'logarithmic',
#  'GridParams':grid,
  'Cosmology':WMAP9,
  'ForestParams':{
#    'FileName':name+'_forest',
    'ForestModel':'Worseck&Prochaska2011',
    'ForestType':'Sightlines', # 'Grid' 'OneToOne'
    'zRange':(0.0,4.5),
    'NumLinesOfSight':2000,
    'Rmin':30000.,
  },
  'QuasarModelParams':{
    'ContinuumParams':{
      'ContinuumModel':'GaussianPLawDistribution',
      'PowerLawSlopes':[(-0.5,0.2),1100,(-0.3,0.2),
                    5700,(-0.78,0.3),10850,(-1.81,0.3),22300,(-1.03,0.3)],
    },
    'EmissionLineParams':{
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
      'DustExtinctionModel':'Exponential E(B-V) Distribution',
      'DustModelName':'SMC',
      'E(B-V)':0.033,
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
}

if __name__=='__main__':
	import sys
	survey,gridType = sys.argv[1:]
	simpars['FileName'] = survey+'_'+gridType
	simpars['GridParams'] = grids[survey+'_'+gridType]
	if gridType == 'LF':
		simpars['ForestParams']['FileName'] = survey+'_LF_forest'
	else:
		simpars['ForestParams']['FileName'] = survey+'_grid_forest'
		# makes the lum and flux grids use the same forest
		simpars['RandomSeed'] = uint32(hash(survey))
	# dust not included for now
	del simpars['QuasarModelParams']['DustExtinctionParams']
	qsoSimulation(simpars,writeFeatures=True)#,forestOnly=True)


