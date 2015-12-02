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
	simpars['QuasarModelParams'] = fiducial_linetweak_model
	qsoSimulation(simpars,writeFeatures=True)#,forestOnly=True)


