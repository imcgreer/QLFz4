

# same as Ross+13,McGreer+13
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

# add exponential dust screening, discussed in Ross+13
expdust_model = {
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
}

