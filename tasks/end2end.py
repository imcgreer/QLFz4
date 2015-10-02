"""
Demonstrate end to end usage of tasks
"""

from simqso import qsoSimulation
from QLFz4.simulations.sdsssims import simpars


# (1) Train Classifier
#         Input)  Training Data Set
#         Output) .pkl of sklearn classifier
#                 list of feature columns used by classifier
# This step is currently manual


# (2) Simulate quasars and their classifiers using simqso
#         Input)  Simulation paramter file: (SED model, flux/redshift distribution)
#         Output) set of quasars (M, z, ugriz)


# !python sdsssims.py
qsoSimulation(simpars, noPhotoMap=True, saveSpectra=True, writeFeatures=True)


# (3) Add variability params to simulated quasars
#         Input) Output from simqso + Training set + list of feature columns used by classifier
#         Ouput) Set of simulated quasars (M, z, ugriz, variability)


# 4) Run simulated quasars though classifier
#         Input) Set of simulated quasars (M, z, ugriz, variability)
#         Ouput) completeness grid: fraction(M, z)


# 5) Run candidates through classifier
#         Input) candidates
#         Output) N(M, z)


# 6) "QuasarSurvey"
#         Input) N(M, z), fraction(M, z), m2M/M2m, area
#         Output) LF