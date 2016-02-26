"""runAddVar.py
Commandline task to read in output of simqso, add variability, and write csv for classification.

python runAddVar.py testnew_flux.fits
"""
import argparse

import matplotlib
import os
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import QLFz4.gmm
import astropy.io.fits as fits
import scipy
import seaborn as sns

DEBUG = False

###
# To do:
# Recreate training sample of quasars to check that no cuts were applied
# Take command line arguments instead of hard coding.
###

DATAROOT = "/Users/yusra/z4ML/"
IN_SIMSFILENAME = '/Users/yusra/SDSS_flux.fits'
OUT_SIMS_W_VAR_FILENAME = 'SDSS_flux_var.csv'


def mag2FluxColor(magColor):
    return np.power(10, -0.4*(magColor))


def flux2MagColor(fluxColor):
    return -2.5*(np.log10(fluxColor))


def boxcox(x, lam, offset):
    if lam == 0:
        return np.log10(x + offset)
    else:
        return (np.power(x + offset, lam) - 1.)/lam


def unBoxcox(y, lam, offset):
    if lam == 0:
        return np.power(10, y) - offset
    else:
        return np.power(y*lam + 1, 1./lam) - offset


#read in empirical data which to generate Gaussian Mixture
M = pd.Series.from_csv(os.path.join(DATAROOT, 'M1450.csv'), index_col=(0), header=None)
z = pd.Series.from_csv(os.path.join(DATAROOT, 'z.csv'), index_col=(0), header=None)
X = pd.DataFrame.from_csv(os.path.join(DATAROOT, 'X_16feat.csv'), index_col=(0),)
X = X.drop('class', 1)


#take only quasars and convert colors to flux ratios
Qraw = pd.DataFrame(X[z > 2.5])
Qraw['z'] = z[z > 2.5]
Qraw['M'] = M[z > 2.5]
for col in ['ug', 'gr', 'ri', 'iz']:
    Qraw[col] = mag2FluxColor(Qraw[col])

#open sims

sims = fits.open(IN_SIMSFILENAME)[1].data
simsDForig = pd.DataFrame()
simsDForig['M'] = sims['M']
simsDForig['z'] = sims['z']
simsDForig['ug'] = sims['synFlux'][:,0]/sims['synFlux'][:,1]
simsDForig['gr'] = sims['synFlux'][:,1]/sims['synFlux'][:,2]
simsDForig['ri'] = sims['synFlux'][:,2]/sims['synFlux'][:,3]
simsDForig['iz'] = sims['synFlux'][:,3]/sims['synFlux'][:,4]

#transform values so that each feature is close to gaussian as possible:
#boxcox transform:  y = (x^lambda - 1) / lambda for lambda > 0 and y = log(x) for lambda = 0
Q = pd.DataFrame()
transformDict = {}
for i, col in enumerate(Qraw.columns):
    transformDict[col] = {}
    transformDict[col]['offset'] = min(-0.000001, np.min(Qraw[col])*1.01)
    if col in simsDForig.columns:
        transformDict[col]['offset'] = min(transformDict[col]['offset'], np.min(simsDForig[col])*1.01)
    transformDict[col]['lam'] = scipy.stats.boxcox_normmax(Qraw[col] - transformDict[col]['offset'])
    Q[col] = boxcox(Qraw[col].values, transformDict[col]['lam'], - transformDict[col]['offset'])

#plot normality
if DEBUG:
    N = len(Q.columns)
    fig = plt.figure(figsize=(12, 8))
    for i, col in enumerate(Q.columns):
        ax = fig.add_subplot(np.ceil(np.sqrt(N)), np.ceil(N/np.ceil(np.sqrt(N))), i+1)
        #stats.probplot(boxcox(Q[col], lam, -offset), dist=stats.norm, plot=ax)
        ax.hist(Q[col].values, alpha=0.5)
        #ax.hist(finalDF[col].dropna().values, alpha = 0.5)
        ax.set_title(col)

varModel = QLFz4.gmm.GMM()
varModel.fit(Q, criteria='aic')

independentVariables = ['M', 'z', 'ug', 'gr', 'ri', 'iz']

simsDF = pd.DataFrame()
for i, col in enumerate(simsDForig.columns):
    simsDF[col] = boxcox(simsDForig[col].values, transformDict[col]['lam'],
                         -transformDict[col]['offset'])

orderedDependentVariables = []
orderedIndependentVariables = []
for col in Q.columns:
    if col in independentVariables:
        orderedIndependentVariables.append(col)
    else:
        orderedDependentVariables.append(col)

newDF = pd.DataFrame()
arr = np.zeros([len(simsDF), 11])
for index, s in simsDF.iterrows():
    new_cen, new_ccovs, new_mc = varModel.cond_dist(s[Q.columns].values)
    conditionalGMM = QLFz4.gmm.gmmFromParams(np.array(new_cen),
                                             np.array(new_ccovs), np.array(new_mc))
    sampled = conditionalGMM.sample()
    arr[index] = sampled[0]

newDF = pd.DataFrame(arr, columns=orderedDependentVariables)
joinedDF = simsDF.join(newDF)[Q.columns]
#transform back to physical units

finalDF = pd.DataFrame()
for i, col in enumerate(joinedDF.columns):
    print col
    finalDF[col] = unBoxcox(joinedDF[col].values, transformDict[col]['lam'],
                            -transformDict[col]['offset'])



if DEBUG:
    tempDF = finalDF.copy(deep=True)
    tempDF['sim'] = False
    temp2 = Qraw.copy(deep=True)
    temp2['sim'] = True
    toPlot = pd.concat((tempDF, temp2))
    p = sns.PairGrid(toPlot, hue='sim')
    p.map_lower(plt.scatter, alpha=0.2)
    p.savefig('simulatedQuasarsNov30.png')

#NOTE: COLORS ARE FLUXE RATIOS!!!!
finalDF.to_csv(OUT_SIMS_W_VAR_FILENAME)

#if __name__=='__main__':
#    parser = argparse.ArgumentParser(description="""Read in output of simqso,
##                                     add variability parameters, and write csv for classification""")
#    parser.add_argument('inFile',
#                         help='input filename: fits table from simqso.qsoSimulation')
#
#    parser.add_argument('outFile',
#                         help='output filename')
#args = parser.parse_args()
