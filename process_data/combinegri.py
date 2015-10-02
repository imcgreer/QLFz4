import pandas as pd
import os
import numpy as np
import argparse
from aggregates import *

def runPatch(root, patchRange):
    lc = {}
    for f in 'gri':
        path = os.path.join(root, f, patchRange)
        lc[f] = pd.DataFrame.from_csv(os.path.join(path, 'annualAvg.csv'),
                                      index_col=(0,1), parse_dates = False)
        keys = lc[f].columns.values
        lc[f].columns = [s+'_%s'%(f) for s in keys]
    df = lc['g'].join(lc['r']).join(lc['i'])
    grouped  = df.groupby(level=0)
    metrics = grouped.apply(computeMetricsGRI)
    patchDir = os.path.join(root, 'gri', patchRange)
    if not os.path.exists(patchDir):
        os.makedirs(patchDir)
    metrics.to_csv(os.path.join(patchDir, 'chi2.csv'))


def combine_gri(gri, griErr, weights = [0.15,0.45,0.4], combine=None):
    """gri: an array of length 3 containing flux arrays (one element per year)
       griErr: array of length 3 containing flux uncertainies (sigma)
    """
    variances = griErr**2
    if combine is not None:
        mfluxes = gri.mean(axis=1)[:,np.newaxis]
        if combine=='subtract_mean':
            fluxes = (gri-mfluxes) + (mfluxes-mfluxes[1])
        elif combine=='divide_mean':
            fluxes = gri/mfluxes-1
            variances /= mfluxes**2
    else:
        fluxes = gri
    ivars = 1/np.sum(variances.T * np.array(weights)**2, axis=1)
    return np.average(fluxes, weights=weights, axis=0), ivars


def gri_chi2(gri, griErr, weights = [0.15,0.45,0.4], combine=None):
    y, yIvars = combine_gri(gri, griErr, weights=weights, combine=combine)
    N = len(y) #added Oct 7th after Ian noticed values too high.
    return np.sum((y-np.mean(y))**2 * yIvars)/(N-1)





def computeMetricsGRI(group):
    return pd.Series({
    'chi2griClippedWMean': gri_chi2(np.array([group['clippedWeightedMean_g'].values,
                                              group['clippedWeightedMean_r'].values,
                                              group['clippedWeightedMean_i'].values]),
                                    np.array([group['clippedWeightedStdCorr_g'].values,
                                              group['clippedWeightedStdCorr_r'].values,
                                              group['clippedWeightedStdCorr_i'].values])
                                    ),
    'chi2griWMean': gri_chi2(np.array([group['clippedWeightedMean_g'].values,
                                       group['clippedWeightedMean_r'].values,
                                       group['clippedWeightedMean_i'].values]),
                             np.array([group['clippedWeightedStdCorr_g'].values,
                                       group['clippedWeightedStdCorr_r'].values,
                                       group['clippedWeightedStdCorr_i'].values])
                            ),
    'chi2griMedian':gri_chi2(np.array([group['median_g'].values,
                                       group['median_r'].values,
                                       group['median_i'].values]),
                             np.array([group['medianSE_g'].values,
                                       group['medianSE_r'].values,
                                       group['medianSE_i'].values])
                            ),
    'chi2griClippedWMeanDivideMean': gri_chi2(np.array([group['clippedWeightedMean_g'].values,
                                                        group['clippedWeightedMean_r'].values,
                                                        group['clippedWeightedMean_i'].values]),
                                              np.array([group['clippedWeightedStdCorr_g'].values,
                                                        group['clippedWeightedStdCorr_r'].values,
                                                        group['clippedWeightedStdCorr_i'].values]),
                                                combine='divide_mean'
                                    ),
    'chi2griWMeanDivideMean': gri_chi2(np.array([group['clippedWeightedMean_g'].values,
                                                 group['clippedWeightedMean_r'].values,
                                                 group['clippedWeightedMean_i'].values]),
                                       np.array([group['clippedWeightedStdCorr_g'].values,
                                                 group['clippedWeightedStdCorr_r'].values,
                                                 group['clippedWeightedStdCorr_i'].values]),
                                       combine='divide_mean'
                                       ),
    'chi2griMedianDivideMean':gri_chi2(np.array([group['median_g'].values,
                                                 group['median_r'].values,
                                                 group['median_i'].values]),
                                       np.array([group['medianSE_g'].values,
                                                 group['medianSE_r'].values,
                                                 group['medianSE_i'].values]),
                                       combine='divide_mean'
                                        )
    })




""" From the original chi2calc.py from 2013
def stack_fluxes(dat,cadence='full',weights=[0.15,0.45,0.4],
                 noivarwt=False,sn2weight=False,combine=None):
    rv = {}
    weights = np.array(weights).reshape(1,3)
    for i,sid in enumerate(dat['objs']['deepSourceId']):
        fluxes = dat[sid][cadence+'_aligned']['fluxes']
        ivars = dat[sid][cadence+'_aligned']['ivars']
        if combine is not None:
            mfluxes = fluxes.mean(axis=0)[np.newaxis,:]
            if combine=='subtract_mean':
                fluxes[:] = (fluxes-mfluxes) + (mfluxes-mfluxes[:,1])
            elif combine=='divide_mean':
                fluxes[:] = fluxes/mfluxes-1
                ivars *= mfluxes**2
        x = {}
        if noivarwt:
            x['fluxes'] = np.average(fluxes,weights=weights[0],axis=1)
            x['ivars'] = 1/np.sum(weights**2/ivars,axis=1)
        elif sn2weight:
            snw = (fluxes*np.sqrt(ivars)).mean(axis=0)
            x['fluxes'] = np.average(fluxes,weights=snw,axis=1)
            x['ivars'] = 1/np.sum(snw**2/ivars,axis=1)
        else:
            x['fluxes'],x['ivars'] = np.average(fluxes,
                                                weights=weights*ivars,
                                                axis=1,returned=True)
        x['fluxes'] = x['fluxes'][:,np.newaxis]
        x['ivars'] = x['ivars'][:,np.newaxis]
        x['mjd'] = dat[sid][cadence+'_aligned']['mjd']
        rv[sid] = x
    return rv

"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory",
                        help = "Input directory containing patch subdirectories containing annualAvg.csv")
    parser.add_argument("--firstPatch", type = int,
                        help = "starting patch in range. Looks for file <filter>/<firstPatch>_<lastPatch>/annualAvg.csv")
    parser.add_argument("--lastPatch", type = int,
                        help = "last patch in range. Looks for file <filter>/<firstPatch>_<lastPatch>/annualAvg.csv")
    args = parser.parse_args()

    if args.firstPatch:
        run(args.directory, args.firstPatch, args.lastPatch)
    else:
        patchRanges = os.listdir(os.path.join(args.directory,'i'))
        for patchRange in patchRanges:
            print patchRange
            runPatch(args.directory, patchRange)

