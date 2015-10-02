from aggregates import *
import os
import pandas as pd
import numpy as np
import argparse




def sigmaClipped(y,yerr,sigmathresh):
    """compute statistics for a sigma-clipped lightcurve
    """
    q25, q50, q75 = np.percentile(y, (25, 50, 75))
    iqrSig = 0.7413*(q75 - q25)
    nsigma = np.abs(y - q50)/iqrSig
    idx = np.where(nsigma <= sigmathresh)
    wm, stdCorr = calcWStdCorrAndMean(y[idx], yerr[idx])
    return {'mean': np.mean(y[idx]),
            'std': np.std(y[idx]),
            'chi2': calcChi2raw(y[idx], yerr[idx]),
            'wm': wm,
            'wStdCorr': stdCorr,
            'N': len(idx[0])}

def calcMedianSE(y, yErr):
    """Compute robust standard error on the median
    """
    N = len(y)
    if N < 2:
        return float(yErr)
    else:
        return 0.7413*(np.percentile(y,75)-np.percentile(y,25))/np.sqrt(N)

def calcMedianAndSE(y, yErr):
    """Compute robust standard error on the median
    """
    N = len(y)
    if N < 2:
        return float(y), float(yErr)
    else:
        q25, q50, q75 = np.percentile(y, (25, 50, 75))
        return q50, 0.7413*(q75-q25)

def mjd2year(t):
    """mjd to Stripe82 Survey year
    """
    t = t - 50925
    return t//365


def computeMetricsByYear(group):
    """Metrics to compute for each object for each year. (annually averaged metrics)
    """
    return pd.Series({
                      'N': group['psfFlux'].count(),
                      'median': group['psfFlux'].median(),
                      'medianSE': calcMedianSE(group['psfFlux'],group['psfFluxErr']),
                      'WeightedMean': calcWM(group['psfFlux'],group['psfFluxErr']),
                      'WeightedStdCorr': calcWStdCorr(group['psfFlux'],group['psfFluxErr']),
                      'clippedWeightedMean': sigmaClippedWM(group['psfFlux'].values,group['psfFluxErr'].values,5),
                      'clippedWeightedStdCorr': sigmaClippedWStdCorr(group['psfFlux'].values,group['psfFluxErr'].values,5),
                       })

def computeMetrics(group):
    """Metrics to compute for each object on full lightcurve
    """
    return pd.Series({'mean': group['psfFlux'].mean(),
                      'N': group['psfFlux'].count(),
                      'chi2': calcChi2raw(group['psfFlux'], group['psfFluxErr']),
                      'median': group['psfFlux'].median(),
                      'sigmaClipped': sigmaClipped(group['psfFlux'].values, group['psfFluxErr'].values, 5),
                      'sigmaG': 0.7413*(np.percentile(group['psfFlux'].values,75)-np.percentile(group['psfFlux'].values,25)),
                      'e50':  group['psfFluxErr'].median(),
                      'e_mean': group['psfFluxErr'].mean(),
                      'WeightedMean': calcWM(group['psfFlux'].values,group['psfFluxErr'].values),
                      'WeightedStdCorr': calcWStdCorr(group['psfFlux'].values,group['psfFluxErr'].values),
                       })


def computeMetricsAggYears(group):
    """Metrics on the annually averaged lightcurves
    """
    return pd.Series({
                      'annAvg_chi2_wm': calcChi2raw(group['WeightedMean'],group['WeightedStdCorr']),
                      'annAvg_chi2_clip_wm': calcChi2raw(group['clippedWeightedMean'],group['clippedWeightedStdCorr']),
                      'annAvg_chi2_median': calcChi2raw(group['median'],group['medianSE']),
                      'annAvg_rms_clip_wm': group['clippedWeightedMean'].std(),
                      'annAvg_e50_clip_wm': group['clippedWeightedStdCorr'].median()
                       })

def computeMetricsByYearLessOp(group):
    """Metrics to compute for each object for each year. (annually averaged metrics)
    """
    return pd.Series({
                      'N': group['psfFlux'].count(),
                      'median': calcMedianAndSE(group['psfFlux'].values,group['psfFluxErr'].values),
                      'WeightedMean': calcWStdCorrAndMean(group['psfFlux'].values,group['psfFluxErr'].values),
                      'clippedWeightedMean': sigmaClippedWStdCorrAndMean(group['psfFlux'].values,group['psfFluxErr'].values,5),
                       })


def run(inDir, outDir,filterName, firstPatch, lastPatch):
    """
    #1) Group by objectId,
    #2) Group by objectId, year and calucate the metrics to generate annually-averaged
    lightcurves, then group by objectId and calculate the metrics.
    #3) Join the full lc metrics and the ann_averaged metrics and print.
    """
    path = os.path.join(outDir, filterName, '%02i_%02i'%(firstPatch, lastPatch))

    if not os.path.exists(path):
        os.makedirs(path)

    filename = '%s%02i_%02i.csv'%(filterName, firstPatch, lastPatch)
    #read data and header separately since processForcedSource printed multiple headers:
    df = pd.read_csv(os.path.join(inDir, filename), index_col=0, comment='d', header=None)
    header = pd.read_csv(os.path.join(inDir, filename), index_col=0, nrows=1)
    df.columns = header.columns
    df.index.name = header.index.name
    df['year'] = pd.Series(mjd2year(df['exposure_time_mid']))
    df = df[df.psfFluxErr != 0]
    grouped = df.groupby('objectId')
    appliedMetrics = grouped.apply(computeMetrics)
    groupedByObjIdYear =  df.groupby(['objectId','year'])
    annAveragedMetrics = groupedByObjIdYear.apply(computeMetricsByYearLessOp)

    annAveragedMetrics['medianSE'] = annAveragedMetrics['median'].apply(lambda x: x[1])
    annAveragedMetrics['median'] = annAveragedMetrics['median'].apply(lambda x: x[0])
    annAveragedMetrics['WeightedStdCorr'] = annAveragedMetrics['WeightedMean'].apply(lambda x: x[1])
    annAveragedMetrics['WeightedMean'] = annAveragedMetrics['WeightedMean'].apply(lambda x: x[0])
    annAveragedMetrics['clippedWeightedStdCorr'] = annAveragedMetrics['clippedWeightedMean'].apply(lambda x: x[1])
    annAveragedMetrics['clippedWeightedMean'] = annAveragedMetrics['clippedWeightedMean'].apply(lambda x: x[0])

    annAveragedMetrics.to_csv(os.path.join(path, 'annualAvg.csv'))
    annAvgGroupedByObjid = annAveragedMetrics.groupby(level=0)
    aggregatedMetrics = annAvgGroupedByObjid.apply(computeMetricsAggYears)
    joined = appliedMetrics.join(aggregatedMetrics)
    #before printing split out the dicts:

    for columnToExpand in ['sigmaClipped', ]:
        keys = joined[columnToExpand].values[0].keys()
        for k in keys:
            joined['%s_%s'%(columnToExpand.replace('&',''),k)] = joined[columnToExpand].apply(lambda x: x[k])

    joined  = joined.drop('sigmaClipped', 1)

    joined.to_csv(os.path.join(path, 'metrics.csv'))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("inputDir",
                        help = "output directory")
    parser.add_argument("outputDir",
                        help = "output directory")
    parser.add_argument("filter",
                        help="u,g,r,i,z or all")
    parser.add_argument("firstPatch", type = int,
                        help = "starting patch in range. Looks for file <filter><firstPatch>_<lastPatch>.dat")
    parser.add_argument("lastPatch", type = int,
                        help = "starting patch in range. Looks for file <filter><firstPatch>_<lastPatch>.dat")


    args = parser.parse_args()

    if args.filter == 'all':
        for filterName in ['i', 'g', 'r', 'z', 'u']:
            run(args.inputDir, args.outputDir, filterName, args.firstPatch, args.lastPatch)
    else:
        run(args.inputDir, args.outputDir, args.filter, args.firstPatch, args.lastPatch)
