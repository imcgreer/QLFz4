from __future__ import division
import pandas as pd
import numpy as np
import argparse
import os
from columnNames import SCIENCE_CCD_EXPOSURE_COLS, FORCEDSOURCE_COLS

# 2.5 / log(10)
FIVEOVER2LOG10  = 1.085736204758129569
# log(10) / 2.5
TWOLOG10_OVER_5 = 0.9210340371976182736
# 10^(-48.6/2.5)
AB_FLUX_SCALE = 3.630780547701013425e-20

FS_COLS_TO_WRITE = ['objectId', 'exposure_id','exposure_time_mid','psfFlux','psfFluxErr']

def hypot(a, b):
    if (a == 0.0) | (b == 0.0):
        return 0.0
    elif np.abs(a) < np.abs(b):
        c = a
        a = b
        b = c
    q = b / a
    return np.abs(a) * np.sqrt(1.0 + q*q)

vhypot = np.vectorize(hypot)

def flux2ab(flux):
    return -2.5 * np.log10(flux) - 48.6

def flux2absigma(flux, fluxsigma):
    return FIVEOVER2LOG10 * fluxsigma / flux

def dn2fluxsigma(dn, dnsigma, fluxmag0, fluxmag0sigma):
    return AB_FLUX_SCALE * vhypot(dn * fluxmag0sigma, dnsigma * fluxmag0)/(fluxmag0*fluxmag0)

def dn2flux(dn, fluxmag0):
    return AB_FLUX_SCALE * dn / fluxmag0

def makePatchRangeString(x, bins):
    x = np.array(x).astype(int)
    inds = np.digitize(x, bins)
    inds[np.where(inds == len(bins))] = len(bins) - 1
    return ['%s_%s'%(bins[i-1], bins[i]) for i in inds]


def write2Partition(df, outDir, filterName, groupByCol='patchRange', colsToWrite= FS_COLS_TO_WRITE):
    groups = df.groupby('patchRange')
    for group in groups:
        #group[1].sort(columns=['objectId', 'exposure_id'], inplace=True)
        group[1][colsToWrite].to_csv(os.path.join(outDir, '%s%s.csv'%(filterName, group[0])), mode='a')




def run(indir, scienceCcdExpFilename, deepSourceFilename,
        outDir, filterName, bins, chunksize):
    sce = pd.read_csv(scienceCcdExpFilename,
                     header=None, index_col=0)
    sce.columns = SCIENCE_CCD_EXPOSURE_COLS
    sce.index.name = 'exposure_id'
    ds = pd.DataFrame.from_csv(deepSourceFilename, index_col=0)
    ds.index = ds.index.astype(int)
    for run in sorted(os.listdir(indir)):
        if run.isdigit():
            path = os.path.join(indir, run, 'DeepForcedSource.csv.gz')
            print path
            compression = 'gzip' if path.endswith('gz') else None
            try:
                sources = pd.read_csv(path, compression=compression, header=None,
                                      chunksize=chunksize, index_col=0)
            except Exception as e:
                print "WARNING: %s: %s"%(path, e)
                continue
            for f in sources:
                print 'reading next %s'%(chunksize)
                f.columns = FORCEDSOURCE_COLS
                f.index.name = 'deepForcedSourceId'
                #filter on flags
                f = f[~(f.flags_pixel_edge | f.flags_pixel_saturated_any | f.flags_pixel_saturated_center |
                      f.flux_sinc_flags | f.flux_psf_flags_psffactor | f.flags_badcentroid)]
                f = f[['objectId', 'exposure_id','exposure_time_mid','flux_psf', 'flux_psf_err']]
                try:
                    f = f.replace({'\N': ''})
                except:
                    pass
                f = f.convert_objects(convert_numeric=True)
                #only take i < 23
                joined = pd.merge(f, ds, left_on='objectId', right_index=True, how='inner')
                #get fluxmag0s
                f2 = pd.merge(joined, sce[['fluxMag0', 'fluxMag0Sigma']], left_on='exposure_id',
                              right_index=True, how='inner')
                f2['psfFlux'] = dn2flux(f2.flux_psf, f2.fluxMag0)
                f2['psfFluxErr'] = dn2fluxsigma(f2.flux_psf, f2.flux_psf_err,
                                                f2.fluxMag0, f2.fluxMag0Sigma)
                f2['xPatch'] = f2.patch.apply(lambda x: x.split(',')[0])
                try:
                    f2['patchRange'] = makePatchRangeString(f2.xPatch, bins)
                except:
                    import pdb; pdb.set_trace()
                write2Partition(f2, outDir, filterName, groupByCol='patchRange', colsToWrite= FS_COLS_TO_WRITE)

if __name__ == "__main__":
    usageExample = """python processForcedSource.py
                      /lsst2/yusra/SDRP-IN2P3/forcedPhot/i
                      /lsst2/yusra/SDRP-IN2P3/ingestProcessed_csv/i_Science_Ccd_Exposure.csv
                      /lsst8/yusra/DeepSourceCsvs/DeepSource_i_lt230_narrow.csv
                     /lsst8/yusra/S13_IN2P3/forcedPhot/ i 100000
]$ python processForcedSource.py /lsst2/daues/forcedPhot_csv_dir/z /lsst2/yusra/SDRP-IN2P3/ingestProcessed_csv/z_Science_Ccd_Exposure.csv /lsst8/yusra/DeepSourceCsvs/DeepSource_i_lt230_narrow.csv /lsst8/yusra/S13_IN2P3/forcedPhot_all z 100000
                   """
    parser = argparse.ArgumentParser(epilog=usageExample)
    parser.add_argument("forcedSourceFolder",
                        help = """Path to folder containing DeepForcedSource Files
                                   Assumes that The folder has the following structure:
                                   ForcedSourceFolder:
                                        1755/
                                           DeepForcedSource.csv.gz
                                        4147/
                                           DeepForcedSource.csv.gzip
                                """)
    parser.add_argument("scienceCcdExposureFilename",
                        help = "Path to Science_Ccd_Exposure.csv file.")
    parser.add_argument("deepSourceFilename",
                        help = "Path to DeepSource file. Only rows that match these objects " +
                               "will be output")
    parser.add_argument("OutputFolder",
                        help = "Path to folder to print resulting csvs")
    parser.add_argument("filterName",
                        help="u, g, r, i, z or combination (e.g. 'gri'",
                        default = 'i')
    parser.add_argument("chunksize", type=int,
                        help="number of rows of DeepSource file to process at a time."
                             "(must fit into into memory)",
                        default = 100000)
    args = parser.parse_args()
    print args
    bins = range(386,159-21,-21)[::-1]

    for f in args.filterName:
        run(args.forcedSourceFolder, args.scienceCcdExposureFilename, args.deepSourceFilename,
            args.OutputFolder, filterName=f, bins=bins, chunksize=args.chunksize)
    print "exiting normally"


