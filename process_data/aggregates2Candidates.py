import pandas as pd
import numpy as np
import argparse

TYPE = {'deepSourceId': 'int64',
        'parentDeepSourceId': 'int64',
        'deepCoaddId': int,
        'ra': float,
        'decl': float,
        'psfMag': float,
        'psfMagSigma': float,
        'tract': int,
        'patch': str,
        'detect_is_primary': int
        }


def transformColumns(d):
    """Calculate computed columns
    """
    A = np.array([5.155, 3.793, 2.751, 2.086, 1.479])
    # deredden
    for i, filt in enumerate('ugriz'):
        d['median_mag_' + filt] = d['median_mag_' + filt] - d['ebv']*A[i]
        d['median_mag_' + filt].loc[d['median_mag_' + filt].isnull()] = \
            d['median_mag_' + filt].max()
    # calc colors
    for blue, red in zip('ugri', 'griz'):
        d['%s%s' % (blue, red)] = d['median_mag_%s' % (blue)] - d['median_mag_%s' % (red)]
    # If modelFlux failed, fall back on instFlux
    d['extendedness'] = 2.5*np.log10(d.modelFlux/d.psfFlux)
    series = 2.5*np.log10(d.instFlux/d.psfFlux)
    d['extendedness'].loc[d['extendedness'].isnull()] = series.loc[d['extendedness'].isnull()]
    # instrinsic sigma
    for f in 'ugriz':
        d['sigma_intrinsic_' + f] = np.sqrt(d['iqrSig_mag_' + f]**2 - d['e50_mag_' + f]**2)
        d['sigma_intrinsic_' + f].loc[d['sigma_intrinsic_' + f].isnull()] = 0
    d['patchx'] = d['patch'].apply(lambda x: x.split(',')[0]).astype(int)
    d['patchy'] = d['patch'].apply(lambda x: x.split(',')[1]).astype(int)
    return d


def getCompress(filename):
    return 'gzip' if filename.endswith('gz') else None


def run(SIDE):
    if SIDE == 'NCSA':
        DEEPSOURCE_NAME = '/lsst8/yusra/DeepSourceCsvs/DeepSourceNCSA_i_lt235.csv.gz'
        PARENTDEEPSOURCE_NAME = '/lsst8/yusra/DeepSourceCsvs/DeepSourceNCSA_i_lt300_narrow.csv.gz'
        UGRIZ_NAME = '/lsst8/yusra/S13Aggregates/lightcurveMetrics/metricsFiles/ugrizMetrics.csv'
        GRI_NAME = '/lsst8/yusra/S13Aggregates/lightcurveMetrics/gri/gri_chi2.csv'
        EBV_NAME = '/lsst8/yusra/DeepSourceCsvs/ebv_NCSA_lt235.dat'
        R_METRICS = '/lsst8/yusra/S13Aggregates/lightcurveMetrics/metricsFiles/r_metrics.csv'
    else:  # SIDE == 'IN2P3'
        DEEPSOURCE_NAME = '/lsst8/yusra/DeepSourceCsvs/DeepSource_i_lt230.csv'
        PARENTDEEPSOURCE_NAME = '/lsst8/yusra/DeepSourceCsvs/DeepSourceIN2P3_i_lt235_narrow_not_primary.csv'
        UGRIZ_NAME = '/lsst8/yusra/S13_IN2P3/agg/metricsFiles/ugrizMetrics.csv'
        GRI_NAME = '/lsst8/yusra/S13_IN2P3/agg/gri/griChi2.csv'
        EBV_NAME = '/lsst8/yusra/DeepSourceCsvs/ebv_ds_lt230'
        R_METRICS = '/lsst8/yusra/S13_IN2P3/agg/metricsFiles/r_metrics.csv'

    df = pd.read_csv(DEEPSOURCE_NAME, compression=getCompress(DEEPSOURCE_NAME), index_col=0)
    agg = pd.read_csv(UGRIZ_NAME, compression=getCompress(UGRIZ_NAME), index_col=0)
    gri = pd.read_csv(GRI_NAME, compression=getCompress(GRI_NAME), index_col=0)
    ebv = pd.read_csv(EBV_NAME, compression=getCompress(EBV_NAME),  sep=' ', index_col=0, header=None,
                      comment='#', names=['ebv', ])
    parent = pd.read_csv(PARENTDEEPSOURCE_NAME,  compression=getCompress(PARENTDEEPSOURCE_NAME),
                         index_col=0, dtype=TYPE,)
    r = pd.read_csv(R_METRICS,  compression=getCompress(R_METRICS),
                    index_col=0)
    r.columns = [rcol + '_r' for rcol in r.columns]

    dfNew = pd.merge(df, agg[agg.columns.difference(df.columns)],
                     left_index=True, right_index=True, how='inner')
    dfNew = pd.merge(dfNew, gri[gri.columns.difference(dfNew.columns)],
                     left_index=True, right_index=True, how='inner')
    d = pd.merge(dfNew, ebv[ebv.columns.difference(dfNew.columns)],
                 left_index=True, right_index=True, how='inner')
    if False:
        cols = r.columns.difference(dfNew.columns)
    else:
        cols = ['sigmaClipped_chi2_r', ]
    d = pd.merge(d, r[cols], left_index=True, right_index=True, how='inner')

    d = transformColumns(d)

    # (1) Cut on survey area. RA boundaries and "holes" around bright stars:
    parent = parent.loc[parent.psfMag < 17.0]
    parent_df = pd.merge(d, parent[['psfMag', 'psfMagSigma', 'ra', 'decl']],
                         how='left', suffixes=['', '_parent'],
                         left_on='parentDeepSourceId', right_index=True)

    # NCSA 5446502 -> 4867183
    # IN2P3 3427395 -> 3130661
    d = parent_df.loc[parent_df.psfMag_parent.isnull()]

    if SIDE == 'NCSA':
        # NCSA 4867183 -> 4666467
        d = d.loc[(d.ra > 321.5) | (d.ra < 10)]

    dview = d.loc[(d.median_mag_i > 17) &
                  (d.median_mag_i < 23.0) &
                  (d.gr > 0.95) &
                  (d.ri < 0.44*d.gr - 0.08) &
                  (d.ri < 0.3*d.gr + 0.03) &
                  (d.iz < 1) &
                  (d.iz > -0.7) &
                  (d.ri > -1) &
                  (d.ri < 0.8)]

    # NCSA ->  110072
    # IN2P3 -> 92342

    dview = dview.loc[(dview.median_mag_u > 22.5) | dview.median_mag_u.isnull()]
    # NCSA -> 108818
    # IN2P3 ->91010
    dview = dview.loc[dview.extendedness < 0.1]
    # NCSA -> 29547
    # IN2P3 ->17304
    dview = dview.loc[(dview.annAvg_chi2_clip_wm_r > 1) &
                      (dview.chi2griClippedWMean > 2)]
    # NCSA -> 12855
    # IN2P3 ->7049

    # transform chi2 to log and fill in nans:
    # To do: change the column names to log chi2
    for col in dview.columns.values:
        # Was in original version but is a no-op if taking log
        # if col.startswith('chi2gri'):
        #    dview[col].loc[dview[col].isnull()] = 0
        if ('chi2' in col):
            dview[col] = np.log10(dview[col])

    dview.to_csv('candidates_%s.csv' % (SIDE))
    return dview


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("side",
                        help="NCSA, IN2P3, BOTH",
                        default="BOTH")

    args = parser.parse_args()
    if args.side == "BOTH":
        ncsa = run("NCSA")
        in2p3 = run("IN2P3")
        # 365 | 366 is the boundary ~ 10deg (use ncsa version of overlap region)
        both = pd.concat([ncsa[ncsa['patchx'] != 365],
                          in2p3[in2p3['patchx'] <= 365]])
        both.to_csv('candidates_merged.csv')
    else:
        run(args.side)
"""
(median_mag_i > 17) & (median_mag_i < 23.0) & (gr > 0.95) & (ri < 0.44*gr -0.08) & (ri < 0.3*gr +0.03)
(extendedness < 0.1)
(iz < 1 & iz > -0.7) &  (median_mag_u > 22.5)
(ri < 0.8 & ri > -1)
annAvg_chi2_clip_wm > 1 & chi2griClippedWMean > 2
"""


