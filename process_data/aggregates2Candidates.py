import pandas as pd
import numpy as np

SIDE = 'IN2P3'
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

if SIDE == 'NCSA':
    DEEPSOURCE_NAME = '/lsst8/yusra/DeepSourceCsvs/DeepSourceNCSA_i_lt235.csv.gz'
    PARENTDEEPSOURCE_NAME = '/lsst8/yusra/DeepSourceCsvs/DeepSourceNCSA_i_lt300_narrow.csv.gz'
    UGRIZ_NAME = '/lsst8/yusra/S13Aggregates/lightcurveMetrics/metricsFiles/ugrizMetrics.csv'
    GRI_NAME = '/lsst8/yusra/S13Aggregates/lightcurveMetrics/gri/gri_chi2.csv'
    EBV_NAME = '/lsst8/yusra/DeepSourceCsvs/ebv_NCSA_lt235.dat'
else:  # SIDE == 'IN2P3'
    DEEPSOURCE_NAME = '/lsst8/yusra/DeepSourceCsvs/DeepSource_i_lt230.csv'
    PARENTDEEPSOURCE_NAME = '/lsst8/yusra/DeepSourceCsvs/DeepSourceIN2P3_i_lt235_narrow_not_primary.csv'
    UGRIZ_NAME = '/lsst8/yusra/S13_IN2P3/agg/metricsFiles/ugrizMetrics.csv'
    GRI_NAME = '/lsst8/yusra/S13_IN2P3/agg/gri/griChi2.csv'
    EBV_NAME = '/lsst8/yusra/DeepSourceCsvs/ebv_ds_lt230'


def getCompress(filename):
    return 'gzip' if filename.endswith('gz') else None

df = pd.read_csv(DEEPSOURCE_NAME, compression=getCompress(DEEPSOURCE_NAME), index_col=0)
agg = pd.read_csv(UGRIZ_NAME, compression=getCompress(UGRIZ_NAME), index_col=0)
gri = pd.read_csv(GRI_NAME, compression=getCompress(GRI_NAME), index_col=0)
ebv = pd.read_csv(EBV_NAME, compression=getCompress(EBV_NAME),  sep=' ', index_col=0, header=None,
                  comment='#', names=['ebv', ])
parent = pd.read_csv(PARENTDEEPSOURCE_NAME,  compression=getCompress(PARENTDEEPSOURCE_NAME),
                     index_col=0, dtype=TYPE,)

dfNew = pd.merge(df, agg[agg.columns.difference(df.columns)],
                 left_index=True, right_index=True, how='inner')
dfNew = pd.merge(dfNew, gri[gri.columns.difference(dfNew.columns)],
                 left_index=True, right_index=True, how='inner')
d = pd.merge(dfNew, ebv[ebv.columns.difference(dfNew.columns)],
             left_index=True, right_index=True, how='inner')


A = np.array([5.155, 3.793, 2.751, 2.086, 1.479])
Adiff = -np.diff(A)

# deredden
for i, filt in enumerate('ugriz'):
    d['median_mag_%s' % (filt)] = d['median_mag_%s' % (filt)] - d['ebv']*A[i]

# calc colors
for blue, red in zip('ugri', 'griz'):
    d['%s%s' % (blue, red)] = d['median_mag_%s' % (blue)] - d['median_mag_%s' % (red)]


#If modelFlux failed, fall back on instFlux
d['extend'] = 2.5*np.log10(d.modelFlux/d.psfFlux)
series =  2.5*np.log10(d.instFlux/d.psfFlux)
d['extend'].loc[d['extend'].isnull()] = series.loc[d['extend'].isnull()]


# (1) Cut on survey area. RA boundaries and "holes" around bright stars:

parent = parent.loc[parent.psfMag < 17.0]
parent_df = pd.merge(d, parent[['psfMag', 'psfMagSigma', 'ra', 'decl']],
                     how='left', suffixes=['', '_parent'],
                     left_on='parentDeepSourceId', right_index=True)

#NCSA 5446502 -> 4867183
#IN2P3 3427395 -> 3130661
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
dview = dview.loc[dview.extend < 0.1]
# NCSA -> 29547
# IN2P3 ->17304
dview = dview.loc[(dview.annAvg_chi2_clip_wm_r > 1) & (dview.chi2griClippedWMean > 2)]
# NCSA -> 12855
# IN2P3 ->7049
dview.to_csv('candidates_%s.csv' % (SIDE))



"""
(median_mag_i > 17) & (median_mag_i < 23.0) & (gr > 0.95) & (ri < 0.44*gr -0.08) & (ri < 0.3*gr +0.03)
(extendedness < 0.1)
(iz < 1 & iz > -0.7) &  (median_mag_u > 22.5)
(ri < 0.8 & ri > -1)
annAvg_chi2_clip_wm > 1 & chi2griClippedWMean > 2


NCSA missing from 2013:
Int64Index([  3509645813486539, 216577406016689313, 216814908044478082,
            216876476400667004, 216982028443193770, 217360259369404035,
            217474605357470359, 217518584748843759, 217562573803886081,
            217580163842448741, 217632935031871530, 217632935031871532],
           dtype='int64', name=u'deepSourceId')
"""


