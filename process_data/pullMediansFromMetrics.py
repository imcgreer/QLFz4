from aggregates import *
import os
import pandas as pd
import numpy as np

path = '/lsst8/yusra/S13_IN2P3/agg/metricsFiles'
lc = {}
for f in 'ugriz':
    filename = os.path.join(path, '%s_metrics.csv'%(f))
    df = pd.DataFrame.from_csv(filename, index_col=(0,), parse_dates = False)
    columns = ['N', 'median', 'sigmaG', 'e_mean', 'e50', 'annAvg_chi2_clip_wm']
    lc[f] = df.reindex(columns=columns)
    lc[f]['median_mag']= flux2ab(lc[f]['median'])
    lc[f]['e50_mag'] = flux2absigma(lc[f]['median'], lc[f]['e50'])
    lc[f]['iqrSig_mag'] = flux2absigma(lc[f]['median'], lc[f]['sigmaG'])
    keys = lc[f].columns.values
    lc[f].columns = [s+'_%s'%(f) for s in keys]


df = lc['u'].join(lc['g']).join(lc['r']).join(lc['i']).join(lc['z'])
df.to_csv(os.path.join(path,'ugrizMetrics.csv'))

