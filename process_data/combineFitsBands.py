from aggregates import *
import argparse
import os
import pandas as pd
import numpy as np

CLASS_DICT = {'GALAXY': 0, 'STAR':0, 'QSO':1 }
def getType(s):
    return CLASS_DICT[s]

def convertToIntClass(df):
    label = df['class']
    del df['class']
    label2 = label.apply(getType)
    return df

def run(path, outFilename='candidates_mag.csv', filters='gri', isTraining=False, keepEachFilter=True):
    classFits = {}
    classList = ["QSO","STAR","GALAXY"] if isTraining else ['ALL',]
    for classLabel in classList:
        fits = {}
        for f in filters:
            fits[f] = pd.DataFrame.from_csv(os.path.join(path, classLabel, 'fits_%s.dat'%(f)),
                                              index_col=(1), parse_dates = False)
            keys = fits[f].columns.values
            fits[f].columns = [s+'_%s'%(f) for s in keys]
        df = fits['g'].join(fits['r'], how='outer').join(fits['i'], how='outer')
        df['class'] = classLabel
        classFits[classLabel] = df

    if isTraining:
        full_df = classFits['QSO'].append(classFits['STAR']).append(classFits['GALAXY'])
    else:
        full_df = classFits['ALL']

    #calc means:
    for key in ['log_geo_mean_tau', 'log_arith_mean_tau', 'std_tau',
                'log_geo_mean_sigma' , 'log_arith_mean_sigma', 'std_sigma',
                'meanPlikePnoise', 'meanPlikePinf']:
        full_df[key+'_gri'] =  full_df[[key+'_g',
                                        key+'_r',
                                        key+'_i']].mean(axis=1)
    #calc sums
    for key in ['log_geo_mean_tau',
                'log_geo_mean_sigma',
                'meanPlikePnoise', 'meanPlikePinf']:
        full_df['sum_' + key+'_gri'] =  full_df[[key+'_g',
                                        key+'_r',
                                        key+'_i']].sum(axis=1)

    full_df['count_gri'] =  full_df[[key+'_g',
                                     key+'_r',
                                     key+'_i']].count(axis=1)

    if not keepEachFilter:
        for key in full_df.columns.values:
                if (key.endswith('_i')) |  (key.endswith('_g')) | (key.endswith('_r')):
                    print key
                    del full_df[key]

    full_df.to_csv(os.path.join(path, outFilename))
    return full_df

if __name__ == "__main__":
    """Example:
       python combineFitsBands.py /lsst8/yusra/S13_IN2P3/candidates_bad_mag \
       gri candidates_mag_test.csv --keepEachFilter
       """
    parser = argparse.ArgumentParser()
    parser.add_argument("path",
                        help = "input/output directory")
    parser.add_argument("filters",
                        help="gri by default")
    parser.add_argument("output_filename",
                        help="candidates_mag.csv")
    parser.add_argument("--isTraining", dest='isTraining', action='store_true',
                        help="path has three subdirectories: QSO, GALAXY, STAR")
    parser.add_argument("--keepEachFilter", dest='keepEachFilter', action='store_true',
                        help="Print individual g, r, i columns in output table too?")

    args = parser.parse_args()
    run(args.path, args.output_filename, args.filters, args.isTraining, args.keepEachFilter)
