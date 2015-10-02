"""File to extract lightcurves
"""
import os
import numpy as np
import pandas as pd
import argparse
from aggregates import *


def run(filename, root, outputRoot, filters):
    d = np.loadtxt(filename, delimiter=',', skiprows = 1,
                   dtype=[('id','i8'),('ra','f8'),('decl','f8')])
    classLabelSet = ['ALL']

    if outputRoot.endswith('/'):
        outputRoot = outputRoot[:-1]

    outputRootMag = outputRoot + '_mag'
    outputRootFlux = outputRoot + '_flux'
    for outputRoot in [outputRootMag, outputRootFlux]:
        for c in classLabelSet:
            path = os.path.join(outputRoot, c)
            if not os.path.exists(path):
                os.makedirs(path)
            for f in filters:
                filterPath = os.path.join(path,f)
                if not os.path.exists(filterPath):
                    os.makedirs(filterPath)

    patchRanges = ['365_386','344_365',
                     '323_344','302_323',
                     '281_302', '260_281',
                     '239_260', '218_239',
                     '197_218', '176_197',
                     '155_176']

    patchRanges = ['365_386','344_365',
                     '323_344','302_323',
                     '281_302', '260_281',
                     '239_260']

    for f in filters:
        print "    filter=%s of %s"%(f, filters)
        for patches in patchRanges:
            print "        patch range %s"%(patches)
            filename = '%s%s.csv'%(f, patches )
            df = pd.DataFrame.from_csv(os.path.join(root,filename), index_col=1)
            df['mag'] = flux2ab(df['psfFlux'])
            df['magErr'] = flux2absigma(df['psfFlux'],df['psfFluxErr'])
            objectIds = set(df.index)
            grouped = df.groupby(level=0)
            for classLabel in classLabelSet:
                print "doing class %s"%(classLabel)
                processHandleMag = open(os.path.join(outputRootMag, classLabel,'process_%s.dat'%(f)), 'a')
                processHandleFlux = open(os.path.join(outputRootFlux, classLabel,'process_%s.dat'%(f)), 'a')
                #join is very slow:
                #joinedGrouped = grouped.join(idf, how='inner')
                #    so do a seek instead:
                for deepSourceId, ra, decl in d:
                    #print type(long(deepSourceId))
                    #break
                    if deepSourceId in objectIds:
                        #print deepSourceId
                        try:
                            g = grouped.get_group(deepSourceId).sort(columns='exposure_time_mid')
                        except KeyError:
                            print 'keyError', deepSourceId
                            continue
                        lcHandleMag = open(os.path.join(outputRootMag, classLabel, f, '%s.dat'%(deepSourceId)), 'wb')
                        lcHandleFlux = open(os.path.join(outputRootFlux, classLabel, f, '%s.dat'%(deepSourceId)), 'wb')
                        N = len(g)
                        for row in zip(g['exposure_time_mid'].values, g['psfFlux'].values, g['psfFluxErr'].values,
                                       g['mag'].values, g['magErr'].values):
                            if np.isfinite(row[3]):
                                print >>lcHandleMag, row[0], row[3], row[4]
                            else:
                                N -= 1
                            print >>lcHandleFlux, row[0], row[1], row[2]
                        print >>processHandleMag, N, ra, decl, 0, '-25.00', "'%s'"%(os.path.join(f, '%s.dat'%(deepSourceId)))
                        print >>processHandleFlux, len(g), ra, decl, 0, '-25.00', "'%s'"%(os.path.join(f, '%s.dat'%(deepSourceId)))
                        lcHandleMag.close()
                        lcHandleFlux.close()
                processHandleMag.close()
                processHandleFlux.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename",
                        help = "filename containing id, z, class, ra, decl")
    parser.add_argument("inputDir",
                        help = "directory containing raw database dumps with format" +
                               "Looks for file <filter><firstPatch>_<lastPatch>.dat")
    parser.add_argument("outputDir",
                        help = "directory to extract lightcurves into")
    parser.add_argument("filter",
                        help="u, g, r,i,z or combination (e.g. 'gri'",
                        default = 'gri')
    args = parser.parse_args()
    run(args.filename, args.inputDir, args.outputDir, args.filter)


"""
note to self:
can delete nans later with:
sed -i -e '/nan/d' */*.dat
"""

