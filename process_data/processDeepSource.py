from __future__ import division
import pandas as pd
import numpy as np
import argparse
from columnNames import DEEPSOURCE_COLS, DEEPCOADD_COLS

# 2.5 / log(10)
FIVEOVER2LOG10  = 1.085736204758129569
# log(10) / 2.5
TWOLOG10_OVER_5 = 0.9210340371976182736
# 10^(-48.6/2.5)
AB_FLUX_SCALE = 3.630780547701013425e-20

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

COLS_TO_EXTRACT_ORIG  = [
   'parent', 'coadd_id', 'coadd_filter_id', 'coord_ra', 'coord_decl',
   'coord_raVar', 'coord_declVar', 'coord_radeclCov',
   'centroid_sdss_x', 'centroid_sdss_y', 'centroid_sdss_xVar', 'centroid_sdss_yVar', 'centroid_sdss_xyCov',
   'flux_psf', 'flux_psf_err', 'flux_sinc', 'flux_sinc_err', 'multishapelet_combo_flux',
   'multishapelet_combo_flux_err', 'flux_gaussian', 'flux_gaussian_err', 'correctfluxes_apcorr',
   'shape_sdss_centroid_x', 'shape_sdss_centroid_y', 'shape_sdss_centroid_xVar', 'shape_sdss_centroid_yVar',
   'shape_sdss_centroid_xyCov', 'shape_sdss_Ixx', 'shape_sdss_Iyy', 'shape_sdss_Ixy', 'shape_sdss_IxxVar',
   'shape_sdss_IyyVar', 'shape_sdss_IxyVar', 'shape_sdss_IxxIyyCov', 'shape_sdss_IxxIxyCov',
   'shape_sdss_IxxIxyCov', 'classification_extendedness', 'flags_negative', 'flags_badcentroid',
   'flags_pixel_edge', 'flags_pixel_interpolated_any', 'flags_pixel_interpolated_center',
   'flags_pixel_saturated_any', 'flags_pixel_saturated_center', 'flux_psf_flags', 'flux_sinc_flags',
   'multishapelet_combo_flux_flags', 'flux_gaussian_flags', 'centroid_sdss_flags', 'shape_sdss_flags', 'detect_is_primary']

COLS_TO_EXTRACT_RENAME = [
    'parentDeepSourceId', 'deepCoaddId', 'filterId', 'ra', 'decl',
    'raVar', 'declVar', 'radeclCov',
    'x', 'y', 'xVar', 'yVar', 'xyCov',
    'psfFlux', 'psfFluxSigma', 'apFlux', 'apFluxSigma', 'modelFlux',
    'modelFluxSigma', 'instFlux', 'instFluxSigma', 'apCorrection',
    'shapeIx', 'shapeIy', 'shapeIxVar', 'shapeIyVar',
    'shapeIxIyCov', 'shapeIxx', 'shapeIyy', 'shapeIxy', 'shapeIxxVar',
    'shapeIyyVar', 'shapeIxyVar', 'shapeIxxIyyCov', 'shapeIxxIxyCov',
    'shapeIyyIxyCov', 'extendedness', 'flagNegative', 'flagBadMeasCentroid',
    'flagPixEdge', 'flagPixInterpAny', 'flagPixInterpCen',
    'flagPixSaturAny', 'flagPixSaturCen', 'flagBadPsfFlux', 'flagBadApFlux',
    'flagBadModelFlux', 'flagBadInstFlux', 'flagBadCentroid', 'flagBadShape','detect_is_primary']

COLS_TO_EXTRACT_NARROW = [
    'parentDeepSourceId', 'deepCoaddId','ra','decl','psfMag','psfMagSigma',
    'tract', 'patch','detect_is_primary']


def run(filenameDeepSource, filenameDeepCoadd, filenameOut,
        chunksize=100000, magLimit=23.0, filterName='i'):
    """Join DeepSource with DeepCoadd files containing zeropoints

    Compute magnitudes and print a wide and narrow joined dataframe.
    Assumes that DeepCoadd file fits in memory.
    """
    compression = 'gzip' if filenameDeepCoadd.endswith('gz') else None
    dc = pd.read_csv(filenameDeepCoadd, compression=compression, header=None, index_col=0)
    dc.columns = DEEPCOADD_COLS
    dc.index.name = 'deepCoaddId'
    #Only need to join tract, patch, and fluxMag0:
    dc = dc[['tract', 'patch', 'fluxMag0', 'fluxMag0Sigma']]

    compression = 'gzip' if filenameDeepSource.endswith('gz') else None
    sources = pd.read_csv(filenameDeepSource,
                          compression=compression, header=None, chunksize=chunksize, index_col=0)
    for chunk in sources:
        chunk.columns  = DEEPSOURCE_COLS
        chunk.index.name = 'deepSourceId'
        #chunk = chunk[chunk.detect_is_primary == 1]
        ds = chunk[COLS_TO_EXTRACT_ORIG]
        ds.columns  = COLS_TO_EXTRACT_RENAME
        ds = ds.replace({'\N': ''})
        ds.index.names = ['deepSourceId']
        ds = ds.convert_objects(convert_numeric=True)
        joined = pd.merge(ds, dc,left_on='deepCoaddId',right_index=True)
        joined['psfFluxCal'] = dn2flux(joined.psfFlux, joined.fluxMag0)
        joined['psfFluxCalSigma'] = dn2fluxsigma(joined.psfFlux, joined.psfFluxSigma,
                                             joined.fluxMag0, joined.fluxMag0Sigma)
        joined['psfMag'] = flux2ab(joined['psfFluxCal'])
        joined['psfMagSigma'] =flux2absigma(joined['psfFluxCal'], joined['psfFluxCalSigma'])
        result = joined[joined['psfMag'] < magLimit]
        result.to_csv(filenameOut + '_%s_lt%i.csv'%(filterName, magLimit*10), mode='a')
        result[COLS_TO_EXTRACT_NARROW].to_csv(filenameOut + '_%s_lt%i_narrow.csv'%(filterName, magLimit*10),
                                              mode='a')

if __name__ == "__main__":
    usageExample = """python processDeepSource.py  
                      /lsst2/yusra/SDRP-IN2P3/coadd_dir/i/DeepSource.csv.gz
                      /lsst2/yusra/SDRP-IN2P3/coadd_dir/i/DeepCoadd.csv.gz
                      /lsst8/yusra/DeepSourceCsvs/DeepSource
                      i 100000 23.0
                   """
    parser = argparse.ArgumentParser(epilog=usageExample)
    parser.add_argument("DeepSourceFilename",
                        help = "Path to DeepSource.csv.gz")
    parser.add_argument("DeepCoaddFilename",
                        help = "Path to DeepCoadd.csv.gz")
    parser.add_argument("OutputFileRoot",
                        help = "Path to fileroot to print result csv." +
                        "For example: /path/to/DeepSourceJoined will print a " +
                        "DeepSourceJoined_lt_235.csv and DeepSourceJoined_lt_235_narrow.csv")
    parser.add_argument("filterName",
                        help="u, g, r, i, z or combination (e.g. 'gri'")
    parser.add_argument("chunksize", type=int,
                        help="number of rows of DeepSource file to process at a time."
                             "(must fit into into memory)",
                        default = 100000)
    parser.add_argument("magLimit", type=float,
                        help="Only print rows with PSF magnitudes less than this value.",
                        default = 23.5)

    args = parser.parse_args()
    print args
    for f in list(args.filterName):
        run(args.DeepSourceFilename, args.DeepCoaddFilename, args.OutputFileRoot,
            chunksize=args.chunksize, magLimit=args.magLimit, filterName=f)

