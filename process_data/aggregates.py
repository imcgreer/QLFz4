import os, sys
import numpy as np
import itertools

C = 2.99792458e14 #um/s
FIVE_OVER_2LOG10 = 1.085736204758129569
TWOLOG10_OVER_5 = 0.9210340371976182736
AB_FLUX_SCALE = 3.630780547701013425e-20

def ab2flux(mag):
    """Compute flux given AB magnitude
    """
    return np.power(10.0, -0.4*(mag + 48.6))

def ab2fluxsigma(mag, magsigma):
    """compute flux sigma given AB mag and AB mag Sigma"""
    mag = np.array(mag)
    magsigma = np.array(magsigma)
    return magsigma * ab2flux(mag) * TWOLOG10_OVER_5;


def flux2absigma(flux, fluxsigma):
    """Compute AB mag sigma given flux and flux sigma"""
    return FIVE_OVER_2LOG10 * fluxsigma / flux;


def flux2ab(flux):
    """Compute AB mag given flux"""
    return -2.5 * np.log10(flux) - 48.6;

def calcWM(y, yerr):
    """Calculate the mean weighted by inverse variance
    """
    N = len(y)
    if N == 1:
        return float(y)
    elif N == 0:
        return np.nan
    else:
        return float(np.add.reduce(y/(yerr*yerr)) / np.add.reduce((1/yerr)*(1/yerr)))

def calcWStd(yerr):
    """Calculate the standard deviation weighted by inverse variance.
    Does not take into account the scatter of the points themselves
    """
    return np.sqrt(1. / np.add.reduce((1/yerr)*(1/yerr)))

def calcWVar(yerr):
    """Calculate the standard deviation weighted by inverse variance.
    Does not take into account the scatter of the points themselves
    """
    return 1. / np.add.reduce((1./(yerr*yerr)))

def calcWStdCorr(y, yerr):
    """Calculate just the standard deviation weighted by inverse variance
    and corrected for the intrinsic scatter.  Weighted mean is an intermediate 
    step and gets calculated anyway and thrown away
    """
    wm, stdCorr = calcWStdCorrAndMean(y, yerr)
    return float(stdCorr)

def calcWStdCorrAndMean(y, yerr):
    """Calculate the weighted mean and standard deviation weighted by inverse variance
    and corrected for the intrinsic scatter
    """
    N = len(y)
    if N  == 1:
        return float(y), float(yerr)
    elif N == 0:
        return np.nan, np.nan
    else:
        wm = calcWM(y, yerr)
        return wm, np.sqrt(calcWVar(yerr) * (1./(N-1)) * \
                       np.add.reduce((y - wm)*(y - wm)/(yerr*yerr)))

def sigmaClippedWM(y, yerr, sigmathresh):
    N = len(y)
    if N ==0:
        return np.nan
    elif N == 1:
        return float(y)
    else:
        q25, q50, q75 = np.percentile(y, (25, 50, 75))
        iqrSig = 0.7413*(q75 - q25)
        nsigma = np.abs(y - q50)/iqrSig
        idx = np.where(nsigma <= sigmathresh)
        return calcWM(y[idx], yerr[idx])

def sigmaClippedWStdCorr(y, yerr, sigmathresh):
    N = len(y)
    if N ==0:
        return np.nan
    elif N == 1:
        return float(yerr)
    else:
        q25, q50, q75 = np.percentile(y, (25, 50, 75))
        iqrSig = 0.7413*(q75 - q25)
        nsigma = np.abs(y - q50)/iqrSig
        idx = np.where(nsigma <= sigmathresh)
        wm, stdCorr = calcWStdCorrAndMean(y[idx], yerr[idx])
        return float(stdCorr)

def sigmaClippedWStdCorrAndMean(y, yerr, sigmathresh):
    N = len(y)
    if N ==0:
        return np.nan
    elif N == 1:
        return float(y), float(yerr)
    else:
        q25, q50, q75 = np.percentile(y, (25, 50, 75))
        iqrSig = 0.7413*(q75 - q25)
        nsigma = np.abs(y - q50)/iqrSig
        idx = np.where(nsigma <= sigmathresh)
        wm, stdCorr = calcWStdCorrAndMean(y[idx], yerr[idx])
        return float(wm), float(stdCorr)

def calcWMmag(mag, magsigma):
    """Compute weighted mean of magnitudes"""
    y =  ab2flux(mag)
    yerr = ab2fluxsigma(mag, magsigma)
    wmflux = calcWM(y, yerr)
    return flux2ab(wmflux)

def calcWMStdCorrMag(mag, magsigma):
    """Compute magnitude weighted mean and corrected standard deviation"""
    y =  ab2flux(mag)
    yerr = ab2fluxsigma(mag, magsigma)
    wm, stdCorr = calcWStdCorrAndMean(y, yerr)
    return flux2ab(wm), flux2absigma(wm, stdCorr)

def calcWMStdMag(mag, magsigma):
    """Compute magnitude weighted mean and standard deviation"""
    y =  ab2flux(mag)
    yerr = ab2fluxsigma(mag, magsigma)
    wm = calcWM(y, yerr)
    std = np.sqrt(calcWStd(yerr))
    return flux2ab(wm), flux2absigma(wm, std)

def calcWChi2(y, yerr):
    N = len(y)
    if N < 2:
        return np.nan
    else:
        wm = calcWM(y, yerr)
        chi2 = np.add.reduce(((y-wm)/yerr)**2)
        return chi2/(N-1)

def calcChi2raw(y, yerr):
    """Compute simple reduced chi2 if more than 1 datapoints present
    """
    N = len(y)
    if N < 2:
        return np.nan
    else:
        chi2 = np.sum(((y-np.mean(y))/yerr)**2)
        return chi2/(N-1)


def calcChi2(y, yerr):
    N = len(y)
    if N < 2:
        return np.nan
    else:
        chi2 = np.add.reduce(((y-np.mean(y))/yerr)**2)
        return chi2/(N-1)

def clip(y, yerr, sigmaThresh):
    """Return flux array with outliers clipped"""
    q25, q50, q75 = np.percentile(y, (25, 50, 75))
    iqrSig = 0.7413*(q75 - q25)
    nsigma = np.abs(y - q50)/iqrSig
    idx = np.where(nsigma <= sigmaThresh)
    return y[idx], yerr[idx]

"""
def combineEpochs(t, m, merr, bins=np.arange(50925, 55305, 365)):
    "Average magnitudes in light curve in provided bins, length of time
       or stripe 82 season"
    groups = np.digitize(t, bins)
    tuples = zip(groups, t, m, merr)
    sorted_tuples = sorted(tuples)
    j = 0
    newtime = []
    newmag = []
    newerr = []
    for j, (key, group) in enumerate(itertools.groupby(sorted_tuples, lambda x: x[0])):
        tups = [g for g in group]
        if len(tups) > 1:
            m, merr = calcWMStdCorrMag([t[2] for t in tups], [t[3] for t in tups])
            meanTime = np.mean([t[1] for t in tups])
            newtime.append(meanTime)
            newmag.append(m)
            newerr.append(merr)
        else:
           newtime.append(tups[0][1])
           newmag.append(tups[0][2])
           newerr.append(tups[0][3])

"""
