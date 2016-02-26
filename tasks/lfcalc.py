#!/usr/bin/env python

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from simqso import sqanalysis,lumfun
import qlffit

class mag2Lum(object):
    def __init__(self,simName,simDir,cosmo):
        self.lcorr = sqanalysis.calcKCorrFromGrid(simName+'_Lum',
                                                  simDir,bandNum=3)
        self.fcorr = sqanalysis.calcKCorrFromGrid(simName+'_Flux',
                                                  simDir,bandNum=3)
        self.cosmo = cosmo
    def __call__(self,m,z,inverse=False):
        if inverse:
            return self.lcorr(m,z) + self.cosmo.distmod(z).value
        else:
            return self.fcorr(m,z) + self.cosmo.distmod(z).value

def qlfz4_selection_function(mag,magErr,flux,fluxErr):
    # no selection criteria - just return True for all
    return (mag[:,3] < 22)

import pandas as pd

def loadSelectionGrid(m_arr, z_arr, fraction, binsM, binsZ):
    """Take luminosity and redshift arrays and return selection
       function array
    """
    import pdb; pdb.set_trace()
    df = pd.DataFrame()
    df['M'] = m_arr
    df['z'] = z_arr
    df['Mbin'] = pd.cut(df.M, binsM)
    df['zbin'] = pd.cut(df.z, binsZ)
    fractionCutDF = fraction #.reset_index()
    fractionCutDF.columns = ['M', 'z', 'fraction']
    joined = pd.merge(df, fractionCutDF, how='left', right_on=('M', 'z'),
                      left_on=('Mbin', 'zbin'))['fraction']
    return joined.fillna(np.inf).values

class SelectionFunction(object):
    def __init__(self,simName,simDir,photo_complete,spec_complete,m2M,
                 sfx='',**kwargs):
        self.m2M = m2M
        if False:
            self.selcrit = qlfz4_selection_function
            self.M_selfun = sqanalysis.calcSelectionFunctionFromGrid(
                                                        simName+'_Lum',
                                                        self.selcrit,simDir,
                                                        **kwargs)
            self.m_selfun = sqanalysis.calcSelectionFunctionFromGrid(
                                                        simName+'_Flux',
                                                        self.selcrit,simDir,
                                                        **kwargs)
        else:
            self.binsM = np.arange(-27.5, -23.5, 0.2)
            self.binsZ = np.arange(3.25, 4.75, 0.1)
            self.fractionDF = pd.read_csv('/Users/yusra/Dropbox/fractionNov30.csv')
            self.m_selfun = lambda m,z: loadSelectionGrid(m - self.m2M(m, z),z,
                                        self.fractionDF,self.binsM,self.binsZ)
            self.M_selfun = lambda M,z: loadSelectionGrid(M, z,
                                         self.fractionDF,self.binsM,self.binsZ)
        self.photo_complete = photo_complete
        self.spec_complete = spec_complete
    def _photo_complete_call(self,m,z,absMag=False):
        if type(self.photo_complete) is float:
            return self.photo_complete
        elif absMag:
            return self.photo_complete(m-self.m2M(m,z,inverse=True),z)
        else:
            return self.photo_complete(m,z)
    def _spec_complete_call(self,m,z,absMag=False):
        if type(self.spec_complete) is float:
            return self.spec_complete
        elif absMag:
            return self.spec_complete(m-self.m2M(m,z,inverse=True),z)
        else:
            return self.spec_complete(m,z)
    def _color_complete_call(self,m,z,absMag=False):
        if absMag:
            #return self.selfun(m,z)
            return self.M_selfun(m,z)
        else:
            #return self.selfun(m-self.m2M(m,z,inverse=True),z)
            return self.m_selfun(m,z)
    def __call__(self,m,z,absMag=False):
        import pdb; pdb.set_trace()
        return self._color_complete_call(m,z,absMag) * \
               self._photo_complete_call(m,z,absMag) * \
               self._spec_complete_call(m,z,absMag)

from scipy.integrate import quad

def get_dVdzdO(cosmo):
    _zfun = lumfun.interp_dVdzdO((3.0,5.0),cosmo)
    dVdzdO = lambda z1,z2: quad(_zfun,z1,z2)[0]
    return dVdzdO

def placeholder_survey():
i_lim = 22.0 #for binned lf
surveyArea = 54.6
simName = 'SDSS'
simDir = '/Users/yusra/'
cosmo = FlatLambdaCDM(H0=70,Om0=0.3)
m2M = mag2Lum(simName,simDir,cosmo)
Nsurvey = pd.read_csv('/Users/yusra/Dropbox/placeholder_north_survey.csv')
dat = Table.read('/Users/yusra/Dropbox/placeholder_north_survey.csv')
datUpdatedRedshifts = pd.read_csv('/Users/yusra/calcZ/goodRedshiftsErrors_redo', sep=' ', header=None)
datUpdatedRedshifts.columns = ('name', 'z+1', 'err')
datUpdatedRedshifts = Table.read('/Users/yusra/Dropbox/placeholder_north_survey.csv')
survey = qlffit.QuasarSurvey(dat['median_mag_i']-2*2.086*dat['ebv'],
                             dat['z'],i_lim,surveyArea,m2M)
photo_complete = 1.0
spec_complete = 1.0
selfun = SelectionFunction(simName,simDir,photo_complete,spec_complete,m2M)
survey.set_selection_function(selfun)
Medges = np.arange(-27.5,-23.0,0.5)
zedges = np.arange(3.75,4.25,0.5)

lfcomp = survey.calcBinnedLF(Medges,zedges,get_dVdzdO(cosmo))
#return lf


#>>> zcenters = 0.5*(zedges[:-1] + zedges[1:])
#>>> zcenters
#array([ 3.35,  3.55,  3.75,  3.95,  4.15,  4.35])
#>>> Mcenters = 0.5*(Medges[:-1] + Medges[1:])

ax = plt.subplot(111)
plt.imshow(lf['counts']/lf['rawCounts'],interpolation='nearest', 
           extent=[zedges[0], zedges[-1], Medges[0],Medges[-1]], origin='lower',
           aspect='auto', cmap='hot')
plt.colorbar()
plt.show()
#ax.set_xticks(zcenters)
#ax.set_yticks(Mcenters)
#ax.set_yticklabels(['%s'%t for t in zcenters])
#ax.set_xticklabels(['%s'%t for t in Mcenters])
#plt.show()
