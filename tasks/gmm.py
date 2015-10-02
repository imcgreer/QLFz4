import copy
import numpy as np
from sklearn.mixture import GMM as skGMM
from scipy.stats import multivariate_normal
DEBUG = False

def gmmFromParams(means, covars, weights):
    nComponents, nDim = means.shape
    x, y, z = covars.shape
    if x != nComponents | y != nDim | z != nDim:
        raise RuntimeError('means shape: %sCovar shape: %s'%(str(means), str(covars.shape)))
    gmm = skGMM(n_components= nComponents, covariance_type='full')
    gmm.means_ = means
    gmm.covars_ = covars
    gmm.weights_ = weights
    return gmm

class GMM():
    """Extend sklearn.mixture.GMM to extract conditional gmm

    sklearn.mixture.GMM does not provide methods for making
    conditionals or marginal GMMs. This class extends
     sklearn.mixture.GMM through encapsulation.
    """
    def __init__(self,  *args, **kwargs):
        self.gmm = skGMM(*args, covariance_type='full', **kwargs)

    def findBestNumComponents(self, X, criteria = 'bic'):
        """Return best number of components.
           criteria must be one of 'bic' or 'aic'
        """
        nRows, nDim = X.shape
        icList = []
        for i in range(nDim//3, nDim):
            gmm = skGMM(n_components = i,  covariance_type= 'full')
            gmm.fit(X)
            getattr(gmm, criteria)
            icList.append((i, getattr(gmm, criteria)(X)))
        idx = np.argmin(np.array(icList)[:,1])
        if DEBUG:
            for i, (nComponents, ic) in enumerate(icList):
                print nComponents, ic, '<--' if i == idx else ''
        return icList[idx][0]

    def fit(self, X, criteria='bic'):
        """Fit data
        """
        nComponents = self.findBestNumComponents(X, criteria = criteria)
        self.gmm.n_components = nComponents
        self.gmm.fit(X)

    def getCovars(self):
        return self.gmm.covars_

    def getMeans(self):
        return self.gmm.means_

    def getWeights(self):
        return self.gmm.weights_

    def getCondDist(self, Y):
        """Return conditional GMM.

        Inputs:  An array of length N_dimensions.
                 Each element corresponds to a dimension.
                 Each non-NaN element will be locked.

        Outputs: A tuple of (means, covars, weights) of the conditional
                 distribution. Dimension order is preserved.

        Adapted from pypr.
        Appendix A in C. E. Rasmussen & C. K. I. Williams, Gaussian Processes
        for Machine Learning, the MIT Press, 2006
        """
        idxNotSet = np.nonzero(np.isnan(Y))[0]
        idxSet = np.nonzero(True - np.isnan(Y))[0]
        idxNew = np.concatenate((idxNotSet, idxSet))
        y = Y[idxSet]
        newMeans = []
        newCovars = []
        fk = []
        for i in range(len(self.gmm.means_)):
            # Make a new co-variance matrix with same ordering
            newCovar = copy.deepcopy(self.gmm.covars_[i])
            newCovar = newCovar[:, idxNew]
            newCovar = newCovar[idxNew, :]
            ua = self.gmm.means_[i][idxNotSet]
            ub = self.gmm.means_[i][idxSet]
            Sbb = newCovar[len(idxNotSet):, len(idxNotSet):]
            L = np.linalg.inv(newCovar)
            Laa = L[0:len(idxNotSet), 0:len(idxNotSet)]
            Lab = L[0:len(idxNotSet), len(idxNotSet):]
            cen = ua - np.dot(np.dot(np.linalg.inv(Laa), Lab), (y - ub))
            cov = np.linalg.inv(Laa)
            newMeans.append(cen)
            newCovars.append(cov)
            # Used for normalizing the mc
            fk.append(multivariate_normal.pdf(Y[idxSet], ub, Sbb))

        # Normalize the weights: p(X|Y) = p(Y,X) / p(Y) using the marginal dist.
        fk = np.array(fk).flatten()
        newWeights = (self.gmm.weights_*fk)
        newWeights = newWeights / np.sum(newWeights)
        return (np.array(newMeans), np.array(newCovars), np.array(newWeights))

