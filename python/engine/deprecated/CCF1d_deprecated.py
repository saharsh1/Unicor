#                   ----------------------------
#                      CCF1d.py (UNICOR class)
#                   ----------------------------
# This file defines the "CCF1d" class. An object of this class stores
# a Spectrum object, saved in the self.spec field and a Template object
# stored in the self.template field.
#
# A Template class has the following methods (under construction):
# CrossCorrelateSpec - multi-order cross correlation.
# CombineCCFs - sums the CCFs of a multi-order spectrum
#
# Dependencies: Spectrum and Template class (of UNICOR).
#               scipy, numpy and pycorrelate.

from scipy import interpolate
import numpy as np
import pycorrelate as pyc
#from astropy.modeling import models, fitting


class CCF1d:

    # =============================================================================
    # =============================================================================
    def __init__(self, **kwargs):
        '''
        Input: All input is optional, and needs to be called along
               with its keyword. Below appears a list of the possible input
               variables.

              spec: spectrum object
              template: template object
              name: str
        '''
        if 'template' and 'spec' in kwargs:
            self.template = kwargs['template']
            self.spec = kwargs['spec']
        else:
            self.template = []
            self.spec = []

    # =============================================================================
    # =============================================================================
    def CrossCorrelateSpec(self, **kwargs):
        '''
        Input: All input is optional, and needs to be called along
               with its keyword. Below appears a list of the possible input
               variables.

               dv - scalar. If the spectrum is not logarithmically
                    evenly-spaced it is resampled, with a stepsize
                    determined according to dv (in km/s). Default 0.1
                VelBounds - scalar. The velocity bounds for the CCF.
                            detemined according to [-VelBound, VelBound]

        Output: self.corr - the CCF values
                self.velocity - the corresponding velocity shift.
        '''
        # Initialize:
        # ----------
        # In case that the spectum is not logarithmically spaced,
        # it must be interpolated to a logarithmically evenly-spaced
        # grid. The parameter of the grid is dv [km/s] (default = 0.1).
        if not self.spec.GridInfo or self.spec.GridInfo['type'] == 'linear':
            if 'dv' in kwargs:
                dv = kwargs['dv']
                self.spec.InterpolateSpectrum(dv, InterpMethod='log')
            else:
                print('No dv input found. Default = 0.1 km/s.')
                dv = 0.1
                self.spec.InterpolateSpectrum(dv, InterpMethod='log')

        # If the data is already logarithmically spaced, then read the dv
        # of the wavelegth grid (used to set the velocity axis of the CCF)
        elif self.spec.GridInfo['type'] == 'log':
            dv = self.spec.GridInfo['delta']

        # The cross correlation is performed on a velociry range defined by
        # the user. The range is converted to a number of CCF lags, using
        # the velocity spacing dv. Default is 100 km/s.
        if 'VelBound' in kwargs:
            VelBound = np.abs(kwargs['VelBound'])
        else:
            print('No velocity bounds input found. Default = 100 km/s.')
            VelBound = 100

        Nlags = VelBound//dv + 1

        corr = []
        vels = []
        SpecVel = []

        # Cross-correlate per order
        # -------------------------
        for I, w in enumerate(self.spec.wv):
            # In order for the CCF to be normalized to the [-1,1] range
            # the signals must be divided by their standard deviation.
            s = self.spec.sp[I]
            sLen = len(s)
            NormFac = np.sqrt(pyc.ucorrelate(s, s, 1)/sLen)
            s = s / NormFac

            # Interpolate the template to the wavelength scale of the
            # observations. We assume here that the template is broadened
            # to match the width of the observed line profiles.
            InterpF = interpolate.interp1d(self.template.model.wv[I],
                                           self.template.model.sp[I],
                                           kind='quadratic')

            GridChoice = np.logical_and(w > np.min(self.template.model.wv[I]),
                                        w < np.max(self.template.model.wv[I]))
            wGrid = np.extract(GridChoice, w)
            spT = InterpF(wGrid)

            # Normazlize the template, in the same manner the observations are
            # normazlized:
            tLen = len(spT)
            NormFac = np.sqrt(pyc.ucorrelate(spT, spT, 1)/tLen)
            spT = spT/NormFac

            # The data and the template are cross-correlated.
            # We assume here that the wavelength scale is
            # logarithmically-evenly spcaed.
            CpositiveLag = pyc.ucorrelate(s, spT, Nlags)/sLen
            CnegativeLag = pyc.ucorrelate(np.flip(s), np.flip(spT), Nlags)/sLen

            # Combine the positive & negative sides of the CCF
            CnegativeLag = np.delete(CnegativeLag, 0)
            CnegativeLag = np.flip(CnegativeLag)
            C = np.concatenate((CnegativeLag, CpositiveLag), axis=0)

            # Calculate the velocity from the lags
            V = dv*(np.arange(0, len(C), 1) - len(C)//2)

            # Find the radial velocity by fitting a parabola to the CCF peaks
            try:
                MaxIndx = np.where(C == max(C))
                # Fit the data using a parabula to find the peak.
                # Currently the number of points to be fitted around
                # the peak is hard-coded to be 4. To be changes in next versions.
                vfit = V[MaxIndx[0][0]-4:MaxIndx[0][0]+4]
                cfit = C[MaxIndx[0][0]-4:MaxIndx[0][0]+4]
                pfit = np.polyfit(vfit, cfit, 2)
                vPeak = -pfit[1]/(2*pfit[0])

            except:
                vPeak = np.NaN

            corr.append(C)
            vels.append(V)
            SpecVel.append(vPeak)

        self.Corr = {'vel': vels, 'corr': corr}
        self.vels = SpecVel
        return self

# =============================================================================
# =============================================================================
    def CombineCCFs(self):
        '''
        Arrange the CCF as a  matrix, where the first index is the order
        number and the second index is the ccf values along the correlation
        lags, i.e., [#order,#lag].

        INPUT: none.
        OUTPUT:
        NOTE: the number of lags is assumed to be identical for all orders
        '''
        # Arrage the correlation matrix
        CorrMat = np.array([[val for val in CorrVals]
                            for CorrVals in self.Corr['corr']],
                           ndmin=2, dtype='float')

        # Read the number of orders in the spectrum
        Nord = CorrMat[:, 0].size

        # Combine the CCFs according to Zucker (2003, MNRAS), section 3.1
        CombinedCorr = 1-(np.prod(1-CorrMat**2, axis=0))**(1/Nord)

        # Return the corresponding velocity grid.
        self.CombinedCorr = {'vel': self.Corr['vel'][0],
                             'corr': CombinedCorr}

        return self

# =============================================================================
# =============================================================================
    def fwhm(x, y):
        """
        Determine full-with-half-maximum of a peaked set of points, x and y.
        Assumes that there is only one peak present in the datasset.
        The function uses a spline interpolation of order k.
        """

        half_max = max(y)/2.0
        s = interpolate.splrep(x, y - half_max)
        roots = interpolate.sproot(s)

        if len(roots) > 2:
            # The dataset appears to have multiple peaks
            return np.nan
        elif len(roots) < 2:
            # No proper peaks were found in the data set
            return np.nan
        else:
            return abs(roots[1] - roots[0])
