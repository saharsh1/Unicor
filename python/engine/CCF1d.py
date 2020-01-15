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

import sys
from scipy import interpolate
import numpy as np
import pycorrelate as pyc


class CCF1d:
    # =============================================================================
    # =============================================================================
    def __init__(self, Corr=[]):
        '''
        Input: All input is optional, and needs to be called along
               with its keyword. Below appears a list of the possible input
               variables.

              spec: A list of spectrum objects
              template: template object
              name: str
        '''
        self.spec = Corr

    # =============================================================================
    # =============================================================================
    def CrossCorrelateSpec(self, spec, template, dv=None, VelBound=100):
        '''
        Input: All input is optional, and needs to be called along
               with its keyword. Below appears a list of the possible input
               variables.

                template - the template for the CCF (see template class)
                dv - scalar. If the spectrum is not logarithmically
                    evenly-spaced it is resampled, with a stepsize
                    determined according to dv (in km/s). Default 0.1
                VelBounds - scalar. The velocity bounds for the CCF.
                            detemined according to [-VelBound, VelBound]
                            Default is 100 km/s

        Output: self.corr - the CCF values
                self.velocity - the corresponding velocity shift.
        '''
        # Initialize:
        # ----------
        c = 299792.458  # The speed of light in km/sec

        # In case that the spectum is not logarithmically spaced,
        # it must be interpolated to a logarithmically evenly-spaced
        # grid. The parameter of the grid is dv [km/s] (default = 0.1).
        if dv is None:
            if not spec.GridInfo or spec.GridInfo['type'] == 'linear':
                dv = 0.1
                spec.InterpolateSpectrum(dv, InterpMethod='log')

            # If the data is already logarithmically spaced, then read the dv
            # of the wavelegth grid (used to set the velocity axis of the CCF)
            elif spec.GridInfo['type'] == 'log':
                dv = spec.GridInfo['delta']
        else:
            if not spec.GridInfo or spec.GridInfo['type'] == 'linear':
                spec.InterpolateSpectrum(dv, InterpMethod='log')
            elif spec.GridInfo['delta'] != dv:
                spec.InterpolateSpectrum(dv, InterpMethod='log')
                spec.GridInfo['delta'] = dv

        # The cross correlation is performed on a velociry range defined by
        # the user. The range is converted to a number of CCF lags, using
        # the velocity spacing dv. Default is 100 km/s.
        Nlags = VelBound//dv
        Nord = len(spec.wv)

        # Calculate the velocity from the lags
        V = dv * np.arange(-Nlags-1, Nlags)

        corr = np.full((len(spec.wv), len(V)), np.nan)
        SpecVel = np.full((len(spec.wv), 1), np.nan)

        for I, w in enumerate(spec.wv):
            # In order for the CCF to be normalized to the [-1,1] range
            # the signals must be divided by their standard deviation.
            s = spec.sp[I]
            sLen = len(s)
            NormFac = np.sqrt(pyc.ucorrelate(s, s, 1)/sLen)
            s = s / NormFac

            # Interpolate the template to the wavelength scale of the
            # observations. We assume here that the template is broadened
            # to match the width of the observed line profiles.
            interpX = np.asarray(template.model.wv[I][:]) * (1-VelBound/c)
            interpY = np.asarray(template.model.sp[I][:])
            InterpF = interpolate.interp1d(interpX,
                                           interpY,
                                           kind='quadratic')

            GridChoice = np.logical_and(w > np.min(interpX),
                                        w < np.max(interpX))
            wGrid = np.extract(GridChoice, w)
            spT = InterpF(wGrid)

            # Normazlize the template, in the same manner the observations are
            # normazlized:
            tLen = len(spT)
            NormFac = np.sqrt(pyc.ucorrelate(spT, spT, 1)/tLen)
            spT = spT/NormFac

            # The data and the template are cross-correlated.
            # We assume that the wavelengths are logarithmically-evenly spcaed
            C = pyc.ucorrelate(spT, s, 2*Nlags+1)/sLen

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

            corr[I, :] = C
            SpecVel[I] = vPeak

            # Show Progress
            n = int((30+1) * float(I) / Nord)
            sys.stdout.write("\rCorrelating: [{0}{1}]".format('#' * n, ' ' * (30 - n)))
        sys.stdout.write("\n")

        self.Corr = {'vel': V, 'corr': corr}
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
#        CorrMat = np.array([[val for val in CorrVals]
#                            for CorrVals in self.Corr['corr']],
#                           ndmin=2, dtype='float')
        CorrMat = self.Corr['corr']
        # Read the number of orders in the spectrum
        Nord = CorrMat.shape[0]

        # Combine the CCFs according to Zucker (2003, MNRAS), section 3.1
        CombinedCorr = 1-(np.prod(1-CorrMat**2, axis=0))**(1/Nord)

        # Return the corresponding velocity grid.
        self.Corr['combined'] = CombinedCorr

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
