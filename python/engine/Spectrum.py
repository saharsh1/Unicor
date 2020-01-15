#                   ------------------------------
#                     Spectrum.py (UNICOR class)
#                   ------------------------------
# This file defines the "Spectrum" class. An object of this class stores
# the measured spectrum (wavelength values and their corresponding flues)
# and some basic information about the measurement: time of observation,
# baricentric velocity correction, etc. See more details under __init__.
# The spectrum can be single- or multi-order.
#
# A Spectrum class stores the following methods:
# InterpolateSpectrum - resamples the spectrum on a linear or logarithmic scale
# TrimEdges - Cuts the edges of the spectrum (mostly to remove zero paddding)
# FilterSpectrum - A Butterworth bandpass filter
# ApplyCosineBell - Apllies a Tuckey window on the data
# RemoveCosmics - Removes outliers that deviate above the measured spectrum
# SpecPreProccess - Calls all the above a suggested order.
#
# Dependencies: numpy, scipy and astropy.

import numpy as np
from scipy import interpolate, signal
from astropy.stats import sigma_clip, mad_std


class Spectrum:

    # =============================================================================
    # =============================================================================
    def __init__(self, wv=[], sp=[],
                 bjd=0, bcv=0, GridInfo=[], name='Anonymous'):
        '''
        Input: All input is optional, and needs to be called along
               with its keyword. Below appears a list of the possible input
               variables.

              wv, sp: lists data from a single observation.
                      Each entry in the list represents a different order.
              name: str, the name of the target.
              bcv: float, barycentric correction to the velocity.
              bjd: barycentric julian day.
        '''

        self.wv = wv
        self.sp = sp
        self.bcv = bcv
        self.bjd = bjd
        self.GridInfo = GridInfo
        self.name = name

# =============================================================================
# =============================================================================
    def SpecPreProccess(self, TrimNum=150, InterpNum=0.5, RemCosmicNum=3,
                        FilterLC=3, FilterHC=0.15, alpha=0.3):
        '''
        This function applies the default NRES preprocessing to the spectra.
        In the case where not all parameters are suppliled,
        default values are assigned and used.
        '''

        print('SpecPreProccess... ', end='')

        self.TrimEdges(TrimNum)
        self.InterpolateSpectrum(InterpNum)
        self.RemoveCosmics(sigma_upper=RemCosmicNum)
        self.FilterSpectrum(lowcut=FilterLC, highcut=FilterHC, order=1)
        self.ApplyCosineBell(alpha=alpha)

        print('Done.')
        return self

    # =============================================================================
    # =============================================================================
    def InterpolateSpectrum(self, delta, **kwargs):
        '''
        This function resamples a list of wavelength and spectra
        so that the wavelength values will be evenly spaced.

        INPUT:  self.wv - list of 1D wavelength arrays.
                self.sp - list of 1D flux arrays, correspond the wv arrays.
                delta   - positive float.
                          If InterpMethod is linear, the step of the evenly
                          spaced array is taken to be delta*(min delta lambda).
                          In this case delta is a unitless parameter.
                          If InterpMethod is log, then the array is resampled
                          to be logarithmically evenly-spaed. The step is
                          taken to be log(1+delta/c), where c is the speed of
                          speed of light in km/sec, and so is delta.
                InterpMethod - string. 'linear' to generate an evenly spaced
                               wavelength vector. 'log' to generate a
                               logarithmically evenly-spaced vector.
        OUTPUT: self.wv - list of evenly spaced wavelegth arrays.
                self.sp - list of flux values, interpolated at wvI.
                self.GridInfo - type of spacing used ('log', 'linear'), the
                                step size delta and it units.

        Notes: The input, wv and sp, are data from a single observation.
               Each entry in the list represents a different order.
               The step size of the interpolated wavelength array is set
               to be f times the minimal (non zero, positive) wavelength
               step in the input data.
               Original data is saved as wvRaw & spRaw.

        Created: Sahar Shahaf, 19/4/2019
        Last Revision: ---
        '''
        # Initialize data
        if 'InterpMethod' in kwargs:
            InterpMethod = kwargs['InterpMethod']
            c = 299792.458  # The speed of light in km/sec
        else:
            InterpMethod = 'linear'

        wvI = []
        spI = []

        # Interpolate per order
        for I, w in enumerate(self.wv):

            # Define the evenly spaced grid
            if InterpMethod == 'linear':
                dlam = delta*np.min(np.diff(w))
                wvI.append(np.arange(min(w), max(w), dlam))
            elif InterpMethod == 'log':
                dlam = np.log(1+delta/c)
                medianDiff = np.median(np.diff(w))
                wvI.append(np.exp(
                           np.arange(np.log(min(w)+medianDiff),
                                     np.log(max(w)-medianDiff),
                                     dlam)))

            # Inpterpolate the flux
            InterpF = interpolate.interp1d(w, self.sp[I], kind='quadratic')
            spI.append(InterpF(wvI[I]))

        # Retrun variables.
        self.wv = wvI
        self.sp = spI
        if InterpMethod == 'log':
            self.GridInfo = {'type': 'log',
                             'delta': delta,
                             'delta_units': 'velocity'}
        else:
            self.GridInfo = {'type': 'linear',
                             'delta': delta,
                             'delta_units': 'wavelegth'}
        return self

    # =============================================================================
    # =============================================================================
    def TrimEdges(self, Ntrim):
        '''
        In some cases the spectra is zero-padded at the edges, which may affect
        the filtering stage. This routine trims the edges of the data.
        The number of trimmed indices is given as Ntrim.
        If Ntrim is not provided, all points that equal to zero,
        are removed.

        INPUT:  self.wv - list of 1D wavelength arrays.
                self.sp - list of 1D flux arrays, correspond the wv arrays.
                Ntrim   - Integer. Number of points to remove from each side
                          of the spectrum.
        OUTPUT: self.wv - list of wavelegth arrays, trimmed.
                self.sp - list of flux values, trimmed.

        Created: Sahar Shahaf, 3/5/2019
        Last Revision: ---
        '''
        # Initialize data
        wvT = []
        spT = []

        # Interpolate per order
        for I, w in enumerate(self.wv):

            Li = next((j for j, x in enumerate(self.sp[I]) if x), None)
            Ri = next((j for j, x in enumerate(np.flip(self.sp[I])) if x), None)

            if Ri > 0:
                s = self.sp[I][Li:-Ri]
                spT.append(s)
                w = w[Li:-Ri]
                wvT.append(w)
            else:
                s = self.sp[I][Li:]
                spT.append(s)
                w = w[Li:]
                wvT.append(w)

        self.wv = wvT
        self.sp = spT
        return self

    # =============================================================================
    # =============================================================================
    def FilterSpectrum(self, lowcut, highcut, order):
        '''
        This function applies a Butterworth fitler to the input
        spectrum, removes the low-pass instrumental repsponse and the high-pass
        instrumental noise.

        INPUT:  self.wv - list of 1D wavelength arrays (evenly sampled).
                self.sp - list of 1D flux arrays, correspond to the wv arrays.
                lowcut  - float. Stopband freq for low-pass filter. Given in
                          units of the minimal frequency (max(w)-min(w))**(-1)
                highcut - float. Stopband freq for the high-pass filter. Given
                          in units of the Nyquist frequency.
                order -   integer. The order of the Butterworth filter.

        OUTPUT: spF - list of filtered flux values at wv.

        Created: Sahar Shahaf, 19/4/2019
        Last Revision: ---
        '''
        # Initialize data
        spF = []

        # Filter the spectrum per order
        for I, w in enumerate(self.wv):
            # the bandpass values are normalized to the Nyquist frequency.
            # here we assume that the spectrum in evenly samples, therefore:
            nyq = 0.5/np.abs(w[1]-w[0])

            # The minimal frequency is set by the full range of the datapoints
            # here we use express the minimal frequency in terms of the Nyquist
            df = (max(w)-min(w))**(-1)/nyq

            # the lowpass filter is provided in terms of the minimal frequency.
            # it is therefore expressed by the Nyquist frequency.
            low = lowcut*df

            # The highpass frequency is given in terms of the Nyquist
            high = highcut

            # Define the Butterworth filter
            b, a = signal.butter(order, [low, high], btype='band')

            # Winsorize the data
            s = self.sp[I]-np.median(self.sp[I])

            # Filter the data
            spF.append(
                signal.filtfilt(b, a, s, padlen=10*max(len(b), len(a))))

        self.sp = spF
        return self

    # =============================================================================
    # =============================================================================
    def ApplyCosineBell(self, alpha=0.1):
        '''
        This function applies a cosine bell (Tukey Window) to the spectrum,

        INPUT:  self.wv - list of 1D wavelength arrays (evenly sampled).
                self.sp - list of 1D flux arrays, correspond to the wv arrays.
                alpha -   (optional) float. Shape parameter of the Tukey win.

        OUTPUT: spC - list of flux values at wv,
                      after the cosine bell was applied

        Created: Sahar Shahaf, 21/4/2019
        Last Revision: ---
        '''
        # Initialize data
        spC = []

        # Interpolate per order
        for I, w in enumerate(self.wv):
            window = signal.tukey(len(w), alpha)
            spC.append(
                    np.multiply(self.sp[I], window))

        self.sp = spC
        return self

# =============================================================================
# =============================================================================
    def RemoveCosmics(self, **kwargs):
        '''
        This function removes positive outliers, i.e., points that
        significantly deviate above the measured spectra. Significant deviation
        is defined as points that deviate more than (sigma_upper)X(1.48mad_std)
        from the main part of the spectrum. Points that deviate below are
        not removed.

        INPUT:  self.wv - list of 1D wavelength arrays (evenly sampled).
                self.sp - list of 1D flux arrays, correspond to the wv arrays.
                sigma_upper - float. Number of sigma beyond which points will
                              be removed. (Optional. default=3).

        OUTPUT: spR - list of flux values at wvR,
                      after exclusion of the outliers.

        Created: Sahar Shahaf, 4/5/2019
        Last Revision: ---
        '''

        # Initialize data
        if 'sigma_upper' in kwargs:
            sigma_upper = kwargs['sigma_upper']
        else:
            sigma_upper = 3

        spR = []
        wvR = []

        # Interpolate per order
        for I, s in enumerate(self.sp):
            filtered_data = sigma_clip(s, sigma_lower=np.Inf,
                                       sigma_upper=sigma_upper,
                                       masked=True, cenfunc='median',
                                       stdfunc=mad_std)
            spRtmp = [x for i, x in enumerate(s) if filtered_data[i]]
            spR.append(spRtmp)
            wvRtmp = [x for i, x in enumerate(self.wv[I]) if filtered_data[i]]
            wvR.append(wvRtmp)

        self.sp = spR
        self.wv = wvR
        return self
