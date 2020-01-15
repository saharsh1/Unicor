#                   ----------------------------
#                    Template.py (UNICOR class)
#                   ----------------------------
# This file defines the "Template" class. An object of this class stores
# a Spectrum object, saved in the self.model filed. This class has some
# additional attributes added to the ones of the Spectrum class.
#
# A Template class stores the following methods (under construction):
# GaussianBroadening - broadens the template with a Gaussian window,
#                      to account for the instrumental broadening.
# RotationalBroadening - broadens the template with a rotaional profile
#                        to account for stellar rotation.
# RetrievePhoenixTemplate - retrieves a model from the PHONIEX FTP.
#
# Dependencies: Spectrum class (of UNICOR).

from engine.Spectrum import Spectrum


class Template:

    # =============================================================================
    # =============================================================================
    def __init__(self, **kwargs):
        '''
        Input: All input is optional, and needs to be called along
               with its keyword. Below appears a list of the possible input
               variables.

              wv, sp: lists data from a single observation.
                      Each entry in the list represents a different order.
              name: str
        '''
        if 'template' in kwargs:
            self.model = kwargs['template']
            if 'vel' in kwargs:
                self.vel = kwargs['vel']
            else:
                self.vel = []
        else:
            self.model = Spectrum()
            self.vel = []

    # =============================================================================
    # =============================================================================
    def GaussianBroadening(self, **kwargs):
        pass

    # =============================================================================
    # =============================================================================
    def RotationalBroadening(self, **kwargs):
        pass

    # =============================================================================
    # =============================================================================
    def RetrievePhoenixTemplate(self, **kwargs):
        pass
