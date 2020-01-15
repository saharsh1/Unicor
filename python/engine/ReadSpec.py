#                   ------------------------------
#                     ReadSpec.py (UNICOR class)
#                   ------------------------------
# This file reads data from various types of instruments and stores it
# into a Spectrum class. Currently, the following spectrum of the following
# instruments are supported: NRES, APOGEE, ...
#
# An alternative to the default data types is to provide a .dat file.
#
# A Spectrum class stores the following methods:

import os
import numpy as np
from astropy.io import fits
from engine.Spectrum import Spectrum
import matplotlib.pyplot as plt
import easygui


def ReadSpec(DataType='dat', FileName=None):
