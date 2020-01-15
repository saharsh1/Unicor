# =======================================================================
# =======================================================================
# In[1]: Imports
# =======================================================================
# =======================================================================
import numpy as np
from astropy.io import fits
from engine.Spectrum import Spectrum
from engine.Template import Template
from engine.CCF1d import CCF1d
import matplotlib.pyplot as plt
from plotly import offline as py
import easygui


# =======================================================================
# =======================================================================
# In[2]: Load Data
# =======================================================================
# =======================================================================
# Select a fits file of NRES spectra to analyze (using easygui package)
NRESspec_fname = easygui.fileopenbox(msg=None, title='Select NRES data:',
                                     default='D:/Dropbox/Tess_TOI/',
                                     filetypes=None, multiple=False)

# Open the fits file using fits module from astropy.
hdul = fits.open(NRESspec_fname)
s = hdul[2].data
w = hdul[6].data*10*(1+-19.17/299792.458)   # NRES gives the data in nm. Convert to Ang.


# =======================================================================
# =======================================================================
# In[2]: Load Template
# =======================================================================
# =======================================================================
# Select a fits file of NRES spectra to analyze (using easygui package)
# NRESspec_fname = easygui.fileopenbox(msg=None, title='Select NRES data:',
#                                     default='D:/Dropbox/Tess_TOI/',
#                                     filetypes=None, multiple=False)
#
# Open the fits file using fits module from astropy.
hdult = hdul # fits.open(NRESspec_fname)
st = hdult[2].data
wt = hdult[6].data*10  # NRES gives the data in nm. Convert to Ang.


# =======================================================================
# =======================================================================
# In[3]: Aaarange data
# =======================================================================
# =======================================================================
sp = Spectrum(wv=w, sp=s)
sp.SpecPreProccess()
sp.TrimEdges(100)
# Set the spectrum to be the template
templatesp = Spectrum(wv=wt, sp=st).SpecPreProccess()

template = Template(template=templatesp)

ccf = CCF1d().CrossCorrelateSpec(spec=sp, template=template, dv=0.15, VelBound=25)
ccf.CombineCCFs()

# =======================================================================
# =======================================================================
# In[4]: Plots...
# =======================================================================
# =======================================================================
# Some nice plots of the spectrum, before and after the analysis
# %config InlineBackend.figure_format = "png"
# py.init_notebook_mode()
plt.figure(figsize=(13, 6))
plt.plot(ccf.Corr['vel'], ccf.Corr['combined'])
plt.grid()
plt.show()
# py.iplot_mpl(plt.gcf())


# =======================================================================
# =======================================================================
# In[3]: Plots...
# =======================================================================
# =======================================================================
# Some nice plots of the spectrum, before and after the analysis
# %config InlineBackend.figure_format = "png"
# py.init_notebook_mode()
# plt.figure(figsize=(10, 6))
# plt.grid()
# for I in [0]:
#    plt.plot(sp.wv[I], sp.sp[I])
#    plt.plot(template.model.wv[I], template.model.sp[I]-200)

# py.iplot_mpl(plt.gcf())

# #In[4]: Correct systematics
# Arrage the correlation matrix
# CorrMat = np.array([[val for val in CorrVals]
#                    for CorrVals in ccf.corr], ndmin=2, dtype=float)
# Read the number of orders in the spectrum
#Nord = CorrMat[:, 0].size
# CombinedCorr = (np.prod(1-np.multiply(CorrMat, CorrMat), axis=0))
