import numpy as np
import os
from astropy.io import fits
from engine.Spectrum import Spectrum
from engine.Template import Template
from engine.CCF1d import CCF1d

# In[]:
files = []
# r=root, d=directories, f = files
for r, d, f in os.walk('D:\\Dropbox\\TESS_TOI\\TOI677\\lsc'):
    for file in f:
        if '.fits' in file:
            files.append(os.path.join(r, file))

# In[]:
obs = []
for fname in files:
    targetname = os.path.basename(fname)[:-4]
    hdul = fits.open(fname)
    s = hdul[2].data
    w = hdul[6].data   # NRES gives the data in nm. Convert to Ang.
    sp = Spectrum(wv=w, sp=s)
    sp.SpecPreProccess()
    sp.TrimEdges(100)
    obs.append(sp)

# In[]:
ccf = []
listNames = []
for I, templatesp in enumerate(obs):
    template = Template(template=templatesp)
    for J in np.arange(I+1, len(obs)-1):
        ccf_tmp = CCF1d().CrossCorrelateSpec(spec=sp, template=template, dv=0.1, VelBound=5)
        ccf_tmp.CombineCCFs()
        ccf.append(ccf_tmp.Corr['combined'])
        listNames.append(str(I) + ',' + str(J))
