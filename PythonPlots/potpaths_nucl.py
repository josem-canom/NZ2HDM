import numpy as np
import pandas as pd
from math import atan

pdatbf = np.array(pd.read_csv('/home/jmc/Documents/1A_DM&BA/Jupyter_Notebooks/instanton_Phi.txt'))
nfds = 4

pdat = np.zeros(shape = (len(pdatbf), nfds))
for i in range(len(pdat)):
    pdat[i,:] = np.fromstring(pdatbf[i,0], dtype=float, sep=' ')

h1path = pdat[:,0]

h2path = pdat[:,1]

a2path = pdat[:,2]

spath = pdat[:,3]

thetapath = np.zeros(len(pdat))
for i in range(len(h2path)):
    thetapath[i] = -atan(a2path[i]/h2path[i])
thetapath[-1] = thetapath[-2]
#print(thetapath[-1])

pdatken = np.array(pd.read_csv('./Jose_nucl.txt', sep='\t', header=None))

#pdatken = np.zeros(shape = (len(pdatbfken), nfds))
#for i in range(len(pdatken)):
#    pdatken[i,:] = np.fromstring(pdatbfken[i,0], dtype=float, sep=' ')

h1pathken = -pdatken[:,0]

h2pathken = -pdatken[:,1]

a2pathken = -pdatken[:,2]

spathken = pdatken[:,3]

thetapathken = np.zeros(len(pdatken))
for i in range(len(h1pathken)-1):
    thetapathken[i] = -atan(a2pathken[i]/h2pathken[i])
thetapathken[-1] = thetapathken[-2]
#print(thetapathken[-1])
