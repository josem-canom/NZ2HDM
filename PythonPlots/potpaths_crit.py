import numpy as np
import pandas as pd
from math import atan

pdat = np.array(pd.read_csv('./Jose_crit.txt', sep='\t', header=None))
nfds = 4

h1path = -pdat[:,0]

h2path = -pdat[:,1]

a2path = -pdat[:,2]

spath = pdat[:,3]

thetapath = np.zeros(len(pdat))
for i in range(len(h1path)-1):
    thetapath[i] = -atan(a2path[i]/h2path[i])
thetapath[-1] = thetapath[-2]
#print(thetapath[-1])
