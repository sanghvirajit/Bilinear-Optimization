import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

x = 0.5

w = pd.read_csv("w.txt", header=None);
df = pd.DataFrame(index=np.arange((w.size)/62), columns=np.arange(1))

n = 15

for i in range(0, 65):
	df.iloc[i, 0] = w.iloc[i*62 + n, 0]
	
df.to_csv("y_experimental_{}.txt".format(x), header=None, index=False)

n = 46

for i in range(0, 65):
	df.iloc[i, 0] = w.iloc[i*62 + n, 0]

df.to_csv("p_experimental_{}.txt".format(x), header=None, index=False)

u = pd.read_csv("u.txt", header=None);
df = pd.DataFrame(index=np.arange((u.size)/31), columns=np.arange(1))

n = 15

for i in range(0, 65):
	df.iloc[i, 0] = u.iloc[i*31 + n, 0]
	
df.to_csv("u_experimental_{}.txt".format(x), header=None, index=False)

u = pd.read_csv("u_projected.txt", header=None);
df = pd.DataFrame(index=np.arange((u.size)/31), columns=np.arange(1))

n = 15

for i in range(0, 65):
	df.iloc[i, 0] = u.iloc[i*31 + n, 0]
	
df.to_csv("u_projected_experimental_{}.txt".format(x), header=None, index=False)
