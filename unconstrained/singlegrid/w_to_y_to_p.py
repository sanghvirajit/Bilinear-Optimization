import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

x = 0.5
T = 129

w = pd.read_csv("w.txt", header=None);
df = pd.DataFrame(index=np.arange((w.size)/126), columns=np.arange(1))

n = 32

for i in range(0, T):
	df.iloc[i, 0] = w.iloc[i*126 + n, 0]
	
df.to_csv("y_experimental_{}.txt".format(x), header=None, index=False)

n = 96

for i in range(0, T):
	df.iloc[i, 0] = w.iloc[i*126 + n, 0]

df.to_csv("p_experimental_{}.txt".format(x), header=None, index=False)

u = pd.read_csv("u.txt", header=None);
df = pd.DataFrame(index=np.arange((u.size)/63), columns=np.arange(1))

n = 32

for i in range(0, T):
	df.iloc[i, 0] = u.iloc[i*63 + n, 0]
	
df.to_csv("u_experimental_{}.txt".format(x), header=None, index=False)
