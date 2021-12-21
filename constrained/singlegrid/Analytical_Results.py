import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import json

mu = 1;

Y = []
Z = []
P = []
U = []
F = []

T = np.linspace(0, 1, 65);
X = np.linspace(0, 1, 33)

print(T)
print(X)

## Save Y

for j in range(len(T)):
	t = T[j]
	for i in range(1, len(X)-1):
		x = X[i]
		y = t * math.sin(math.pi * x)
		Y.append(y)
	
with open('y.txt', 'w') as f:
    for item in Y:
        f.write("%s\n" % item)

## Save P

for j in range(len(T)):
	t = T[j]
	for i in range(1, len(X)-1):
		x = X[i]
		p = mu * (t-1) * math.sin(math.pi * x)
		P.append(p)

with open('p.txt', 'w') as f:
    for item in P:
        f.write("%s\n" % item)

## Save U
        
for j in range(len(T)):
	t = T[j]
	for i in range(1, len(X)-1):
		x = X[i]
		u = -t * (t-1) * math.sin(math.pi * x) * math.sin(math.pi * x)
		U.append(u)

with open('u_unconst.txt', 'w') as f:
    for item in U:
        f.write("%s\n" % item)
        
## Save Z

for j in range(len(T)):
	t = T[j]
	for i in range(1, len(X)-1):
		x = X[i]
		z = mu * math.sin(math.pi * x) - mu * math.pi * math.pi * (t-1) * math.sin(math.pi * x) - mu * t * (t-1)**2 * math.sin(math.pi * x) * math.sin(math.pi * x) * math.sin(math.pi * x) + t * 	   	               math.sin(math.pi * x)
		Z.append(z)

with open('z_unconst.txt', 'w') as f:
    for item in Z:
        f.write("%s\n" % item)

## Save F
        
for j in range(len(T)):
	t = T[j]
	for i in range(1, len(X)-1):
		x = X[i]
		f = - math.sin(math.pi * x) - math.pi * math.pi * t * math.sin(math.pi * x) - t**2 * (t-1) * math.sin(math.pi * x) * math.sin(math.pi * x) * math.sin(math.pi * x)
		F.append(f)

with open('f_unconst.txt', 'w') as f:
    for item in F:
        f.write("%s\n" % item)
        
##############################            
u_a = -0.1
u_b = 0.1
U = []

for i in range(0, len(P)):
	U.append(-1.0 * (1/mu) * P[i] * Y[i])
 		
	if(U[i] > u_b):
		U[i] = u_b;        
	if(U[i] < u_a):
		U[i] = u_a;
 
with open('u_const.txt', 'w') as f:
    for item in U:
        f.write("%s\n" % item)
        
## Save Z
Z1 = []

for j in range(len(T)):
	t = T[j]
	for i in range(1, len(X)-1):
		x = X[i]
		z1 = mu * math.sin(math.pi * x) - mu * math.pi**2 * (t-1) * math.sin(math.pi * x)
		Z1.append(z1)
		  	
for i in range(0, len(P)):

	z = Z1[i] + U[i] * P[i] + Y[i]
	Z.append(z)

with open('z_const.txt', 'w') as f:
    for item in Z:
        f.write("%s\n" % item)

## Save F
F1 = []
     
for j in range(len(T)):
	t = T[j]
	for i in range(1, len(X)-1):
		x = X[i]
		f1 = - math.sin(math.pi * x) - math.pi**2 * t * math.sin(math.pi * x)
		F1.append(f1)
		  		
for i in range(0, len(P)):

	f = F1[i] + U[i] * Y[i]
	F.append(f)  				

with open('f_const.txt', 'w') as f:
    for item in F:
        f.write("%s\n" % item)
