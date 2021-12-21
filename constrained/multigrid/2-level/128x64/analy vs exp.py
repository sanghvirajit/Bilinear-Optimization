import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

mu = 1;
x = 0.5;

T = np.linspace(0, 1, 129)

#PLOT Y
Y = []

for j in range(len(T)):
	t = T[j]
	y = t * math.sin(math.pi * x)
	Y.append(y)
		
y4 = pd.read_csv("y_experimental_{}.txt".format(x), header=None)

fig = plt.gcf()
plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

plt.plot(T, y4, label="y final")
plt.plot(T, Y, label="y analytical")

plt.xlabel("Time")
plt.ylabel("y")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical Y_{}.png".format(x))
plt.show()

#PLOT P
P = []

for j in range(len(T)):
	t = T[j]	
	p = mu * (t-1) * math.sin(math.pi * x)
	P.append(p)

p4 = pd.read_csv("p_experimental_{}.txt".format(x), header=None)

fig = plt.gcf()

plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

plt.plot(T, p4, label="p final")
plt.plot(T, P, label="p analytical")

plt.xlabel("Time")
plt.ylabel("p")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical P_{}.png".format(x))
plt.show()

#Plot control U
u_a = -0.1
u_b = 0.1
U = []

for i in range(0, len(P)):
  
  U.append(-1.0 * (1/mu) * P[i] * Y[i])
    
  if(U[i] > u_b):
  	U[i] = u_b;        
	
  if(U[i] < u_a):
  	U[i] = u_a;

u1 = pd.read_csv("u_const_experimental_{}.txt".format(x), header=None)
u2 = pd.read_csv("u_unconst_experimental_{}.txt".format(x), header=None)

fig = plt.gcf()

plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

plt.plot(T, U, label="u const analytical")
plt.plot(T, u1, label="u const exp")
plt.plot(T, u2, label="u unconst exp")

plt.xlabel("Time")
plt.ylabel("u")

plt.legend()
plt.grid()

fig.savefig("analytical U_const_{}.png".format(x))
plt.show()
