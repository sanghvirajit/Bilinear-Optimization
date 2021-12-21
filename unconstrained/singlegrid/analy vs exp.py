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

#PLOT U
U = []

for j in range(len(T)):
	t = T[j]
	u = -t * (t-1) * math.sin(math.pi * x) * math.sin(math.pi * x)
	U.append(u)

u1 = pd.read_csv("u_experimental_{}.txt".format(x), header=None)

fig = plt.gcf()

plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

plt.plot(T, U, label="u unconst analytical")
plt.plot(T, u1, label="u unconst exp")

plt.xlabel("Time")
plt.ylabel("u")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical U_unconst_{}.png".format(x))
plt.show()
