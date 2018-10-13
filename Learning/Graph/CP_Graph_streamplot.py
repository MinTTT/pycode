import numpy as np
import matplotlib.pyplot as plt


Y, X = np.mgrid[0:1:10j, 0:1:10j]
U = Y - X
V = (1/(1+X)) - Y
speed = np.sqrt(U*U + V*V)

Protein = np.arange(0, 1, 0.05)

Y_1 = Protein
Y_2 = (1/(1+Protein))

plt.figure(figsize=(9, 8))
#  Varying density along a streamline
streamline = plt.streamplot(X, Y, U, V, density=[3, 3],color=U, cmap='coolwarm')
plt.colorbar(streamline.lines)
plt.plot(Protein, Y_1)
plt.plot(Protein, Y_2)
plt.xlim(0,1)
plt.show()