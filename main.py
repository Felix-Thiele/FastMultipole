import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from FastMultipole import FastMultipole
from BruteForce import BruteForce
from BarnesHut import BarnesHut

bodies = [[1,2,1], [4.2,4.2,1], [1.2,1,1], [2.7,2.3,1], [0.3,2.6,1], [2.8,0.1,1], [0.8, 1.8, 1]]
FM = FastMultipole(np.array(bodies), 30)
BF = BruteForce(np.array(bodies))
BH = BarnesHut(np.array(bodies))

pot, forces = FM.calc_forces()
bf_forces = BF.calc_forces()
bh_forces = BH.calc_forces()

plt.xlim(-1,5)
plt.ylim(-1,5)

for b in bodies:
    plt.plot(b[0], b[1], 'kx')
    plt.quiver(b[0], b[1], np.real(forces[tuple(b)]), np.imag(forces[tuple(b)]), color='r', scale=20 , headlength=7, headwidth=7)
    plt.quiver(b[0], b[1], np.real(bf_forces[tuple(b)]), np.imag(bf_forces[tuple(b)]), color='b', scale=20, headlength=5, headwidth=5)
    plt.quiver(b[0], b[1], np.real(bh_forces[tuple(b)]), np.imag(bh_forces[tuple(b)]), color='y', scale=20 , headlength=3, headwidth=3)

r_patch = mpatches.Patch(color='r', label='FastMultipole')
b_patch = mpatches.Patch(color='b', label='BruteForce')
y_patch = mpatches.Patch(color='y', label='BarnesHut')
plt.legend(handles=[r_patch, b_patch, y_patch])
plt.show()