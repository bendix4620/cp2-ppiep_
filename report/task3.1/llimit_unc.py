'''
Analyze infinite volume and continuity limit of the ground state
Part 2.5, optional: Estimating a systematical error

This is done by fitting the limit <x^2>(N -> infinity) for only even and only
odd values of N and calculating the difference of these two limits.

The data file to calculate the difference with:'''
filename = "llimit/data_1.000000e-03.txt"

import os, warnings
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
import scipy.optimize._optimize

os.chdir(os.path.dirname(__file__))
print(os.getcwd())

def geteps2(filename):
    eps2 = filename.split("_")[-1].split(".txt")[0]
    return float(eps2)

def func(x, lim, a, b, c=1):
    '''f(x, lim, a, b, c) = lim + exp( -a(x - b)^c )'''
    return lim + np.exp(-a*(x-b)**c)

# initial fit parameters
p0 = [9.47187310e-01, 3.32122110e-01, 7.98516607e+01] # 1.0e-03

name = os.path.basename(filename).split(".txt")[0]

# read in data
n, x2 = np.genfromtxt(filename).T

# exclude data far from limit
m = np.argmax(x2) + int(len(x2) / 10) + 3
if m >= len(x2)-1:
    print("WARNING: No asymptotic region found")

unc = np.zeros(len(x2)) + 1e-3 # approximated uncertainty

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    # even and odd fit
    popte, pcove = fit(func, n[m  ::2], x2[m  ::2], sigma=unc[m  ::2], absolute_sigma=True, maxfev=100000, p0=p0)
    popto, pcovo = fit(func, n[m+1::2], x2[m+1::2], sigma=unc[m+1::2], absolute_sigma=True, maxfev=100000, p0=p0)

# save results numerically
with open(f"llimit/plots/double{name}.txt", "w") as export:
    export.write(f"even:\npopt\n{popte}\nstd\n{np.sqrt(np.diag(pcove))}\npcov\n{pcove}\nfunction\n{func.__doc__}\n")
    export.write(f"odd:\npopt\n{popto}\nstd\n{np.sqrt(np.diag(pcovo))}\npcov\n{pcovo}\nfunction\n{func.__doc__}\n")
    export.write(f"difference:\n{np.abs(popte[0]-popto[0])}")

# save results graphically
x = np.linspace(n[m], np.max(n), 1000)
fig, ax = plt.subplots()
line = ax.errorbar(n, x2, yerr=unc, fmt="x", ecolor="k", capsize=2, elinewidth=0.5, label="data")
ax.plot(x, func(x, *popte), label="fit even")
ax.plot(x, func(x, *popto), label="fit odd")
plt.ylim((0.945, 0.95))
plt.xlim((99.5, 121))
plt.xticks(np.arange(100, 121, 2))
plt.xlabel("$N$")
plt.ylabel("$<(x/r)^2>$")
plt.legend()
plt.savefig(f"llimit/plots/double{name}.pdf")
plt.close()
