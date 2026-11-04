'''
Analyze infinite volume and continuity limit of the ground state
Part 3: Fitting the limit <x^2>(N -> infinity, eps^2 -> 0)

Data is expected in the llimit/plots directory as data_{eps^2 value}.txt
'''

import os
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit

os.chdir(os.path.join(os.path.dirname(__file__), "llimit/plots"))

def geteps2(filename):
    eps2 = filename.split("_")[-1].split(".txt")[0]
    return float(eps2)

# get files
files = glob("data_*.txt")
# prepare data arrays
eps2 = np.empty(len(files))
lim = np.empty(len(files))
ulim = np.empty(len(files))

for i, filename in enumerate(files):
    # read limit and uncertainty from file
    with open(filename, "r") as inline:
        # get lim
        line = inline.readline()
        line = inline.readline()
        l = line.replace("[", "").split(" ")
        for e in l:
            if "." in e:
                s = e
                break

        # get unc(lim)
        line = inline.readline()
        line = inline.readline()
        ul = line.replace("[", "").split(" ")
        for e in ul:
            if "." in e:
                us = e
                break

    # save values
    eps2[i] = geteps2(filename)
    lim[i] = float(s)
    ulim[i] = np.sqrt(float(us)**2 + 2.4e-4**2) # add uncertainty from even-odd-difference


def f(x, lim, a):
    '''f(x) = lim + ax^2'''
    return lim + a*x*x

# sort data
order = np.argsort(eps2)
eps2 = eps2[order]
lim = lim[order]
ulim = ulim[order]

# exclude data far from limit
slic = slice(0, -7, 1)
# fit
popt, pcov = fit(f, eps2[slic], lim[slic], sigma=ulim[slic], absolute_sigma=True, maxfev=100000, p0=[9.46930422e-01, 1])
# save results numerically
with open("eps2.txt", "w") as outline:
    outline.write(f"popt\n{popt}\nstd\n{np.sqrt(np.diag(pcov))}\npcov\n{pcov}\nfunction\n{f.__doc__}\n")
    outline.write(f"limit estimate\n{lim[0]}+-{ulim[0]}")

# save results graphically
fig, ax = plt.subplots()
line = ax.errorbar(eps2, lim, yerr=ulim, fmt="x", elinewidth=1, capsize=2, ecolor="k", label="data")
x = np.linspace(np.min(eps2), np.max(eps2[:]), 1000)
ax.plot(x, f(x, *popt), color=line[0].get_color(), label="fit")
plt.xlabel("$\\epsilon^2$")
plt.ylabel("$<(x/r)^2>(L=\\infty)$")
plt.legend()
plt.ylim((0.9465, 0.9493))
plt.xscale("log")
plt.savefig("eps2.pdf", format="pdf", bbox_inches="tight")
plt.show()
