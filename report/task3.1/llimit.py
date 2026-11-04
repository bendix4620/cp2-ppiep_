'''
Analyze infinite volume and continuity limit of the ground state
Part 2: Fitting the limit <x^2>(N -> infinity)

Data is expected in the llimit directory as data_{eps^2 value}.txt
'''

import os, warnings
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
import scipy.optimize._optimize

os.chdir(os.path.dirname(__file__))

# fit function
def func(x, lim, a, b, c=1):
    '''f(x, lim, a, b, c) = lim + exp( -a(x - b)^c )'''
    return lim + np.exp(-a*(x-b)**c)

# The initial values for each value of eps^2, must be in order of the files
p0 = [
    [9.48609420e-01, 7.23897890e-01, 2.35343508e+01], # 1.0e-02 (eps^2 values)
    [9.47187310e-01, 3.32122110e-01, 7.98516607e+01], # 1.0e-03
    [9.46988512e-01, 1.07314737e-01, 2.53475060e+02], # 1.0e-04
    [9.46948565e-01, 3.38225305e-02, 8.01053614e+02], # 1.0e-05
    [9.46945731e-01, 1.06898207e-02, 2.53294411e+03], # 1.0e-06
    [9.46946222e-01, 3.37767605e-03, 8.00864420e+03], # 1.0e-07
    [9.46945849e-01, 1.06641490e-03, 2.53182677e+04], # 1.0e-08
    [9.47315490e-01, 4.09651650e-01, 5.90543949e+01], # 1.7e-03
    [9.46994144e-01, 1.42840832e-01, 1.90013366e+02], # 1.7e-04
    [9.46949849e-01, 4.51289619e-02, 6.00764919e+02], # 1.7e-05
    [9.46948156e-01, 1.42728143e-02, 1.89983482e+03], # 1.7e-06
    [9.46946585e-01, 4.50762680e-03, 6.00644019e+03], # 1.7e-07
    [9.46941850e-01, 5.26897710e-03, 2.05194201e+04], # 1.7e-08
    [9.47421220e-01, 4.82112360e-01, 4.33052458e+01], # 3.1e-03
    [9.47056430e-01, 1.98029170e-01, 1.43475832e+02], # 3.1e-04
    [9.46945745e-01, 5.90592260e-02, 4.49024300e+02], # 3.1e-05
    [9.46948044e-01, 1.89782189e-02, 1.42394070e+03], # 3.1e-06
    [9.46947752e-01, 6.01712950e-03, 4.50500268e+03], # 3.1e-07
    [9.46941850e-01, 5.26897710e-03, 1.55194201e+04], # 3.1e-08
    [9.48120780e-01, 6.49222180e-01, 3.25034153e+01], # 5.1e-03
    [9.47071840e-01, 2.46336050e-01, 1.06272885e+02], # 5.1e-04
    [9.46955702e-01, 7.94953894e-02, 3.37271892e+02], # 5.1e-05
    [9.46942527e-01, 2.51393173e-02, 1.06654336e+03], # 5.1e-06
    [9.46949510e-01, 8.03402491e-03, 3.37900155e+03], # 5.1e-07
    [9.46941850e-01, 5.26897710e-03, 1.18194201e+04]  # 5.1e-08
]
# p0 = [[9.47, 1, 1, 1]]*25

files = glob("llimit/data_*.txt")
for i, filename in enumerate(files):
    # exclude any folders or plots
    if "plot" in filename or "exclude" in filename:
        continue
    print("Processing", filename)
    name = os.path.basename(filename).split(".txt")[0]

    # Read in data
    n, x2 = np.genfromtxt(filename).T

    # exclude data far from the limit for the fits
    m = np.argmax(x2) + int(len(x2) / 10)
    if m >= len(x2)-1:
        print("WARNING: No asymptotic region found")

    # the uncertainty approximated by altering precision parameters of the simulation
    unc = np.zeros(len(x2)) + 1e-3

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        # fit
        popt, pcov = fit(func, n[m:], x2[m:], sigma=unc[m:], absolute_sigma=True, maxfev=100000, p0=p0[i])

    # save results numerically
    with open(f"llimit/plots/{name}.txt", "w") as export:
        export.write(f"popt\n{popt}\nstd\n{np.sqrt(np.diag(pcov))}\npcov\n{pcov}\nfunction\n{func.__doc__}")

    # save results graphically
    x = np.linspace(n[m], np.max(n), 1000)
    fig, ax = plt.subplots()
    line = ax.errorbar(n, x2, yerr=unc, fmt="x", ecolor="k", capsize=2, elinewidth=0.5, label="data")
    ax.plot(x, func(x, *popt), label="fit", color=line[0].get_color())
    plt.xlabel("$N$")
    plt.ylabel("$<(x/r)^2>$")
    plt.legend()
    plt.savefig(f"llimit/plots/{name}.pdf")
    plt.close()
