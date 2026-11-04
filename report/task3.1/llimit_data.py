'''
Analyze infinite volume and continuity limit of the ground state
Part 1: Calculating the observable <x^2>

Data is expected to be structured in folders, all simulation for one value of
eps^2 grouped in a folder:

report
 |> this_file.py
 '> llimit
     |> llimit_{eps^2 value 1}
     |> llimit_{eps^2 value 2}
     |> llimit_{eps^2 value 3}
     '> llimit_{eps^2 value ...}
         |> groundstate_v2_N*   The files with simulated data
         '> groundstate_v2_N*   Keep the original filenames ans -types!

Data will be saved in the llimit directory as data_{eps^2 value}.txt for the
next script to read
'''

import os
import numpy as np
from glob import glob

os.chdir(os.path.dirname(__file__))

def getdtype(n, d):
    '''Create the dtype of a line in the output files'''
    return np.dtype([
        ("norm", float),
        ("e", float),
        ("ekin", float),
        ("epot", float),
        ("pos", float, (d, )),
        ("upos", float, (d, )),
        # ("mom", float, (d, )), # only exists in older versions
        # ("umom", float, (d, )), # only exists in older versions
        ("wf", complex, (n**d, ))
        ])

def converter(x):
    '''Complex number converter. The default converter from numpy has a bug'''
    return complex(x.decode())

def getconverters(n, d):
    '''Return the converters dictionary'''
    # return {4+4*d + i : converter for i in range(n**d)} # was required in older versions
    return {4+2*d + i : converter for i in range(n**d)}

def getparams(filename):
    '''Extract parameters N, D, mu and eps2 from the filename'''
    str, eps2 = filename.split("_epssq")
    str, mu = str.split("_mu")
    str, d = str.split("_D")
    str, n = str.split("_N")
    return int(n), int(d), float(mu), float(eps2)

def geteps2(folder):
    '''Extract parameter eps2 from the folder name'''
    str = os.path.basename(folder)
    eps2 = str.split("_")[-1]
    return float(eps2)


folders = glob("llimit/llimit_*")
for folder in folders:
    # exclude anything that is not a folder with data inside
    if os.path.isfile(folder) or "plot" in folder or "exclude" in folder:
        continue
    print("Processing folder", folder)

    # get eps^2 and all files simulated with this value
    eps2 = geteps2(folder)
    files = glob(f"llimit/llimit_{eps2:.6e}/groundstate_v2_N*")

    # prepare arrays for expectation value and simulation range of N
    x2arr = np.empty(len(files))
    narr = np.empty_like(x2arr)

    for i, filename in enumerate(files):
        n, d, _, _ = getparams(filename)

        # retrieve data
        data = np.genfromtxt(filename,
            dtype = getdtype(n, d),
            converters = getconverters(n, d)
        )

        # position operator
        pos = np.arange(-n/2, n/2, 1)
        # calculate and save expectation value
        exp_pos2 = np.sum( np.real(data["wf"] * np.conj(data["wf"])) * pos * pos )
        x2arr[i] = exp_pos2 * eps2
        narr[i] = n

    # sort and save results
    order = np.argsort(narr)
    np.savetxt(f"llimit/data_{eps2:.6e}.txt", np.transpose([narr[order], x2arr[order]]))
