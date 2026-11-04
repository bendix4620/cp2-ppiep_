
import os
import numpy as np
from glob import glob

def getdtype(n, d):
    return np.dtype([
        ("norm", float),
        ("e", float),
        ("ekin", float),
        ("epot", float),
        ("pos", float, (d, )),
        ("upos", float, (d, )),
        #("mom", float, (d, )),
        #("umom", float, (d, )),
        ("wf", complex, (n**d, ))
        ])

def converter(x):
    return complex(x.decode())

def getconverters(n, d):
    return {4+2*d + i : converter for i in range(n**d)}

def getparams(filename):
    str, steps = filename.split("_steps")
    str, tau = str.split("_tau")
    str, method = str.split("_method")
    str, eps2 = str.split("_epssq")
    str, mu = str.split("_mu")
    str, d = str.split("_D")
    str, n = str.split("_N")
    
    
    return int(n), int(d), float(mu), float(eps2), float(tau), int(method), int(steps)

def geteps2(folder):
    str = os.path.basename(folder)
    eps2 = str.split("q")[-2:][1]
    return float(eps2)


folders = glob("data/*")
for folder in folders:
    if os.path.isfile(folder) or "plot" in folder or "exclude" in folder:
        continue
    print("Processing folder", folder)
    
    eps2 = geteps2(folder)
    
    files = glob(f"data/epssq{eps2}/evolution_v2*")
    
    epotarr = np.empty(len(files))
    ekinarr = np.empty_like(epotarr)
    narr = np.empty_like(epotarr)
  

    for i, filename in enumerate(files):
        
        #print("Processing", filename)

        n, d, _, _, _, _, _ = getparams(filename)

        #data = np.genfromtxt(filename,
        #    dtype = getdtype(n, d),
        #    converters = getconverters(n, d)
        #)


        epot = np.genfromtxt(filename, usecols= 3)
        #exp_pos2 = np.sum( np.real(data["wf"] * np.conj(data["wf"])) * pos * pos )
        #x2arr[i] = exp_pos2 * eps2
        epotarr[i] = epot[-1]
        narr[i] = n 
        
        ekin = np.genfromtxt(filename, usecols=2)
        ekinarr[i] = ekin[-1]
    
    np.savetxt(f"x_vs_N_epssq{eps2}.txt", np.transpose([narr, epotarr, ]))
    np.savetxt(f"Ekin_vs_N_epssq{eps2}.txt", np.transpose([narr, ekinarr, ]))

print('successfull')


# %%

# %%
