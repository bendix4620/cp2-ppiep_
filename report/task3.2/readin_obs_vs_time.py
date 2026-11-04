
import os
import numpy as np
from glob import glob



def gettau(file):
    str = os.path.basename(file)
    tau = str.split("_tau")[-2:][1]
    tau = tau.split("_")[-2:][0]
    return float(tau)

def geteps2(file):
    str = os.path.basename(file)
    eps2 = str.split("q")[-2:][1]
    eps2 = eps2.split("_method")[-2:][0]
    return float(eps2)


def getN(file):
    str = os.path.basename(file)
    N = str.split("_N")[-2:][1]
    N = N.split("_D1")[-2:][0]
    return int(N)

    
N = 1000
eps2 = 1.6e-05
tau = 0.01    
files = glob(f"data/tau_var/evolution_v2*")


for i, filename in enumerate(files):
    epot_t = np.genfromtxt(filename, usecols=3)
    ekin_t = np.genfromtxt(filename, usecols=2)
    tau = gettau(filename)
    time = np.empty(len(epot_t))
    for j in range(len(epot_t)):
        time[j] = j*tau
    
    
    np.savetxt(f"tau{tau}.txt", np.transpose([time, epot_t, ekin_t, ]))
    
files = glob(f"data/N_var/evolution_v2*")
tau = 0.01

for i, filename in enumerate(files):
    epot_t = np.genfromtxt(filename, usecols=3)
    ekin_t = np.genfromtxt(filename, usecols=2)
    N = getN(filename)
    time = np.empty(len(epot_t))
    for j in range(len(epot_t)):
        time[j] = j*tau
    
    
    np.savetxt(f"N{N}.txt", np.transpose([time, epot_t, ekin_t, ]))
    
files = glob(f"data/eps_var/evolution_v2*")
tau = 0.01

for i, filename in enumerate(files):
    epot_t = np.genfromtxt(filename, usecols=3)
    ekin_t = np.genfromtxt(filename, usecols=2)
    eps2 = geteps2(filename)
    time = np.empty(len(epot_t))
    for j in range(len(epot_t)):
        time[j] = j*tau
    
    
    np.savetxt(f"eps2{eps2}.txt", np.transpose([time, epot_t, ekin_t, ]))
    
    
   
   

    
    
