#%%
from ast import Str
import os
import numpy as np
from glob import glob
import matplotlib.pyplot as plt


def gettau(filename):
    str = os.path.basename(filename)
    str = str.split("tau")
    tau = str[1].split(".txt")[-2:][0]
    return float(tau)

def geteps2(filename):
    str = os.path.basename(filename)
    str = str.split("eps2")
    eps2 = str[1].split(".txt")[-2:][0]
    return float(eps2)

def getN(filename):
    str = os.path.basename(filename)
    str = str.split("N")
    N = str[1].split(".txt")[-2:][0]
    return int(N)

D = 1
#------------------------------------------------------------------
#variing tau:
#------------------------------------------------------------------

files = glob(f"data/tau*")
eps2 = 1.6e-05
N = 1000
tau = []


ekin = []
x = []
time = []

for i,file in enumerate(files): 
    tau.append(gettau(file))
    ekin.append(np.genfromtxt(file, dtype=float, usecols=2))
    x.append(np.genfromtxt(file, dtype=float, usecols=1))
    time.append(np.genfromtxt(file, dtype=float, usecols=0))

diff_t001t01 = np.abs(x[0][-1] - x[3][-1])
diff_t001t05 = np.abs(x[0][-1] - x[2][-1])
diff_t001t1 = np.abs(x[0][-1] - x[1][-1])

print('V(x)\nVariing tau differences are: \n0.001-0.01 = {} \n0.001-0.05 = {} \n0.001-0.1 = {}'.format(diff_t001t01, diff_t001t05, diff_t001t1))

diff_t001t01 = np.abs(ekin[0][-1] - ekin[3][-1])
diff_t001t05 = np.abs(ekin[0][-1] - ekin[2][-1])
diff_t001t1 = np.abs(ekin[0][-1] - ekin[1][-1])


print('Ekin(x)\nVariing tau differences are: \n0.001-0.01 = {} \n0.001-0.05 = {} \n0.001-0.1 = {}'.format(diff_t001t01, diff_t001t05, diff_t001t1))
#-------------------------------------------
#x

fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(16,6))
fig.text(.5, .93, 'τ → 0, N = {0}, $Ɛ^2$ = {1}'.format(N, eps2), fontsize=18, ha = 'center')
fig.tight_layout(pad=6.0)

ax1.set_title('(a)', fontsize=16, ha = 'right')
ax2.set_title('(b)', fontsize=16, ha = 'right')
ax1.tick_params(axis= 'both', labelsize= 14)
ax2.tick_params(axis= 'both', labelsize= 14)

ax1.plot( time[0], x[0], color = 'darkblue', label = 'τ = {0}'.format(tau[0]))
ax1.plot( time[3], x[3], color = 'orange', label = 'τ = {0}'.format(tau[3]))
ax1.plot( time[2], x[2], color = 'firebrick', label = 'τ = {0}'.format(tau[2]))
ax1.plot( time[1], x[1], color = 'darkgreen', label = 'τ = {0}'.format(tau[1]))


ax1.legend(fontsize = 14)
ax1.set_xlabel('time', fontsize = 16)
ax1.set_ylabel('V(x)', fontsize = 16)


ax2.plot( time[0], x[0], color = 'darkblue', label = 'τ = {0}'.format(tau[0]))
ax2.plot( time[3], x[3], color = 'orange', label = 'τ = {0}'.format(tau[3]))
ax2.plot( time[2], x[2], color = 'firebrick', label = 'τ = {0}'.format(tau[2]))
ax2.plot( time[1], x[1], color = 'darkgreen', label = 'τ = {0}'.format(tau[1]))


ax2.legend(fontsize = 14)
ax2.set_xlabel('time', fontsize = 16)
ax2.set_ylabel('V(x)', fontsize = 16)
ax2.set_ylim(1.6164,1.6172)
ax2.set_xlim(9.815,9.818)
plt.savefig(f"figures/V(x)_timeev_vartau.pdf")
plt.show()

#-------------------------------------------
#Ekin
    
fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(16,6))
fig.text(.5, .93, 'τ → 0, N = {0}, $Ɛ^2$ = {1}'.format(N, eps2), fontsize=18, ha = 'center')
fig.tight_layout(pad=6.0)

ax1.set_title('(a)', fontsize=16, ha = 'right')
ax2.set_title('(b)', fontsize=16, ha = 'right')
ax1.tick_params(axis= 'both', labelsize= 14)
ax2.tick_params(axis= 'both', labelsize= 14)

ax1.plot( time[0], ekin[0], color = 'darkblue', label = 'τ = {0}'.format(tau[0]))
ax1.plot( time[3], ekin[3], color = 'orange', label = 'τ = {0}'.format(tau[3]))
ax1.plot( time[2], ekin[2], color = 'firebrick', label = 'τ = {0}'.format(tau[2]))
ax1.plot( time[1], ekin[1], color = 'darkgreen', label = 'τ = {0}'.format(tau[1]))


ax1.legend(fontsize = 14)
ax1.set_xlabel('time', fontsize = 16)
ax1.set_ylabel('E$_{kin}$', fontsize = 16)



ax2.plot( time[0], ekin[0], color = 'darkblue', label = 'τ = {0}'.format(tau[0]))
ax2.plot( time[3], ekin[3], color = 'orange', label = 'τ = {0}'.format(tau[3]))
ax2.plot( time[2], ekin[2], color = 'firebrick', label = 'τ = {0}'.format(tau[2]))
ax2.plot( time[1], ekin[1], color = 'darkgreen', label = 'τ = {0}'.format(tau[1]))


ax2.legend(fontsize = 14)
ax2.set_xlabel('time', fontsize = 16)
ax2.set_ylabel('E$_{kin}$', fontsize = 16)
ax2.set_xlim(9.815,9.818)
ax2.set_ylim(0.42425, 0.4254)
plt.savefig(f"figures/E_kin_timeev_vartau.pdf")
plt.show()

#------------------------------------------------------------------
#variing N:
#------------------------------------------------------------------

files = glob(f"data/N*")
eps2 = 1.6e-05
tau = 0.01
N = []


ekin = []
x = []
time = []

for i,file in enumerate(files): 
    N.append(getN(file))
    ekin.append(np.genfromtxt(file, dtype=float, usecols=2))
    x.append(np.genfromtxt(file, dtype=float, usecols=1))
    time.append(np.genfromtxt(file, dtype=float, usecols=0))


diff_N1200N1000 = np.abs(x[3][-1] - x[1][-1])
diff_N1200N800 = np.abs(x[3][-1] - x[2][-1])
diff_N1200N600 = np.abs(x[3][-1] - x[0][-1])

print('V(x)\nVariing N differences are: \nN1200N1000 = {} \nN1200N600 = {} \n0.001-0.1 = {}'.format(diff_N1200N1000, diff_N1200N800, diff_N1200N600))

diff_N1200N1000 = np.abs(ekin[3][-1] - ekin[1][-1])
diff_N1200N800 = np.abs(ekin[3][-1] - ekin[2][-1])
diff_N1200N600 = np.abs(ekin[3][-1] - ekin[0][-1])

print('Ekin\nVariing N differences are: \nN1200N1000 = {} \nN1200N600 = {} \n0.001-0.1 = {}'.format(diff_N1200N1000, diff_N1200N800, diff_N1200N600))
#-------------------------------------------
#x
fig, ((ax1, ax2, ax3)) = plt.subplots(1, 3, figsize=(24, 6))
fig.text(.5, .93, 'N → ∞, τ = {0}, $Ɛ^2$ = {1}'.format(tau, eps2), fontsize=21, ha = 'center')
fig.tight_layout(pad=6.0)


ax1.set_title('(a)', fontsize=19, ha = 'right')
ax2.set_title('(b)', fontsize=19, ha = 'right')
ax3.set_title('(c)', fontsize=19, ha = 'right')
ax1.tick_params(axis= 'both', labelsize= 17)
ax2.tick_params(axis= 'both', labelsize= 17)
ax3.tick_params(axis= 'both', labelsize= 17)

ax1.plot( time[0], x[0], color = 'darkgreen', label = 'N = {0}'.format(N[0]))
ax1.plot( time[2], x[2], color = 'firebrick', label = 'N = {0}'.format(N[2]))
ax1.plot( time[1], x[1], color = 'orange', label = 'N = {0}'.format(N[1]))
ax1.plot( time[3], x[3], color = 'darkblue',  linestyle = 'dotted',label = 'N = {0}'.format(N[3]))


ax1.legend(fontsize = 17)
ax1.set_xlabel('time', fontsize = 19)
ax1.set_ylabel('V(x)', fontsize = 19)


ax2.plot( time[0], x[0], color = 'darkgreen', label = 'N = {0}'.format(N[0]))
ax2.plot( time[2], x[2], color = 'firebrick', label = 'N = {0}'.format(N[2]))
ax2.plot( time[1], x[1], color = 'orange', label = 'N = {0}'.format(N[1]))
ax2.plot( time[3], x[3], color = 'darkblue',  linestyle = 'dotted',label = 'N = {0}'.format(N[3]))


ax2.legend(fontsize = 17)
ax2.set_xlabel('time', fontsize = 19)
ax2.set_ylabel('V(x)', fontsize = 19)
ax2.set_ylim(1.42,1.63)
ax2.set_xlim(9.815,9.823)


ax3.plot( time[0], x[0], color = 'darkgreen', label = 'N = {0}'.format(N[0]))
ax3.plot( time[2], x[2], color = 'firebrick', label = 'N = {0}'.format(N[2]))
ax3.plot( time[1], x[1], color = 'orange', label = 'N = {0}'.format(N[1]))
ax3.plot( time[3], x[3], color = 'darkblue',  linestyle = 'dotted',label = 'N = {0}'.format(N[3]))


ax3.legend(fontsize = 17)
ax3.set_xlabel('time', fontsize = 19)
ax3.set_ylabel('V(x)', fontsize = 19)
ax3.set_ylim(1.616,1.618)
ax3.set_xlim(9.815,9.823)



plt.savefig(f"figures/V(x)_timeev_varN.pdf")
plt.show()

#-------------------------------------------
#Ekin
    
fig, ((ax1, ax2, ax3)) = plt.subplots(1, 3, figsize=(24,6))
fig.text(.5, .93, 'N → ∞, τ = {0}, $Ɛ^2$ = {1}'.format(tau, eps2), fontsize=21, ha = 'center')
fig.tight_layout(pad=6.0)

ax1.set_title('(a)', fontsize = 19, ha = 'right')
ax2.set_title('(b)', fontsize = 19, ha = 'right') 
ax3.set_title('(c)', fontsize = 19, ha = 'right')
ax1.tick_params(axis= 'both', labelsize= 17)
ax2.tick_params(axis= 'both', labelsize= 17)
ax3.tick_params(axis= 'both', labelsize= 17)


ax1.plot( time[0], ekin[0], color = 'darkgreen', label = 'N = {0}'.format(N[0]))
ax1.plot( time[2], ekin[2], color = 'firebrick', label = 'N = {0}'.format(N[2]))
ax1.plot( time[1], ekin[1], color = 'orange', label = 'N = {0}'.format(N[1]))
ax1.plot( time[3], ekin[3], color = 'darkblue',  linestyle = 'dotted', label = 'N = {0}'.format(N[3]))


ax1.legend(fontsize = 17)
ax1.set_xlabel('time', fontsize = 19)
ax1.set_ylabel('E$_{kin}$', fontsize = 19)


ax2.plot( time[0], ekin[0], color = 'darkgreen', label = 'N = {0}'.format(N[0]))
ax2.plot( time[2], ekin[2], color = 'firebrick', label = 'N = {0}'.format(N[2]))
ax2.plot( time[1], ekin[1], color = 'orange', label = 'N = {0}'.format(N[1]))
ax2.plot( time[3], ekin[3], color = 'darkblue',  linestyle = 'dotted',label = 'N = {0}'.format(N[3]))


ax2.legend(fontsize = 17)
ax2.set_xlabel('time', fontsize = 19)
ax2.set_ylabel('E$_{kin}$', fontsize = 19)
ax2.set_xlim(9.815,9.823)
ax2.set_ylim(0.41, 0.625)


ax3.plot( time[0], ekin[0], color = 'darkgreen', label = 'N = {0}'.format(N[0]))
ax3.plot( time[2], ekin[2], color = 'firebrick', label = 'N = {0}'.format(N[2]))
ax3.plot( time[1], ekin[1], color = 'orange', label = 'N = {0}'.format(N[1]))
ax3.plot( time[3], ekin[3], color = 'darkblue',  linestyle = 'dotted',label = 'N = {0}'.format(N[3]))


ax3.legend(fontsize = 17)
ax3.set_xlabel('time', fontsize = 19)
ax3.set_ylabel('E$_{kin}$', fontsize = 19)
ax3.set_xlim(9.815,9.823)
ax3.set_ylim(0.424, 0.426)


plt.savefig(f"figures/E_kin_timeev_varN.pdf")
plt.show()

#------------------------------------------------------------------
#variing eps:
#------------------------------------------------------------------

files = glob(f"data/eps2*")
eps2 = []
N = [32000, 1150,11500, 1000] 
tau = 0.01


ekin = []
x = []
time = []

for i,file in enumerate(files): 
    eps2.append(geteps2(file))
    ekin.append(np.genfromtxt(file, dtype=float, usecols=2))
    x.append(np.genfromtxt(file, dtype=float, usecols=1))
    time.append(np.genfromtxt(file, dtype=float, usecols=0))

diff_e08_e07 = np.abs(x[0][-1] - x[2][-1])
diff_e08_1e05 = np.abs(x[0][-1] - x[1][-1])
diff_e08_4e05 = np.abs(x[0][-1] - x[3][-1])

print('V(x)\nVariing eps differences are: \ne08_e07 = {} \ne08_1e05 = {} \ne08_4e05 = {}'.format(diff_e08_e07, diff_e08_1e05, diff_e08_4e05))

diff_e08_e07 = np.abs(ekin[0][-1] - ekin[2][-1])
diff_e08_1e05 = np.abs(ekin[0][-1] - ekin[1][-1])
diff_e08_4e05 = np.abs(ekin[0][-1] - ekin[3][-1])

print('Ekin\nVariing eps differences are: \ne08_e07 = {} \ne08_1e05 = {} \ne08_4e05 = {}'.format(diff_e08_e07, diff_e08_1e05, diff_e08_4e05))
#-------------------------------------------
# #x
fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(16,6))
fig.text(.5, .93, '$Ɛ^2$ → 0, N → ∞, τ = {0}'.format(tau), fontsize=18, ha = 'center')
fig.tight_layout(pad=6.0)

ax1.set_title('(a)', fontsize=14, ha = 'right')
ax2.set_title('(b)', fontsize=14, ha = 'right')


ax1.plot( time[0], x[0], color = 'darkblue', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[0], N[0]))
ax1.plot( time[2], x[2], color = 'orange', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[2], N[2]))
ax1.plot( time[1], x[1], color = 'firebrick', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[1], N[1]))
ax1.plot( time[3], x[3], color = 'darkgreen', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[3], N[3]))

ax1.legend(fontsize = 14)
ax1.set_xlabel('time', fontsize = 16)
ax1.set_ylabel('V(x)', fontsize = 16)
ax1.tick_params(axis= 'both', labelsize= 14)



ax2.plot( time[0], x[0], color = 'darkblue', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[0], N[0]))
ax2.plot( time[2], x[2], color = 'orange', linestyle = 'dashed', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[2], N[2]))
ax2.plot( time[1], x[1], color = 'firebrick', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[1], N[1]))
ax2.plot( time[3], x[3], color = 'darkgreen', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[3], N[3]))

ax2.legend(fontsize = 14)
ax2.set_xlabel('time', fontsize = 16)
ax2.set_ylabel('V(x)', fontsize = 16)
ax2.set_ylim(1.6161,1.6183)
ax2.set_xlim(9.815,9.823)
ax2.tick_params(axis= 'both', labelsize= 14)




plt.savefig(f"figures/V(x)_timeev_vareps.pdf")
plt.show()

#-------------------------------------------
#Ekin
    
fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(16,6))
fig.text(.5, .93, '$Ɛ^2$ → 0, N → ∞, τ = {0}'.format(tau), fontsize=18, ha = 'center')
fig.tight_layout(pad=6.0)

ax1.set_title('(a)', fontsize=16, ha = 'right')
ax2.set_title('(b)', fontsize=16, ha = 'right')

ax1.plot( time[0], ekin[0], color = 'darkblue', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[0], N[0]))
ax1.plot( time[2], ekin[2], color = 'orange', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[2], N[2]))
ax1.plot( time[1], ekin[1], color = 'firebrick', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[1], N[1]))
ax1.plot( time[3], ekin[3], color = 'darkgreen', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[3], N[3]))

ax1.legend(fontsize = 14)
ax1.set_xlabel('time', fontsize = 16)
ax1.set_ylabel('E$_{kin}$', fontsize = 16)
ax1.tick_params(axis= 'both', labelsize= 14)



ax2.plot( time[0], ekin[0], color = 'darkblue', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[0], N[0]))
ax2.plot( time[2], ekin[2], color = 'orange', linestyle = 'dashed', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[2], N[2]))
ax2.plot( time[1], ekin[1], color = 'firebrick', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[1], N[1]))
ax2.plot( time[3], ekin[3], color = 'darkgreen', label = '$Ɛ^2$ = {0}, N = {1}'.format(eps2[3], N[3]))

ax2.legend(fontsize = 14)
ax2.set_xlabel('time', fontsize = 16)
ax2.set_ylabel('E$_{kin}$', fontsize = 16)
ax2.set_xlim(9.815,9.823)
ax2.set_ylim(0.424, 0.426)
ax2.tick_params(axis= 'both', labelsize= 14)

plt.savefig(f"figures/E_kin_timeev_vareps.pdf")
plt.show()



# %%
i = x[0][0]*eps2[0]

print(i)
# %%
