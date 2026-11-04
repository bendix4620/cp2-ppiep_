#%%

from operator import itemgetter
import numpy as np
from glob import glob
import matplotlib.pyplot as plt


def getparams(filename):
    str = filename.split("q")
    eps2 = str[1].split(".txt")[0]
 
    return eps2

 
D = 1


#fixed tau = 0.01
ekin_files = glob(f"Ekin*")
x_files = glob(f"x*")

tau = 0.01
eps2 = np.zeros(len(ekin_files))
data_ekin = []
data_x = []
ekin = np.zeros(len(ekin_files))
x = np.zeros(len(ekin_files))
n = np.zeros(len(ekin_files))


#x

for i,file in enumerate(x_files):  
    eps2[i] = getparams(file)
    data_x.append(np.genfromtxt(file, dtype=float, delimiter=" "))

data_x = np.array(data_x)
n = data_x[:,:,0]
x = data_x[:,:,1]
n_index = np.zeros_like(n)
x_sorted = np.zeros_like(x)

for i in range(len(n)):
    n_index[i] = np.argsort(n[i,:])
    n[i]=np.sort(n[i,:])
    for j in range(len(n[i])):
        x_sorted[i,j] = x[i,int(n_index[i][j])]
    n[i][-1] = 31500
    
fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(16,6))
fig.text(.5, .93, '$Ɛ^2$ → 0, N → ∞, τ = 0.01', fontsize=18, ha = 'center')
fig.tight_layout(pad=6.0)
ax1.set_title('(a)', fontsize=16, ha = 'right')
ax2.set_title('(b)', fontsize=16, ha = 'right')
ax1.tick_params(axis= 'both', labelsize= 14)
ax2.tick_params(axis= 'both', labelsize= 14)


ax1.plot( n[0], x_sorted[0], color = 'darkblue', label = '$Ɛ^2$ = {0}'.format(eps2[0]))
ax1.plot( n[1], x_sorted[1],  color = 'orange',label = '$Ɛ^2$ = {0}'.format(eps2[1]))
ax1.plot( n[2], x_sorted[2],  color = 'firebrick', label = '$Ɛ^2$ = {0}'.format(eps2[2]))
ax1.plot( n[3], x_sorted[3], color = 'darkgreen',label = '$Ɛ^2$ = {0}'.format(eps2[3]))
ax1.legend(fontsize=14)
ax1.set_xlabel('N',fontsize=16)
#ax1.set_ylim(-0.00001, 0.0001)
ax1.set_ylabel('V(x)',fontsize=16)


ax2.hlines(y=x_sorted[0,-1], xmin = 31000, xmax = 31500, label = '$Ɛ^2$ = {0}'.format(eps2[0]), color = 'darkblue',)
ax2.hlines(y=x_sorted[1,-1], xmin = 31000, xmax = 31500, label = '$Ɛ^2$ = {0}'.format(eps2[1]), color = 'orange')
ax2.hlines(y=x_sorted[2,-1],  xmin = 31000, xmax = 31500, label = '$Ɛ^2$ = {0}'.format(eps2[2]), color = 'firebrick')
ax2.hlines(y=x_sorted[3,-1],  xmin = 31000, xmax = 31500, label = '$Ɛ^2$ = {0}'.format(eps2[3]), color = 'darkgreen')
ax2.legend(fontsize=14)
ax2.set_xlabel('N', fontsize=16)
ax2.set_ylabel('V(x)', fontsize=16)
plt.savefig(f"figures/fixed_tau_x.pdf")
plt.show()

#-----------------------------------------------------------
#kinetic energy 


n = np.zeros(len(ekin_files))

for i,file in enumerate(ekin_files):  
    eps2[i] = getparams(file)
    #data[i] = np.genfromtxt(file, dtype=float, delimiter=" ")
    data_ekin.append(np.genfromtxt(file, dtype=float, delimiter=" "))

data_ekin = np.array(data_ekin)

n = data_ekin[:,:,0]
ekin = data_ekin[:,:,1]
n_index = np.zeros_like(n)
ekin_sorted = np.zeros_like(ekin)


for i in range(len(n)):
    n_index[i] = np.argsort(n[i,:])
    n[i]=np.sort(n[i,:])
    for j in range(len(n[i])):
        ekin_sorted[i,j] = ekin[i,int(n_index[i][j])]
    n[i][-1] = 31500


fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(16,6))
fig.text(.5, .93, '$Ɛ^2$ → 0, N → ∞, τ = 0.01', fontsize=18, ha = 'center')
fig.tight_layout(pad=6.0)
ax1.set_title('(a)', fontsize=16, ha = 'right')
ax2.set_title('(b)', fontsize=16, ha = 'right')
ax1.tick_params(axis= 'both', labelsize= 14)
ax2.tick_params(axis= 'both', labelsize= 14)


ax1.plot( n[0], ekin_sorted[0], color = 'darkblue', label = '$Ɛ^2$ = {0}'.format(eps2[0]))
ax1.plot( n[1], ekin_sorted[1], color = 'orange', label = '$Ɛ^2$ = {0}'.format(eps2[1]))
ax1.plot( n[2], ekin_sorted[2], color = 'firebrick', label = '$Ɛ^2$ = {0}'.format(eps2[2]))
ax1.plot( n[3], ekin_sorted[3], color = 'darkgreen', label = '$Ɛ^2$ = {0}'.format(eps2[3]))

ax1.legend(fontsize=14)
ax1.set_xlabel('N',fontsize=16)
ax1.set_ylabel('E$_{kin}$',fontsize=16)
ax1.set_ylim(0.3,0.45)




ax2.hlines(y=ekin_sorted[0,-1], xmin = 31000, xmax = 31500, label = '$Ɛ^2$ = {0}'.format(eps2[0]), color = 'darkblue')
ax2.hlines(y=ekin_sorted[1,-1], xmin = 31000, xmax = 31500, label = '$Ɛ^2$ = {0}'.format(eps2[1]), color = 'orange')
ax2.hlines(y=ekin_sorted[2,-1],  xmin = 31000, xmax = 31500, label = '$Ɛ^2$ = {0}'.format(eps2[2]), color = 'firebrick')
ax2.hlines(y=ekin_sorted[3,-1],  xmin = 31000, xmax = 31500, label = '$Ɛ^2$ = {0}'.format(eps2[3]), color = 'darkgreen')
ax2.legend(fontsize=14)
ax2.set_xlabel('N',fontsize=16)
ax2.set_ylabel('E$_{kin}$',fontsize=16)
plt.savefig(f"figures/fixed_tau_Ekin.pdf")
plt.show()



# %%

print(x_sorted[0,0]*eps2[0])
# %%
