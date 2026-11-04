#%%
import numpy as np 
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os



name_raw = "evolution_N{}_D{}_mu{:.6f}_epssq{:.6f}_method{}_tau*_steps1"
N = 200
D = 1
mu = 20
eps2 = 0.001
name = name_raw.format(N, D, mu, eps2, "{}")
names = {
    "eu" : name.format(0),
    "cn" : name.format(1),
    "ss" : name.format(2)
}

differ = {
    "eu" : [None],
    "cn" : [None],
    "ss" : [None]
}
taus = {
    "eu" : [0],
    "cn" : [0],
    "ss" : [0]
}

dtype = np.dtype([
    ("norm", float),
    ("e", float),
    ("ekin", float),
    ("epot", float),
    ("pos", float, (D, )),
    ("unc", float, (D, )),
    ("wf", complex, (N**D, ))
])

def converter(x):
    return complex(x.decode())
converters = {4+2*D + i : converter for i in range(N**D)}



filename = "evolution_N200_D1_mu20.000000_epssq0.001000_method1_tau0.100000_steps10000"
data = np.genfromtxt(filename, dtype=dtype, converters=converters)

D = 1
#ani_data_list =  data["wf"]
ani_data_list = [np.abs(i) for i in  data["wf"]]

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure(figsize=(8, 6), dpi=300)
ax = fig.add_subplot(1,1,1)





#print(ani_data_list[1])
# animation function.  This is called sequentially
x_array = np.arange(-N+1,N)
def animate(i):
    x = x_array
    #y = np.sin(2 * np.pi * (x - 0.01 * i))
    y = np.concatenate([np.flip(ani_data_list[i][1:]),ani_data_list[i][:]]) 
    #y = ani_data_list[i][:]
    ax.clear()
    ax.set_ylim([0,0.100])
    ax.set_xlim([-21, 21])
    ax.set_xlabel('n')
    ax.set_ylabel('$\psi(n)$')
    ax.scatter(x, y)
 


# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, func = animate, frames = 10000 ,interval = 10)

# save the animation as an mp4.  
anim.save('timeevolution_animation.mp4', extra_args=['-vcodec', 'libx264'])

plt.show()





# %%

# %%
