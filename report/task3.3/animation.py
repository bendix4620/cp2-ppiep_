import os
os.chdir(os.path.dirname(__file__))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


N = 12000
D = 1

def convert(x):
    return complex(x.decode())

def get_array(filename, **kwargs):
    # dtype = np.dtype([
    #     ("norm", float),
    #     ("E", float),
    #     ("K", float),
    #     ("V", float),
    #     ("x", float, (D, )),
    #     ("x2", float, (D, )),
    #     ("wave", complex, (N**D, ))
    # ])

    cols = np.arange(4+2*D, 4+2*D + N**D, 1)
    conv = {c : convert for c in cols}

    return np.abs(np.genfromtxt(
        filename, usecols=cols, converters=conv, dtype=complex, **kwargs
    ))


def animate(frame):
    ax.clear()
    ax.plot(frame)
    ax.set_ylim(0, max)
    return

data = get_array("data/N12000_epssp1e-7.txt", max_rows=None, skip_header=0)
max = np.max(data)
fig, ax = plt.subplots()
anim = animation.FuncAnimation(fig, func=animate, frames=data[::100])
anim.save('animation_big.mp4', fps=144)
