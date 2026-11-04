import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit

os.chdir(os.path.dirname(__file__))
print(os.getcwd())
files = glob("tau2zero/evolution_v2_N*")

N=200; D=1

def getmethod(x):
    _, x = x.split("_method")
    x = x.split("_")[0]
    if x == "0":
        return "eu"
    elif x == "1":
        return "cn"
    elif x == "2":
        return "ss"
    else:
        raise ValueError("Could not understand method={}".format(x))

def gettau(x):
    _, x = x.split("_tau")
    x = x.split("_")[0]
    return float(x)

def getsteps(x):
    _, x = x.split("_steps")
    x = x.split("_")[0]
    return int(x)

states = {
    "eu" : [],
    "cn" : [],
    "ss" : []
}
taus = {
    "eu" : [],
    "cn" : [],
    "ss" : []
}

dtype = np.dtype([
    ("norm", float),
    ("e", float),
    ("ekin", float),
    ("epot", float),
    ("pos", float, (D, )),
    ("upos", float, (D, )),
    ("mom", float, (D, )),
    ("umom", float, (D, )),
    ("wf", complex, (N**D, ))
])

def converter(x):
    return complex(x.decode())
converters = {4+4*D + i : converter for i in range(N**D)}

for filename in files:
    print("Processing", filename)
    method = getmethod(filename)
    tau = gettau(filename)

    data = np.genfromtxt(filename, converters=converters, dtype=dtype)

    states[method].append(data["wf"][-1, :])

    taus[method].append(tau)

def f(x, a, b, c):
    return b*(np.exp(a*x)) + c

fig, ax = plt.subplots()

sort = {}
for method in ("eu", "cn", "ss"):
    sort[method] = np.argsort(taus[method])

for m0, m1 in ( ("eu", "cn"), ("eu", "ss"), ("cn", "ss") ):
    s0 = sort[m0]
    w0 = np.array(states[m0])[s0]

    s1 = sort[m1]
    w1 = np.array(states[m1])[s1]

    t = np.array(taus[m0])[s0]
    v = np.linalg.norm(w0 - w1, axis=1)
    np.savetxt(f"tau2zero/data_{m0}-{m1}.txt", np.transpose([t, v]))
    line, = plt.plot(t, v, "x", label = f"{m0}-{m1}")

    # popt, pcov = fit(f, t, v, maxfev=10000)
    # plt.plot(t, f(t, *popt), color=line.get_color(), linewidth = 1)

ax.invert_xaxis()
plt.xlabel("time step $\\tau$")
plt.ylabel("difference in wavefunctions between 2 integrators")
plt.legend()
plt.xscale("log")
# plt.yscale("log")
plt.savefig("tau2zero/01tau2zero.pdf", format="pdf", bbox_inches="tight")
plt.close()
