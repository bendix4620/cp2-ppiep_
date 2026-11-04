import numpy as np
import matplotlib.pyplot as plt
import ast

content = np.genfromtxt(
    "../build/groundstate_N200_D1_mu20.000000_epssq0.001000",
    dtype = np.complex128, comments = "#", converters={0: lambda s: complex(s.decode().replace('+-', '-'))})

E = np.abs(content[0])
psi = content[1:]

def pot(n):
    mu = 20
    eps2 = 0.001
    return mu/8 * (eps2 * (n)**2 - 1 )**2

n = np.arange(len(psi)) - len(psi)/2
plt.figure()
plt.plot(n, np.abs(psi), label="wavefunction, E = {0:.3f}".format(E))
# plt.hlines(E, n[0], n[-1], label="energy", color="orange")
p = pot(n)
plt.plot(n, p*(1.3 * np.max(psi)/np.max(p)), label="scaled potential")

# plt.ylim(0, 1)
plt.legend(loc="upper center")
plt.xlabel("position n")
plt.ylabel("absolute value")
plt.savefig("groundstate.pdf", bbox_inches="tight", format="pdf")
