import os
os.chdir(os.path.dirname(__file__))
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit


N = 5000
D = 1
E0 = 4.866819e-1
tau = "0.100000"
total = "100000"

def get_array(filename, **kwargs):
    return np.genfromtxt(filename, dtype=complex, **kwargs)

data = get_array(f"data/melement_tau{tau}_total{total}.txt", max_rows=None)
plt.plot(np.arange(1000)*float(tau), np.abs(data[:1000]))
plt.ylabel("$M(k \\cdot \\tau)$")
plt.xlabel("$k \\cdot \\tau$")
# plt.savefig(f"data/M_tau{tau}.pdf", format="pdf", bbox_inches="tight")
plt.show()

f = np.abs(np.fft.ifft(data))
x = np.fft.fftfreq(len(data), float(tau))

plt.plot(x, f)
plt.xlim(-1, 2)
plt.ylim(0, 100)
plt.ylabel("$M(\\omega)$")
plt.xlabel("$\\omega$")
# plt.savefig(f"data/peaks_overview_tau{tau}_total{total}.pdf", format="pdf", bbox_inches="tight")
plt.show()

p, _ = find_peaks(f, height=0.0001, prominence=0.000001, wlen=20)
xp = x[p]
print(xp, len(xp))
fp = f[p]
plt.plot(x, f)
for i in range(len(p)):
    plt.plot([xp[i]]*2, [0, fp[i]], color="red")
# plt.plot(xp, fp, "rx")
plt.ylim(3.3e-4, 1)
plt.yscale("log")
plt.xlim(-0.1, 1.5)
plt.ylabel("$M(\\omega)$")
plt.xlabel("$\\omega$")
# plt.savefig(f"data/peaks_tau{tau}_total{total}.pdf", format="pdf", bbox_inches="tight")
plt.show()

plt.plot([E0, *(xp + E0)], "x")
plt.xlabel("n")
plt.ylabel("$E_n$")
# plt.savefig(f"data/energy_levels_tau{tau}_total{total}.pdf", format="pdf", bbox_inches="tight")
# np.savetxt(f"data/energy_levels_tau{tau}_total{total}.txt", [E0, *(xp + E0)])
plt.show()
