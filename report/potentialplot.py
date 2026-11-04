#%%
import matplotlib.pyplot as plt
import numpy as np 

epslow = 0.001
epshigh = 0.002
mu = 20

def f(x, eps2):
    f = mu/8 * (eps2 * x**2 - 1 )**2 
    return f

x_ar = np.arange(-50,51)
y_ar1 = []
y_ar2 = []

for i in range(-50,51):
    y_ar1.append(f(i, epslow))
    y_ar2.append(f(i, epshigh))

plt.figure(figsize=(7, 5))
plt.plot(x_ar, y_ar1, label = '$\epsilon$ = 0.001')
plt.plot(x_ar, y_ar2, label = '$\epsilon$ = 0.002')
plt.ylim(0,5)
plt.legend()
plt.xlabel('x')
plt.ylabel('V(x) = $\mu / 8 (\epsilon^2$ x$^2 - 1)^2 $')
plt.show()
# %%

# %%
