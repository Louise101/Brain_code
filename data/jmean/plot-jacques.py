# !/home/lewis/anaconda3/bin/python
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pylab as plt
import numpy as np
from scipy.optimize import curve_fit


def getR2(xdata, ydata, f, *opts):
    """https://en.wikipedia.org/wiki/Coefficient_of_determination
    """
    residuals = ydata - f(xdata, *opts)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata - np.mean(ydata))**2)
    r_squared = 1. - (ss_res / ss_tot)
    return r_squared


def jacfunc(x, c1, c2, k1, k2):
    return c1 * np.exp(-(k1 * x) / delta) - c2 * np.exp(-(k2 * x) / delta)


width = 6.510  # inches
height = width / 1.618

plt.rc('font', family='serif', serif='Times')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('axes', labelsize=12)


# setup grid
fig, ax = plt.subplots()
fig.subplots_adjust(left=.08, bottom=.12, right=.95, top=.965)

file = "data/jmean/validation-630.dat"

delta = 0.261
xs = np.linspace(0, 2., 1000)
x, y = np.loadtxt(file, unpack=True)

ax.scatter(x[::6], y[::6], marker="x", label="MCRT data - 630nm")
popt, pconv = curve_fit(jacfunc, x, y, p0=[6.27, 1.18, 1., 14.4])
r2 = getR2(x, y, jacfunc, *popt)
ax.plot(xs, jacfunc(xs, *popt), label=r"MCRT fit $r^2$={:07.5f} - 630nm".format(r2), color="red", linestyle=":")
ax.plot(xs, jacfunc(xs, 6.27, 1.18, 1., 14.4), label="Jacques - 630nm", color="green")

print(*popt)
print(6.27, 1.18, 1., 14.4)


file = "data/jmean/validation-420.dat"

delta = 0.047
xs = np.linspace(0, 2., 1000)
x, y = np.loadtxt(file, unpack=True)

ax.scatter(x[::3], y[::3], marker="x", label="MCRT data - 420nm")
popt, pconv = curve_fit(jacfunc, x, y, p0=[5.76, 1.31, 1., 10.2])
r2 = getR2(x, y, jacfunc, *popt)
ax.plot(xs, jacfunc(xs, *popt), label=r"MCRT fit $r^2$={:07.5f} - 420nm".format(r2), color="red", linestyle=":")
ax.plot(xs, jacfunc(xs, 5.76, 1.31, 1., 10.2), label="Jacques - 420nm", color="blue")

print(*popt)
print(5.76, 1.31, 1., 10.2)

ax.set_xlabel("Depth/cm")
ax.set_ylabel(r"Fluence rate $\Psi/\Psi_0$ [-]")
ax.set_xlim(-.01, 1.)
ax.set_ylim(10**-9,)
ax.legend()
fig.set_size_inches(width, height)
fig.savefig('validation-jacques-test.pdf')
plt.show()
