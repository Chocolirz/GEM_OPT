import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import numpy as np
from scipy.optimize import curve_fit
from tqdm import tqdm

params = params = {
    'axes.labelsize': 21,
    'font.size': 20,
    'font.family': 'sans-serif', 
    'font.serif': 'Arial', 
    'legend.fontsize': 18,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18, 
    'axes.labelpad': 15,
    
    'figure.figsize': [10,8], # value in inches based on dpi of monitor
    'figure.dpi': 105.5, # My monitor has a dpi of around 105.5 px/inch
   
    'axes.grid': False,
    'grid.linestyle': '-',
    'grid.alpha': 0.25,
    'axes.linewidth': 1,
    'figure.constrained_layout.use': True,
   
   
    # Using Paul Tol's notes:
    'axes.prop_cycle': 
        mpl.cycler(color=['#4477aa', # blue
                          '#ee6677', # red/pink
                          '#228833', # green
                          '#aa3377', # purple
                          '#66ccee', # cyan
                          '#ccbb44', # yellow
                          '#bbbbbb', # grey
                          'k']),
       
       # Pick either the cycler above, or the cycler below:
        
       # (mpl.cycler(color=['#4477aa', # blue
       #                     '#ee6677', # red/pink
       #                     '#228833', # green
       #                     '#aa3377', # purple
       #                     '#66ccee', # cyan
       #                     '#ccbb44', # yellow
       #                     '#bbbbbb', # grey
       #                     ]) + 
       #   mpl.cycler(linestyle=['-', # solid
       #                         '--', # dashed
       #                         ':', # dotted
       #                         '-.', # dash dot
       #                         (0, (3, 1, 1, 1, 1, 1)), # narrow dash dot dot
       #                         (0, (1, 2, 5, 2, 5, 2)), # dash dash dot
       #                         (0, (5, 2.5, 1, 2.5, 1, 2.5)), # dash dot dot
       #                         ])), 
       
    'lines.linewidth': 2.5,
   
    'image.cmap': 'turbo',#'PuBuGn',
   
    'legend.frameon': False,
   
    'legend.edgecolor': 'white',
   
    'legend.framealpha': 0.5,
} 

plt.rcParams.update(params)

log_tick_format = ticker.FuncFormatter(lambda y,pos: f'{y:.{int(np.maximum(-np.log10(y),0)):1d}f}')

def fitFunction(x, A, B, C, D):
    return A*x**6 + B*x**4 + C*x**2 + D



plt.figure(figsize=(8,8))
plt.xlabel(r'$R$ [um]')
plt.ylabel(r'# of electrons per avalanche')
plt.yscale('log')
plt.xlim(0, 110)

data = np.loadtxt('rdist_triple_LH_corrected.txt', skiprows=1)
bins_left = data[:,0]
bins_right = data[:,1]
bins_width = bins_left[0] - bins_right[0]
entries = data[:,2]

bins_centre = []
for i in range(len(bins_left)):
    bins_centre.append((bins_left[i] + bins_right[i]) * 1e4 / 2)

bins_centre_new = np.append(bins_centre, 105)
sigma = np.ones_like(entries) * 10
sigma = np.append(sigma, 1e-6)
entries_new = np.append(entries, 0)

#popt, pcov = curve_fit(fitFunction, bins_centre, entries, p0 = (-1e-6, -5e-6, -4e-2, -6, 1e6), sigma = sigma)
popt, pcov = curve_fit(fitFunction, bins_centre, entries, p0 = (-5e-6, -4e-2, -6, 1e6))

print(popt)

x_plot = np.linspace(0, 210 / 2, 1000)
y_plot = popt[0]*x_plot**6 + popt[1]*x_plot**4 + popt[2]*x_plot**2 + popt[3]

plt.plot(bins_centre, entries / 1000, 'o', color = 'black', label = 'Data, random')
plt.plot(x_plot, y_plot / 1000, color = 'orange', label = 'Fitted, random')

data = np.loadtxt('rzdist_triple_LH_corrected.txt', skiprows=1)
bins_left = data[:,0]
bins_right = data[:,1]
bins_width = bins_left[0] - bins_right[0]
entries = data[:,2]

bins_centre = []
for i in range(len(bins_left)):
    bins_centre.append((bins_left[i] + bins_right[i]) * 1e4 / 2)

bins_centre_new = np.append(bins_centre, 105)
sigma = np.ones_like(entries) * 10
sigma = np.append(sigma, 1e-6)
entries_new = np.append(entries, 0)

#popt, pcov = curve_fit(fitFunction, bins_centre, entries, p0 = (-1e-6, -5e-6, -4e-2, -6, 1e6), sigma = sigma)
popt, pcov = curve_fit(fitFunction, bins_centre, entries, p0 = (-5e-6, -4e-2, -6, 1e6))

print(popt)

x_plot = np.linspace(0, 210 / 2, 1000)
ratio = y_plot / (popt[0]*x_plot**6 + popt[1]*x_plot**4 + popt[2]*x_plot**2 + popt[3])
y_plot = popt[0]*x_plot**6 + popt[1]*x_plot**4 + popt[2]*x_plot**2 + popt[3]

plt.plot(bins_centre, entries / 1000, 'o', color = 'gray', label = 'Data, centre')
plt.plot(x_plot, y_plot / 1000, color = 'blue', label = 'Fitted, centre')


plt.grid(True, which='major')
plt.grid(True, which='minor')
plt.legend()
plt.savefig('rdist_triple_LH_corrected_mixed.png', dpi = 300)
plt.show()

# -1.77470958e-07  1.72758442e-03 -2.62086413e+01  3.77759747e+05
# normalisation constant 3.03939e+07
# CDF easily computed, however,
# there is no algebraic formula for the inverse of general polynomials of degree 5 or higher (Abel-Ruffini theorem).
# Use numerical root finding
# For N electrons, divide r into 1000 parts, dP = 0.5 * [f(r) + f(r+dr)] dr
# n(r) = N * dP
# for n(r) times, write r. 
'''
r_list = np.linspace(0, 105, 106)
f_r = fitFunction(r_list, *popt) / 3.03939e+07
dr = 210 / (2*106)
with open(f'Ne_r_single.txt', 'w') as file:
    for i in tqdm(range(len(r_list))):
        n_r = dr * f_r[i]
        file.write(str(r_list[i]) + '   ' + str(n_r) +'\n')
'''

plt.figure(figsize=(8,8))
plt.xlabel('Radius [um]')
plt.ylabel('Random / Centre')
plt.plot(x_plot, ratio, color = '#4477aa')
plt.hlines(y=np.average(ratio), xmin = 0, xmax = 105, linestyle='--', label=f'Average {round(np.average(ratio),4)}', color = '#ee6677')
plt.grid(True)
plt.legend()
plt.savefig('./rdist_triple_LH_corrected_ratio.png', dpi = 300)
plt.show()