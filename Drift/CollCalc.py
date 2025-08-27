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

PlotCollRatio = True

if PlotCollRatio:
    data = np.loadtxt('collection.txt', skiprows=1)
    voltage = data[:,0]
    radius = data[:,1]
    z_start = data[:,2]
    coll_ratio = data[:,3]
    print(np.average(coll_ratio[1:5]))

    plt.figure(figsize=(8,8))
    plt.xlabel('Starting $z$ [mm]')
    plt.ylabel('Probability of collection [%]')
    plt.xlim(1,9)
    plt.ylim(97.75,98.75)
    plt.plot(z_start[:5], coll_ratio[:5], 'o', color = 'red', label = '170 um')
    plt.hlines(y = np.average(coll_ratio[1:5]), xmin = 1, xmax = 9, linestyle='--', color = 'red')
    plt.plot(z_start[5:], coll_ratio[5:], 'o', color = 'blue', label = '210 um')
    plt.hlines(y = np.average(coll_ratio[6:]), xmin = 1, xmax = 9, linestyle='--', color = 'blue')
    plt.grid(True)
    plt.legend()
    plt.savefig('./img/Collection_Ratio.png', dpi = 300)
    plt.show()