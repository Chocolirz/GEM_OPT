import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import numpy as np
from scipy.optimize import curve_fit

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

data = np.loadtxt('./PhotoGain.txt', skiprows=1)
diameter = data[:,0]
voltages = data[:,1]
NoPP = data[:,2]
NoPG = data[:,3]
Gain = data[:,4]

def Theory(x, A, B, C):
    return A*np.exp(x*B + C)

def Fit(x, y, Theory):
    popt, pcov = curve_fit(Theory, x, y, p0 = [100, 0.05, -5])
    return popt[0], popt[1], popt[2]

PlotPhotoGain = False
PlotPhE = True

if PlotPhotoGain:
    #A, B, C = Fit(diameter, NoPG, Theory)
    x = np.linspace(0, 1000, 100)

    plt.figure(figsize=(8,8))
    plt.title('Gain of Photons')
    plt.xlabel('Diameter [um]')
    plt.ylabel('Gain [photons]')
    #plt.yscale('log')
    #plt.plot(x, Theory(x, A, B, C), color = 'red', label=f'Fitted $Ae^(Bx+C)$ \n A={A}, B={B}, C={C}')
    plt.plot(diameter, NoPG, 'o', color = 'blue', label='Simulated')
    plt.legend()
    plt.savefig('./img/PhotoGain.png', dpi = 300)
    plt.show()

if PlotPhE:
    plt.figure(figsize=(8,8))
    plt.xlabel('Gain [#Electrons]')
    plt.ylabel('Gain [#Photons]')
    #plt.yscale('log')
    #plt.plot(x, Theory(x, A, B, C), color = 'red', label=f'Fitted $Ae^(Bx+C)$ \n A={A}, B={B}, C={C}')
    plt.plot(Gain, NoPG, 'o', color = 'blue', label='Simulated')
    for i in range(len(Gain)):
        plt.text(Gain[i] + 2, NoPG[i] + 2, f'{int(diameter[i])}'+r'${\rm \mu m}$', fontsize = 'x-small')
    plt.xlim((110,450))
    plt.legend()
    plt.savefig('./img/PhE.png', dpi = 300)
    plt.show()