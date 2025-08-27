import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import numpy as np

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

PlotTC = False
if PlotTC:
    # This is for pure CF4
    p = [50, 60, 70, 80, 90, 100] # in Torr
    A = [19.2, 17.6, 16.0, 14.5, 13.2, 11.6] # cm^-1 Torr^-1
    A_err = []
    B = [465, 435, 405, 377, 356, 325] # V cm^-1 Torr^-1
    B_err = []

    # Enter field or voltage
    V = 560 # V
    d = 0.057 # cm
    E = V / d

    # Calculate the coefficient
    alpha = []
    gain = []
    aop = []
    Eop = []
    for i in range(len(A)):
        alpha.append(p[i] * A[i] * np.exp(-B[i] * p[i] / E))
        gain.append(np.exp(alpha[i] * d))
        aop.append(alpha[i]/p[i])
        Eop.append(E/p[i])

    plt.figure(figsize=(8,8))
    plt.xlabel(r'Pressure [Torr]')
    plt.ylabel(r'$\alpha$ [cm$^{-1}$]')
    plt.plot(p, alpha, 'o', color = 'red')
    plt.show()


    plt.figure(figsize=(10,8))
    plt.title(f'Gain as a function of pressure')
    plt.xlabel(r'Pressure [Torr]')
    plt.ylabel(r'Gain')
    plt.plot(p, gain, 'o', color = 'red')
    plt.show()

    plt.figure(figsize=(10,8))
    plt.title(r'$\alpha/p$ vs $E/p$ plot')
    plt.xlabel(r'$E/p$ [V/(cm Torr)]')
    plt.ylabel(r'$\alpha/p$ [cm^{-1} Torr^{-1}]')
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(Eop, aop, 'o', color = 'red')
    plt.show()

PlotCS = True
if PlotCS:
    p = 6666.12 # 50 Torr, pressure
    k = 1.3806488e-23 # Boltzmann's constant
    T = 293.15 # K, temperature
    mfp_ion = [] # mean free path
    mfp_tot = []

    data = np.loadtxt('./build/cs.txt', encoding="cp1252", skiprows=61)
    colors = [plt.cm.rainbow(i / data.shape[1]) for i in range(data.shape[1])]

    plt.figure(figsize=(12,8))
    plt.xlabel(r"$E_{\rm electrons}$ [eV]")
    plt.ylabel(r"$\sigma$ of different types [cm$^2$]")
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim((1.5e-2,100))
    
    for col in range(1, data.shape[1]):
        x = data[:, 0]
        y = data[:, -col]
        plt.plot(x, y, color = colors[col])

    plt.grid(True)
    plt.savefig('./img/ElectronCrossSection.png', dpi=300)
    plt.show()