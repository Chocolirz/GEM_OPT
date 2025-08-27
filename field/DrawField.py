import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import numpy as np
import matplotlib.patches as mpatches

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

PlotFieldr = False
PlotFieldz = True
PlotFieldEdge = False
PlotFieldD = False

def moving_average(y, window_size):
    return np.convolve(y, np.ones(window_size)/window_size, mode='same')

if PlotFieldz:
    plt.figure(figsize=(8,8))
    for i in range(54, 61):
        v = 10*i
        data = np.loadtxt(f'efield_z_{v}V_210um.txt', skiprows=1)
        z = data[:, 0]
        Ex = data[:, 1]
        Ey = data[:, 2]
        Ez = data[:, 3]
        plt.plot(z, Ez, label=f'{v} V')
    plt.xlabel('z [cm]')
    plt.ylabel('E [V/cm]')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.savefig('./img/eFieldz.png', dpi = 300)
    plt.show()

if PlotFieldr:
    plt.figure(figsize=(8,8))
    for i in range(15, 26):
        v = 10*i
        data = np.loadtxt(f'efield_r_{v}um.txt', skiprows=1)
        r = data[:, 0]
        Ex = data[:, 1]
        Ey = data[:, 2]
        Ez = data[:, 3]
        plt.plot(r, Ez, label=f'{v} um')
    plt.xlabel('r [cm]')
    plt.ylabel('E [V/cm]')
    plt.yscale('linear')
    plt.grid(True)
    plt.legend()
    plt.savefig('./img/eFieldr.png', dpi = 300)
    plt.show()

if PlotFieldEdge:
    plt.figure(figsize=(8,8))
    for i in range(15, 26):
        v = 10*i
        data = np.loadtxt(f'efield_z_edge_{v}um.txt', skiprows=1)
        z = data[:, 0]
        Ex = data[:, 1]
        Ey = data[:, 2]
        Ez = data[:, 3]
        Emag = []
        for i in range(len(z)):
            Emag.append(np.sqrt(Ex[i]**2 + Ey[i]**2 + Ez[i]**2))
        plt.plot(z, Emag, label=f'{v} um')
    plt.xlabel('z [cm]')
    plt.ylabel('E [V/cm]')
    plt.grid(True)
    plt.legend()
    plt.savefig('./img/eFieldzEdge.png', dpi = 300)
    plt.show()

if PlotFieldD:
    data_above = np.loadtxt('efieldAbove.txt', skiprows=1)
    diameters = data_above[:,0]
    field_above = data_above[:,1]
    data_below = np.loadtxt('eFieldBelow.txt', skiprows=1)
    field_below = data_below[:,1]
    plt.figure(figsize=(8,8))
    plt.xlabel('Diameter [um]')
    plt.ylabel('Drift Field [V/cm]')
    plt.plot(diameters, field_above, 'o')
    plt.grid(True)
    plt.savefig('./img/eFieldAbove.png', dpi = 300)
    plt.show()

    plt.figure(figsize=(8,8))
    plt.xlabel('Diameter [um]')
    plt.ylabel('Induction Field [V/cm]')
    plt.plot(diameters, field_below, 'o')
    plt.grid(True)
    plt.savefig('./img/eFieldBelow.png', dpi = 300)
    plt.show()