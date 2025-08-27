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

plt.rcParams.update(params)

log_tick_format = ticker.FuncFormatter(lambda y,pos: f'{y:.{int(np.maximum(-np.log10(y),0)):1d}f}')

PlotRBot = True
if PlotRBot:
    data = np.loadtxt('rBotRelation.txt', skiprows=2)
    rTop = data[:, 0]
    rBot = data[:, 1]
    gain = data[:, 2]

    plt.figure(figsize=(8, 8))
    for i in range(len(rTop)):
        if rTop[i] == 150:
            color = '#4477aa'  # blue
        elif rTop[i] == 160:
            color = '#ee6677'  # red/pink
        elif rTop[i] == 170:
            color = '#228833'  # green
        elif rTop[i] == 180:    
            color = '#aa3377'  # purple
        elif rTop[i] == 190:
            color = '#66ccee'  # cyan
        elif rTop[i] == 200:
            color = '#ccbb44'  # yellow
        elif rTop[i] == 210:
            color = '#bbbbbb'  # grey
        plt.plot(rBot[i], gain[i], 'o', color = color)
    #plt.title(f'Gain dependence for {rTop[0]} um r_top')
    plt.xlabel('$r_{bottom}$ [um]')
    plt.ylabel('Gain')
    plt.yscale('log')
    plt.grid(True)
    patch_150 = mpatches.Patch(color='#4477aa', label='150 um')
    patch_160 = mpatches.Patch(color='#ee6677', label='160 um')
    patch_170 = mpatches.Patch(color='#228833', label='170 um')
    patch_180 = mpatches.Patch(color='#aa3377', label='180 um')
    patch_190 = mpatches.Patch(color='#66ccee', label='190 um')
    patch_200 = mpatches.Patch(color='#ccbb44', label='200 um')
    patch_210 = mpatches.Patch(color='#bbbbbb', label='210 um')
    plt.legend(handles=[patch_150, patch_160, patch_170], loc = 'lower right')
    plt.savefig(f'./img/rBotRelation.png', dpi=300)
    plt.show()