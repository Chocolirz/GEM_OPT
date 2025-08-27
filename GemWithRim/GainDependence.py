import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
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

data = np.loadtxt('./rGain.txt', skiprows=1)
divide = 1
r = data[:, 0]  # radius in cm
gain = data[:, 1]  # gain
anode_ratio = data[:, 2]  # anode ratio
cu_ratio = data[:, 3]  # production ratio
glass_ratio = data[:, 4]  # rim in cm
att_hole_ratio = data[:, 5]
att_drift_ratio = data[:, 6]
production = data[:, 7]
rim = data[:, 8]
att_ratio = []
att_number = []
hit_ratio = []
hit_number = []
total_number = []
para_gain = 21000
para_att = 218
for i in range(len(r)):
    att_ratio.append(att_hole_ratio[i] + att_drift_ratio[i])
    att_number.append(att_ratio[-1]*(gain[i]/anode_ratio[i]))
    hit_ratio.append(cu_ratio[i] + glass_ratio[i])
    hit_number.append(hit_ratio[-1]*(gain[i]/anode_ratio[i]))
    total_number.append(gain[i]/(0.01 *anode_ratio[i]))


PlotRG = False
if PlotRG:
    plt.figure(figsize=(8, 8))
    plt.title('Gain as a function of radius')
    plt.xlabel('Diameter [um]')
    plt.ylabel('Gain')
    for i in range(len(r)):
        if rim[i] == 0:
            plt.plot(r[i], gain[i], 'o', color='blue', label=f'{r[i]} um')
    plt.grid(True)
    plt.savefig('./img/Gain_vs_Radius.png', dpi=300)
    plt.show()

PlotRatio = False
if PlotRatio:
    plt.figure(figsize=(8,8))
    plt.title("Ratio as a function of radius")
    plt.xlabel('Diameter [um]')
    plt.ylabel('Ratio')
    plt.plot(r, anode_ratio, 'o', color = 'blue', label = "Anode")
    plt.plot(r, hit_ratio, 'o', color = 'red', label = 'Hit')
    plt.plot(r, att_ratio, 'o', color = 'green', label = "Attachment")
    plt.grid(True)
    plt.legend()
    plt.savefig('./img/Ratio_vs_Radius.png', dpi = 300)
    plt.show()

PlotNumber = False
if PlotNumber:
    plt.figure(figsize=(8,8))
    plt.xlabel('Diameter [um]')
    plt.ylabel('# Electrons')
    plt.yscale('log')
    plt.plot(r[divide:], gain[divide:], 'o', color = 'blue', label = "Anode")
    plt.plot(r[divide:], hit_number[divide:], 'o', color = 'red', label = 'Hit')
    plt.plot(r[divide:], att_number[divide:], 'o', color = 'green', label = "Attachment")
    plt.plot(r[divide:], total_number[divide:], 'o', color = 'purple', label = 'Total')
    plt.plot(r[:divide], gain[:divide], 's', color = 'blue')
    plt.plot(r[:divide], hit_number[:divide], 's', color = 'red')
    plt.plot(r[:divide], att_number[:divide], 's', color = 'green')
    plt.plot(r[:divide], total_number[:divide], 's', color = 'purple')
    plt.plot([280], [para_gain], 'o', color = 'blue')
    plt.plot([280], [para_att], 'o', color = 'green')
    plt.grid(True)
    plt.legend()
    plt.savefig('./img/Number_vs_Radius.png', dpi = 300)
    plt.show()

PlotdVGain = True
if PlotdVGain:
    datadVGain = np.loadtxt("./dVGain.txt", skiprows=1)
    radius = datadVGain[:, 0]
    voltages = datadVGain[:, 1]
    Gains = datadVGain[:, 2]
    field_above = datadVGain[:, 3]
    field_below = datadVGain[:, 4]

    def expFit(x,A,B,C):
        return A*np.exp(B*x + C)
    
    popt170, pcov170 = curve_fit(expFit, voltages[:7], Gains[:7], p0=(1, 0.05, -30))
    popt210, pcov210 = curve_fit(expFit, voltages[7:], Gains[7:], p0=(1, 0.05, -30))
    v_list = np.linspace(540, 600, 600)
    
    plt.figure(figsize=(8,8))
    plt.xlabel('Voltages [V]')
    plt.ylabel('Gain')
    #plt.yscale('log')
    plt.plot(v_list, expFit(v_list, *popt170), label = f'$y={round(popt170[0],2)}exp({round(popt170[1],4)} x +{round(popt170[2],1)})$', color = '#4477aa')
    plt.plot(v_list, expFit(v_list, *popt210), label = f'$y={round(popt210[0],2)}exp({round(popt210[1],4)} x +{round(popt210[2],1)})$', color = '#ee6677')
    plt.plot(voltages[:7], Gains[:7], 'o', label = '170um', color = '#4477aa')
    plt.plot(voltages[7:], Gains[7:], 'o', label = '210um', color = '#ee6677')
    plt.grid(True)
    plt.legend()
    plt.savefig('./img/Voltages_vs_Gain.png', dpi = 300)
    plt.show()

    plt.figure(figsize=(8,8))
    plt.xlabel('Voltages [V]')
    plt.ylabel('Field')
    plt.plot(voltages[:7], field_above[:7], label = 'Drift, 170um')
    plt.plot(voltages[7:], field_above[7:] - 25, label = 'Drift - 25V/cm, 210um')
    plt.plot(voltages[:7], field_below[:7]-200, label = 'Induction - 190V/cm, 170um')
    plt.plot(voltages[7:], field_below[7:]-225, label = 'Induction - 210V/cm, 210um')
    plt.yscale('log')
    plt.grid(True, which='major')
    plt.grid(True, which='minor')
    plt.legend()
    plt.savefig('./img/Field_vs_Voltage_Drift_and_Induction.png', dpi = 300)
    plt.show()