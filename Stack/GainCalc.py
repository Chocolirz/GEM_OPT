import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
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



    # ITO Drift, sigma = 0.022 cm, gain = 32561 / 1e5



PlotThreeGEMsAligned = True
if PlotThreeGEMsAligned:
    # Your given drift results
    drift_results_aligned = np.array([
        4911, 
        2859, 2859, 2859, 2859, 2859, 2859, 
        902.7, 902.7, 902.7, 902.7, 902.7, 902.7, 
        511, 511, 511, 511, 511, 511, 
        82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 82.3, 
        27.2, 27.2, 27.2, 27.2, 27.2, 27.2,
        4.7, 4.7, 4.7, 4.7, 4.7, 4.7, 
        3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 3.3, 
        0.7, 0.7, 0.7, 0.7, 0.7, 0.7
    ]) / 1e5

    # Lattice generation
    a = 0.0280
    R_max = 5 * a
    points = []
    max_rows = int(R_max / (math.sqrt(3)/2 * a)) + 2
    max_cols = int(R_max / a) + 2

    for row in range(-max_rows, max_rows + 1):
        y = row * (math.sqrt(3)/2) * a
        for col in range(-max_cols, max_cols + 1):
            x = col * a + (a/2 if row % 2 else 0)
            r = math.hypot(x, y)
            if r <= R_max + 1e-12:
                points.append((x, y, r))

    points.sort(key=lambda p: round(p[2], 12))
    positions = [(x, y) for x, y, _ in points]

    # First to second layer probability map
    prob_map = {pos: p for pos, p in zip(positions, drift_results_aligned)}

    # Convolution to get 1st → 3rd layer probability
    third_layer_map = {}
    for p1, prob1 in zip(positions, drift_results_aligned):
        for p2, prob2 in zip(positions, drift_results_aligned):
            x_new = p1[0] + p2[0]
            y_new = p1[1] + p2[1]
            third_layer_map[(x_new, y_new)] = third_layer_map.get((x_new, y_new), 0) + prob1 * prob2



    ### Move on the the ITO



    # Parameters
    FWHM = 0.044  # cm
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))  # Gaussian sigma from FWHM

    # Grid for continuous map
    grid_size = 500  # number of points along each axis
    x_range = (-0.22, 0.22)  # adjust as needed
    y_range = (-0.22, 0.22)
    x = np.linspace(*x_range, grid_size)
    y = np.linspace(*y_range, grid_size)
    X, Y = np.meshgrid(x, y)
    Z_new = np.zeros_like(X)

    # Superpose Gaussian contributions from each lattice point
    for (x0, y0), P in third_layer_map.items():
        Z_new += 84 * 84.32 * 84.32 * P * np.exp(-((X-x0)**2 + (Y-y0)**2)/(2*sigma**2))

    # Plot continuous color map
    plt.figure(figsize=(6,6))
    im = plt.imshow(Z_new, extent=[x_range[0], x_range[1], y_range[0], y_range[1]],
                    origin='lower', cmap='flag') # norm=LogNorm() if needed
    plt.colorbar(im, label='Number of Electrons')
    #plt.xlim((0.05, 0.075))
    #plt.ylim((0.05, 0.075))
    plt.xlabel("x [cm]")
    plt.ylabel("y [cm]")
    plt.show()

    # Compute std
    # Normalize Z_new to make it a proper probability density
    Z_norm = Z_new / np.sum(Z_new)

    # Compute expected values (means)
    x_mean = np.sum(X * Z_norm)
    y_mean = np.sum(Y * Z_norm)

    # Compute variances
    x_var = np.sum(((X - x_mean)**2) * Z_norm)
    y_var = np.sum(((Y - y_mean)**2) * Z_norm)

    # Standard deviations
    x_std = np.sqrt(x_var)
    y_std = np.sqrt(y_var)

    print(f"Standard deviation in x: {x_std:.6f} cm")
    print(f"Standard deviation in y: {y_std:.6f} cm")



# --------------------------------------------------------------------------------------------------------------------------------------------




PlotThreeGEMsMisaligned = False
if PlotThreeGEMsMisaligned:
    # Your given drift results
    drift_results_misaligned = np.array([
        4512.5, 3358, 1872.75, 4512.5, 1376, 3358, 1872.75, # r = 1/2, s3/2, s7/2, 1/2, 1.5, s3/2, s7/2
        769, 1872.75, 155, 1872.75, 155, 769, # r = s13/2, s7/2, s19/2, s7/2, s19/2, s13/2
        769, 122.25, 1376, 122.5, 769, 122.25, # s13/2, s21/2, 1.5, 2.5, s13/2, s21/2
        89.5, 57.25, 122.25, 18, 155, 16.75, 155, 16.75, 122.25, 18, 89.5, 57.25, 
        # 1.5s3, s31/2, s21/2, s37/2, s19/2, s39/2, s19/2, s39/2, s21/2, s37/2, 1.5s3, s31/2
        57.25, 7.5, 122.5, 5, 57.25, 7.5, # s31/2, s43/2, 2.5, 3.5, s31/2, s43/2
        5, 18, 0, 18, 0, 5, # 3.5, s37/2, _, s37/2, _, 3.5
        5, 0.75, 7.5, 0, 16.75, 0, 16.75, 0, 7.5, 0, 5, 0.75, 
        # 3.5, s57/2, s43/2, _, s39/2, _,           s39/2, _, s43/2, _, 3.5, s57/2
        0.75, 0, 5, 0, 0.75, 0 # s57/2, _, 3.5, _, s57/2, _
    ]) / 1e5

    # Lattice generation
    a = 0.0280
    R_max = 4 * a
    points = []
    max_rows = int(R_max / (math.sqrt(3)/2 * a)) + 2
    max_cols = int(R_max / a) + 2

    for row in range(-max_rows, max_rows + 1):
        y = row * (math.sqrt(3)/2) * a
        for col in range(-max_cols, max_cols + 1):
            x = col * a + (a/2 if row % 2 else 0)
            r = math.hypot(x, y)
            if r <= R_max + 1e-12:
                points.append((x, y, r))

    points.sort(key=lambda p: round(p[2], 12))
    positions = [(x + 0.5*a, y) for x, y, _ in points]
    
    # First to second layer probability map
    prob_map = {pos: p for pos, p in zip(positions, drift_results_misaligned)}

    # Convolution to get 1st → 3rd layer probability
    third_layer_map = {}
    for p1, prob1 in zip(positions, drift_results_misaligned):
        for p2, prob2 in zip(positions, drift_results_misaligned):
            x_new = p1[0] + p2[0]
            y_new = p1[1] + p2[1]
            third_layer_map[(x_new, y_new)] = third_layer_map.get((x_new, y_new), 0) + prob1 * prob2






    ### Move on the the ITO



    # Parameters
    FWHM = 0.044  # cm
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))  # Gaussian sigma from FWHM

    # Grid for continuous map
    grid_size = 500  # number of points along each axis
    x_range = (-0.22, 0.22)  # adjust as needed
    y_range = (-0.22, 0.22)
    x = np.linspace(*x_range, grid_size)
    y = np.linspace(*y_range, grid_size)
    X, Y = np.meshgrid(x, y)
    Z_new = np.zeros_like(X)

    # Superpose Gaussian contributions from each lattice point
    for (x0, y0), P in third_layer_map.items():
        Z_new += 84 * 84.32 * 84.32 * P * np.exp(-((X-x0)**2 + (Y-y0)**2)/(2*sigma**2))

    # Plot continuous color map
    plt.figure(figsize=(6,6))
    im = plt.imshow(Z_new, extent=[x_range[0], x_range[1], y_range[0], y_range[1]],
                    origin='lower', cmap='viridis') # norm=LogNorm() if needed
    plt.colorbar(im, label='Number of Electrons')
    plt.xlabel("x [cm]")
    plt.ylabel("y [cm]")
    plt.show()


    # Compute std
    # Normalize Z_new to make it a proper probability density
    Z_norm = Z_new / np.sum(Z_new)

    # Compute expected values (means)
    x_mean = np.sum(X * Z_norm)
    y_mean = np.sum(Y * Z_norm)

    # Compute variances
    x_var = np.sum(((X - x_mean)**2) * Z_norm)
    y_var = np.sum(((Y - y_mean)**2) * Z_norm)

    # Standard deviations
    x_std = np.sqrt(x_var)
    y_std = np.sqrt(y_var)

    print(f"Standard deviation in x: {x_std:.6f} cm")
    print(f"Standard deviation in y: {y_std:.6f} cm")