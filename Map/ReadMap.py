import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import numpy as np
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

# define some constants
rlim = 0.0210 / 2.;  # cm, radius of the GEM hole
upperGap = 0.2;  # cm
lowerGap = 0.2;  # cm
copperThick = 0.0002;  # cm, thickness of the copper plates
gap = 0.057;  # cm, GEM is 570 um thick
pitch = 0.03;  # distance between GEM holes
xmin = -10. * pitch; 
xmax = 10. * pitch;  # cm
ymin = -10. * pitch;
ymax = 10. * pitch;  # cm
zmin = 0.;
zmax = 3. * gap + 3. * lowerGap + upperGap + (3. * 2. * copperThick);  # cm, three GEMs
tsg = 2. * gap + 2. * lowerGap + 4. * copperThick;  # Height of the top of the second GEM.

# Lattice generation
a = 0.0280
R_max = 9 * a
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
positions_inv = [(x + 0.5*pitch, y) for x, y, _ in points]


# x, y, z, theta, phi of the scintillation photons 
# (that can reach the top of next GEM)
# were stored in a file called scint_positions.txt
# read it, 
data = np.loadtxt('./scint_positions.txt', skiprows=1)
x = data[:, 0]
y = data[:, 1]
z = data[:, 2] - tsg # make it the relative position, counting from the tsg
theta = data[:, 3]
phi = data[:, 4]
r = np.sqrt(x**2 + y**2)

# and check the starting positions
CheckStartingPositions = False
if CheckStartingPositions:
    print(z[:100])
    plt.figure(figsize=(8,8))
    plt.xlabel('r [cm]')
    plt.ylabel('z [cm]')
    plt.plot(x[:1000], z[:1000], 'o')
    plt.grid(True)
    plt.show()

# Define a function that judges if the photon is going to pass a circle centred at (xc,yc,zc) with radius rlim:
def through_circle(xc, yc, zc, rlim, x, y, z, theta, phi):
    dxy = (z - zc) * np.tan(np.pi - theta) # angle in radians to match np.sin()
    dx = dxy * np.cos(phi)
    dy = dxy * np.sin(phi)
    dr = np.sqrt((x + dx - xc)**2 + (y + dy - yc)**2)
    if dr < rlim:
        return True
    else:
        return False
    
# Define a function that judges if the photon is going to pass a hole centred at (xc, yc, zc-(gap/2)) with radius rlim:
def through_hole(xc, yc, zc, rlim, x, y, z, theta, phi):
    thick = 0.0574
    if theta > 2.79:
        if through_circle(xc, yc, zc, rlim, x, y, z, theta, phi) == True:
            if through_circle(xc, yc, zc - thick, rlim, x, y, z, theta, phi) == True:
                return True
            else:
                return False
        else:
            return False
    else:
        return False

# Define a function that judges if the photon is going to pass the bottom layer of GEM
def pass_GEM(positions, zc, rlim, x, y, z, theta, phi):
    x_list = [x for x, _ in positions]
    y_list = [y for _, y in positions]
    judge = 0
    for i in range(len(x_list)):
        if through_hole(x_list[i], y_list[i], zc, rlim, x, y, z, theta, phi) == True:
            judge += 1
        else:
            pass
    if judge == 1:
        return True
    else:
        return False
    
# Define a function that counts the number of photons passing the bottom layer of GEM
# Here x, y, z, theta, phi are lists
def pass_GEM_counter(positions, zc, rlim, x_pos, y_pos, z_pos, theta_pos, phi_pos):
    count = 0
    indices = []
    for i in tqdm(range(len(x_pos))):
        if pass_GEM(positions, zc, rlim, x_pos[i], y_pos[i], z_pos[i], theta_pos[i], phi_pos[i]) == True:
            count += 1
            indices.append(i)
        else:
            pass
    return count, indices




######------------------------------------------------------------------------------
'''
PlotBFG = False
if PlotBFG:
    x_bfg = x + (z-lowerGap)*np.tan(np.pi-theta)*np.cos(phi)
    y_bfg = y + (z-lowerGap)*np.tan(np.pi-theta)*np.sin(phi)
    r_bfg = np.sqrt(x_bfg**2 + y_bfg**2)
    weights = 1/r_bfg

    plt.figure(figsize=(8,8))
    plt.axis('equal')
    plt.xlabel('$x$ [cm]')
    plt.ylabel('$y$ [cm]')
    plt.scatter(x_bfg, y_bfg)
    #plt.hist2d(x + z*np.tan(np.pi-theta)*np.cos(phi), y + z*np.tan(np.pi-theta)*np.sin(phi), bins=[1000,1000], range=[[-1,1],[-1,1]], cmap='plasma')
    plt.xlim(-0.05, 0.05)
    plt.ylim(-0.05, 0.05)
    plt.grid(True)
    plt.savefig('./img/bfg_photon_map.png', dpi = 300)
    plt.show()

    plt.figure(figsize=(8,8))
    plt.xlabel(r'$R$ [cm]')
    plt.ylabel('Entries')
    plt.hist(r_bfg, bins = 200, range=[0,0.05], weights = weights, histtype='step')
    plt.xlim((0,0.05))
    plt.yscale('log')
    plt.grid(True)
    plt.savefig('./img/r_bfg.png', dpi = 300)
    plt.show()
    

PlotTSG = True
if PlotTSG:
    x_tsg = x + z*np.tan(np.pi-theta)*np.cos(phi)
    y_tsg = y + z*np.tan(np.pi-theta)*np.sin(phi)
    r_tsg = np.sqrt(x_tsg**2 + y_tsg**2)
    weights = 1/r_tsg

    plt.figure(figsize=(8,8))
    plt.axis('equal')
    plt.xlabel('$x$ [cm]')
    plt.ylabel('$y$ [cm]')
    plt.scatter(x_tsg, y_tsg)
    #plt.hist2d(x + z*np.tan(np.pi-theta)*np.cos(phi), y + z*np.tan(np.pi-theta)*np.sin(phi), bins=[1000,1000], range=[[-1,1],[-1,1]], cmap='plasma')
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.grid(True)
    plt.savefig('./img/tsg_photon_map.png', dpi = 300)
    plt.show()

    plt.figure(figsize=(8,8))
    plt.xlabel(r'$R$ [cm]')
    plt.ylabel('Entries')
    plt.hist(r_tsg, bins = 200, range=[0,1.5], weights = weights, histtype='step')
    plt.yscale('log')
    plt.xlim((0,1.5))
    plt.grid(True)
    plt.savefig('./img/r_tsg.png', dpi = 300)
    plt.show()



### Execute the programme
ttg_count, ttg_indices = pass_GEM_counter(positions_inv, 0, rlim, x, y, z, theta, phi)
print(ttg_count, '/', len(x), ' photons passes the next GEM')

### Make some plots

# Plot the position those passing the second GEM at the top of the third GEM
x_ttg = []
y_ttg = []
theta_ttg = []

for j in ttg_indices:
    x_ttg.append(x[j] + (z[j] + gap + (2*copperThick) + lowerGap) * np.tan(np.pi - theta[j]) * np.cos(phi[j]))
    y_ttg.append(y[j] + (z[j] + gap + (2*copperThick) + lowerGap) * np.tan(np.pi - theta[j]) * np.sin(phi[j]))
    theta_ttg.append(theta[j])

PlotTTG = True
if PlotTTG:
    plt.figure(figsize=(8,8))
    plt.axis('equal')
    plt.xlabel('$x$ [cm]')
    plt.ylabel('$y$ [cm]')
    plt.scatter(x_ttg, y_ttg, label='TTG')
    plt.legend()
    plt.grid(True)
    plt.savefig('./img/ttg_photon_map_inv.png', dpi = 300)
    plt.show()

    plt.figure(figsize=(8,8))
    plt.xlabel(r'$\theta$ [rad]')
    plt.ylabel('Entries')
    plt.hist(theta_ttg, bins=45, histtype='step')
    plt.grid(True)
    plt.savefig('./img/ttg_theta_distribution_inv.png', dpi = 300)
    plt.show()


x_ptg = [x[i] for i in ttg_indices] # photons passed the second GEM
y_ptg = [y[i] for i in ttg_indices]
z_ptg = [z[i] for i in ttg_indices]
theta_ptg = [theta[i] for i in ttg_indices]
phi_ptg = [phi[i] for i in ttg_indices]


ito_count, ito_indices = pass_GEM_counter(positions, 0-gap-lowerGap-(2*copperThick), rlim, x_ptg, y_ptg, z_ptg, theta_ptg, phi_ptg)
print(ito_count, '/', len(x_ptg), 'photons further pass the third GEM and reaches ITO.')

# Plot the position those passing the second GEM at the top of the third GEM
x_ito = []
y_ito = []
theta_ito = []

# Note that ito_indices are indices in _ptg, not the original x, therefore we need to call the ttg_indices[j]-th component of x
for j in ito_indices:
    x_ito.append(x[ttg_indices[j]] + (z[ttg_indices[j]] + 2*gap + 2*lowerGap) * np.tan(np.pi - theta[ttg_indices[j]]) * np.cos(phi[ttg_indices[j]]))
    y_ito.append(y[ttg_indices[j]] + (z[ttg_indices[j]] + 2*gap + 2*lowerGap) * np.tan(np.pi - theta[ttg_indices[j]]) * np.sin(phi[ttg_indices[j]]))
    theta_ito.append(theta[ttg_indices[j]])

PlotITO = True
if PlotITO:
    plt.figure(figsize=(8,8))
    plt.axis('equal')
    plt.xlabel('$x$ [cm]')
    plt.ylabel('$y$ [cm]')
    plt.scatter(x_ito, y_ito, label='ITO')
    plt.legend()
    plt.grid(True)
    plt.savefig('./img/ITO_photon_map_inv.png', dpi = 300)
    plt.show()

    plt.figure(figsize=(8,8))
    plt.xlabel(r'$\theta$ [rad]')
    plt.ylabel('Entries')
    plt.hist(theta_ito, bins=45, histtype='step')
    plt.grid(True)
    plt.savefig('./img/ITO_theta_distribution_inv.png', dpi = 300)
    plt.show()

    plt.figure(figsize=(8,8))
    plt.xlabel(r'$\theta$ [rad]')
    plt.ylabel('Entries')
    plt.hist(theta_ito, bins=45, histtype='step', label='ITO')
    plt.hist(theta_ttg, bins=45, histtype='step', label='TTG')
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.savefig('./img/theta_distribution_combined_inv.png', dpi = 300)
    plt.show()
'''