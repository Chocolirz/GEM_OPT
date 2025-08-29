import ReadMap as map
import ThirdLayerNumber as tln
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

x_bfg = map.x + (map.z-map.lowerGap)*np.tan(np.pi-map.theta)*np.cos(map.phi)
y_bfg = map.y + (map.z-map.lowerGap)*np.tan(np.pi-map.theta)*np.sin(map.phi)
theta_bfg = map.theta
phi_bfg = map.phi

# Change this to map.positions_inv
ttg_count, ttg_indices = map.pass_GEM_counter(map.positions_inv, 0, map.rlim, map.x, map.y, map.z, map.theta, map.phi)
print(ttg_count, '/', len(map.x), ' photons passes the next GEM')

# Plot the position those passing the second GEM at the top of the third GEM
x_bsg = []
y_bsg = []
theta_bsg = []
phi_bsg = []

for j in ttg_indices:
    x_bsg.append(map.x[j] + (map.z[j] + map.gap + (2*map.copperThick)) * np.tan(np.pi - map.theta[j]) * np.cos(map.phi[j]))
    y_bsg.append(map.y[j] + (map.z[j] + map.gap + (2*map.copperThick)) * np.tan(np.pi - map.theta[j]) * np.sin(map.phi[j]))
    theta_bsg.append(map.theta[j])
    phi_bsg.append(map.phi[j])

x_ptg = [map.x[i] for i in ttg_indices] # photons passed the second GEM
y_ptg = [map.y[i] for i in ttg_indices]
z_ptg = [map.z[i] for i in ttg_indices]
theta_ptg = [map.theta[i] for i in ttg_indices]
phi_ptg = [map.phi[i] for i in ttg_indices]


ito_count, ito_indices = map.pass_GEM_counter(map.positions, 0-map.gap-map.lowerGap-(2*map.copperThick), map.rlim, x_ptg, y_ptg, z_ptg, theta_ptg, phi_ptg)
print(ito_count, '/', len(x_ptg), 'photons further pass the third GEM and reaches ITO.')

# Plot the position those passing the second GEM at the top of the third GEM
x_btg = []
y_btg = []
theta_btg = []
phi_btg = []

# Note that ito_indices are indices in _ptg, not the original x, therefore we need to call the ttg_indices[j]-th component of x
for j in ito_indices:
    x_btg.append(map.x[ttg_indices[j]] + (map.z[ttg_indices[j]] + 2*map.gap + map.lowerGap) * np.tan(np.pi - map.theta[ttg_indices[j]]) * np.cos(map.phi[ttg_indices[j]]))
    y_btg.append(map.y[ttg_indices[j]] + (map.z[ttg_indices[j]] + 2*map.gap + map.lowerGap) * np.tan(np.pi - map.theta[ttg_indices[j]]) * np.sin(map.phi[ttg_indices[j]]))
    theta_btg.append(map.theta[ttg_indices[j]])
    phi_btg.append(map.phi[ttg_indices[j]])

total_distribution_x = []
total_distribution_y = []
total_distribution_theta = []
total_distribution_phi = []

# first layer - use btg
for _ in range(1):
    total_distribution_x.extend(x_btg)
    total_distribution_y.extend(y_btg)
    total_distribution_theta.extend(theta_btg)
    total_distribution_phi.extend(phi_btg)
print('Generated first layer. ')
# second layer - use bsg
for _ in range(15):
    total_distribution_x.extend(x_bsg)
    total_distribution_y.extend(y_bsg)
    total_distribution_theta.extend(theta_bsg)
    total_distribution_phi.extend(phi_bsg)
print('Generating second layer...')
for _ in range(9):
    total_distribution_x.extend([x + map.pitch for x in x_bsg])
    total_distribution_x.extend([x - map.pitch for x in x_bsg])
    total_distribution_x.extend([x + 0.5*map.pitch for x in x_bsg])
    total_distribution_x.extend([x - 0.5*map.pitch for x in x_bsg])
    total_distribution_x.extend([x + 0.5*map.pitch for x in x_bsg])
    total_distribution_x.extend([x - 0.5*map.pitch for x in x_bsg])
    total_distribution_y.extend([y for y in y_bsg])
    total_distribution_y.extend([y for y in y_bsg])
    total_distribution_y.extend([y + 0.5*np.sqrt(3)*map.pitch for y in y_bsg])
    total_distribution_y.extend([y + 0.5*np.sqrt(3)*map.pitch for y in y_bsg])
    total_distribution_y.extend([y - 0.5*np.sqrt(3)*map.pitch for y in y_bsg])
    total_distribution_y.extend([y - 0.5*np.sqrt(3)*map.pitch for y in y_bsg])
    for _ in range(6):
        total_distribution_theta.extend(theta_bsg)
        total_distribution_phi.extend(phi_bsg)
for _ in range(2):
    total_distribution_x.extend([x for x in x_bsg])
    total_distribution_x.extend([x for x in x_bsg])
    total_distribution_x.extend([x + 1.5*map.pitch for x in x_bsg])
    total_distribution_x.extend([x - 1.5*map.pitch for x in x_bsg])
    total_distribution_x.extend([x + 1.5*map.pitch for x in x_bsg])
    total_distribution_x.extend([x - 1.5*map.pitch for x in x_bsg])
    total_distribution_y.extend([y + np.sqrt(3)*map.pitch for y in y_bsg])
    total_distribution_y.extend([y - np.sqrt(3)*map.pitch for y in y_bsg])
    total_distribution_y.extend([y + 0.5*np.sqrt(3)*map.pitch for y in y_bsg])
    total_distribution_y.extend([y + 0.5*np.sqrt(3)*map.pitch for y in y_bsg])
    total_distribution_y.extend([y - 0.5*np.sqrt(3)*map.pitch for y in y_bsg])
    total_distribution_y.extend([y - 0.5*np.sqrt(3)*map.pitch for y in y_bsg])
    for _ in range(6):
        total_distribution_theta.extend(theta_bsg)
        total_distribution_phi.extend(phi_bsg)
print('Generated second layer.')
# third layer - use bfg
for _ in range(tln.photon_values[0]):
    total_distribution_x.extend(x_bfg)
    total_distribution_y.extend(y_bfg)
    total_distribution_theta.extend(theta_bfg)
    total_distribution_phi.extend(phi_bfg)
print('Generating third layer...')
for _ in tqdm(range(tln.photon_values[1])):
    total_distribution_x.extend([x + map.pitch for x in x_bfg])
    total_distribution_x.extend([x - map.pitch for x in x_bfg])
    total_distribution_x.extend([x + 0.5*map.pitch for x in x_bfg])
    total_distribution_x.extend([x - 0.5*map.pitch for x in x_bfg])
    total_distribution_x.extend([x + 0.5*map.pitch for x in x_bfg])
    total_distribution_x.extend([x - 0.5*map.pitch for x in x_bfg])
    total_distribution_y.extend([y for y in y_bfg])
    total_distribution_y.extend([y for y in y_bfg])
    total_distribution_y.extend([y + 0.5*np.sqrt(3)*map.pitch for y in y_bfg])
    total_distribution_y.extend([y + 0.5*np.sqrt(3)*map.pitch for y in y_bfg])
    total_distribution_y.extend([y - 0.5*np.sqrt(3)*map.pitch for y in y_bfg])
    total_distribution_y.extend([y - 0.5*np.sqrt(3)*map.pitch for y in y_bfg])
    for _ in range(6):
        total_distribution_theta.extend(theta_bfg)
        total_distribution_phi.extend(phi_bfg)
print('Generating third layer...')
for _ in tqdm(range(tln.photon_values[7])):
    total_distribution_x.extend([x for x in x_bfg])
    total_distribution_x.extend([x for x in x_bfg])
    total_distribution_x.extend([x + 1.5*map.pitch for x in x_bfg])
    total_distribution_x.extend([x - 1.5*map.pitch for x in x_bfg])
    total_distribution_x.extend([x + 1.5*map.pitch for x in x_bfg])
    total_distribution_x.extend([x - 1.5*map.pitch for x in x_bfg])
    total_distribution_y.extend([y + np.sqrt(3)*map.pitch for y in y_bfg])
    total_distribution_y.extend([y - np.sqrt(3)*map.pitch for y in y_bfg])
    total_distribution_y.extend([y + 0.5*np.sqrt(3)*map.pitch for y in y_bfg])
    total_distribution_y.extend([y + 0.5*np.sqrt(3)*map.pitch for y in y_bfg])
    total_distribution_y.extend([y - 0.5*np.sqrt(3)*map.pitch for y in y_bfg])
    total_distribution_y.extend([y - 0.5*np.sqrt(3)*map.pitch for y in y_bfg])
    for _ in range(6):
        total_distribution_theta.extend(theta_bfg)
        total_distribution_phi.extend(phi_bfg)
print('Generating third layer...')
for _ in range(tln.photon_values[13]):
    total_distribution_x.extend([x + 2*map.pitch for x in x_bfg])
    total_distribution_x.extend([x - 2*map.pitch for x in x_bfg])
    total_distribution_x.extend([x + map.pitch for x in x_bfg])
    total_distribution_x.extend([x - map.pitch for x in x_bfg])
    total_distribution_x.extend([x + map.pitch for x in x_bfg])
    total_distribution_x.extend([x - map.pitch for x in x_bfg])
    total_distribution_y.extend([y for y in y_bfg])
    total_distribution_y.extend([y for y in y_bfg])
    total_distribution_y.extend([y + np.sqrt(3)*map.pitch for y in y_bfg])
    total_distribution_y.extend([y + np.sqrt(3)*map.pitch for y in y_bfg])
    total_distribution_y.extend([y - np.sqrt(3)*map.pitch for y in y_bfg])
    total_distribution_y.extend([y - np.sqrt(3)*map.pitch for y in y_bfg])
    for _ in range(6):
        total_distribution_theta.extend(theta_bfg)
        total_distribution_phi.extend(phi_bfg)
print('Generated third layer')

print('Generating plots...')
h, xedges, yedges = np.histogram2d(
    total_distribution_y, total_distribution_x,
    bins=[100, 100],
    range=[[-0.1, 0.1], [-0.1, 0.1]]
)

plt.figure(figsize=(10,10))
plt.xlabel('$x$ [cm]')
plt.ylabel('$y$ [cm]')
im = plt.imshow(h/100, origin='lower', aspect='equal', extent=[-0.1,0.1,-0.1,0.1])
plt.colorbar(im, label='Number of photons')
plt.savefig('./img/event_photon_map_inv.png', dpi = 300)
plt.show()

print("Writing file:")
np.savetxt("x_misaligned.txt", total_distribution_x)
print("1/4")
np.savetxt("y_misaligned.txt", total_distribution_y)
print("2/4")
np.savetxt("theta_misaligned.txt", total_distribution_theta)
print("3/4")
np.savetxt("phi_misaligned.txt", total_distribution_phi)
print("Complete!")