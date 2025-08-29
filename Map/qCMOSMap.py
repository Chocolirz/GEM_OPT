import numpy as np
#import matplotlib.pyplot as plt
from tqdm import tqdm
import ReadMap as map

# Import data for every photon
print("Loading photon distribution on the plane of focus...")
x = np.loadtxt('x_misaligned.txt')
print('1/4')
y = np.loadtxt('y_misaligned.txt')
print("2/4")
theta = np.loadtxt('theta_misaligned.txt')
print('3/4')
phi = np.loadtxt('phi_misaligned.txt')
print("Data successfully loaded. ")

# Initialise position:
x_0 = np.linspace(0,4,9) # cm
y_0 = np.linspace(0,4,9) # cm
z_new = -2*(map.gap + 2*map.copperThick)-map.lowerGap

# with MOD 170 mm and minimum angle of view 30.3 degrees (gives 9.8 cm range of view), 
# the lens can accept every photon from the GEM plate, 
# as long as it goes through the aperture. 

# Determine if the photon is going to pass through the aperture
# we use the through_circle function again
# Lens surface is 7.71 cm below BTG
# radius regarded as 27.5 mm

#through_aperture_indices = []
ratio = []
through_aperture_number = 0
for x_shift in x_0:
    x_new = x + x_shift
    for y_shift in y_0:
        y_new = y + y_shift
        for i in tqdm(range(len(x))):
            if map.through_circle(0, 0, -2*(map.gap + 2*map.copperThick)-map.lowerGap-7.71, 2.75, x_new[i], y_new[i], z_new, theta[i], phi[i]):
                #through_aperture_indices.append[i]
                through_aperture_number += 1
        print(y_shift, through_aperture_number/len(x))
        ratio.append([x_shift, y_shift, through_aperture_number/len(x)])
    print(x_shift)

print(ratio)
np.savetxt('ratio_mis.txt', ratio)
'''
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ratio = np.array(ratio)  # make it a numpy array for convenience
ax.scatter(ratio[:,0], ratio[:,1], ratio[:,2], marker='o')

ax.set_xlabel(r'$x$ [cm]')
ax.set_ylabel(r'$y$ [cm]')
ax.set_zlabel(r'Collection ratio')

plt.show()
'''