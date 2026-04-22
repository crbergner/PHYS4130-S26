# import libraries
import numpy as np
from matplotlib import pyplot as plt
import random as ra
import matplotlib.animation as an

# functions for main code ############################################################
def move(px, py):
    direction = ra.randint(0, 3)
    if direction == 0: py += 1
    elif direction == 1: py -= 1
    elif direction == 2: px -= 1
    elif direction == 3: px += 1
    return px, py

def sticking(px, py, grid, stickiness):
    if np.any(grid[px-1:px+2, py-1:py+2] > 0):
        if ra.random() <= stickiness:
            return True
    return False

# pick angle where particles are not spawning
def balanced_spawn_angle(counts, num_bins):
    # lower count = higher weight
    weights = 1.0 / (counts + 0.1)
    probabilities = weights / np.sum(weights)
    bin_idx = np.random.choice(range(num_bins), p=probabilities)
    bin_width = (2 * np.pi) / num_bins
    return ra.uniform(bin_idx * bin_width, (bin_idx + 1) * bin_width)

# ####################################################################################

# initial variables
n = 200
stickiness = 0.9
grid = np.zeros((n, n))
num_p = 500
r_max = 0

center = n // 2
grid[center, center] = 1
spawn = n // 2 - 2
delete = n // 2 + 10
num_bins = 100 # more bins, more angles
angle_counts = np.zeros(num_bins)

for i in range(num_p):

    theta = balanced_spawn_angle(angle_counts, num_bins)
    px = int(center + spawn * np.cos(theta))
    py = int(center + spawn * np.sin(theta))

    stuck = False
    while not stuck:
        px, py = move(px, py)
        if px < 0 or px >= n or py < 0 or py >= n:
            theta = balanced_spawn_angle(angle_counts, num_bins)
            px = int(center + spawn * np.cos(theta))
            py = int(center + spawn * np.sin(theta))
            continue

        dx, dy = px - center, py - center
        if (dx**2 + dy**2) > delete**2:
            theta = balanced_spawn_angle(angle_counts, num_bins)
            px = int(center + spawn * np.cos(theta))
            py = int(center + spawn * np.sin(theta))
            continue

        if sticking(px, py, grid, stickiness):
            grid[px, py] = i + 1
            stuck = True

            # updating the angle
            angle = np.arctan2(py - center, px - center)
            if angle < 0: angle += 2 * np.pi
            bin_idx = int(angle / (2 * np.pi) * num_bins) % num_bins
            angle_counts[bin_idx] += 1

            # spawn and delete radius
            dist = np.sqrt((px - center)**2 + (py - center)**2)
            if dist > r_max:
                r_max = dist
                spawn = r_max + 15
                delete = spawn + 15

# plot
print("")
plt.figure(figsize=(8,8))
plt.imshow(grid, cmap='magma')
plt.title(f"Balanced DLA Fractal (Particles: {num_p})")
plt.axis('off')
plt.show()

# ####################################################################################
# EXTENSION 2 ############################################################
# ###################################################################################

def move_2Dt(px, py):
    # aimless wandering
    direction = ra.randint(0, 5)
    if direction == 0:
      px += 1
    elif direction == 1:
      px -= 1
    elif direction == 2:
      py += 1
    elif direction == 3:
      py -= 1
    elif direction == 4: # diagonal 1
      py += 1
      px += (1 if py % 2 == 1 else -1)
    elif direction == 5: # diagonal 2
      py -= 1
      px += (1 if py % 2 == 1 else -1)

    return px, py

def sticking_2Dt(px, py, grid, stickiness):
  if py % 2 == 1:
    dx = 1
  else:
    dx = -1

  neighbors = [(px + 1, py), (px - 1, py), (px, py + 1), (px, py - 1), (px + dx, py + 1), (px + dx, py - 1)]
  for i, j in neighbors:
        if 0 <= i < grid.shape[0] and 0 <= j < grid.shape[1]:
            if grid[i, j] > 0:
                return ra.random() <= stickiness
  return False

# ##################################################################################

# resetting inital variables
n = 200
grid = np.zeros((n, n))
center = n // 2
grid[center, center] = 1
num_p = 500
# spawn, delete = n // 2 - 10, n // 2 + 2

num_bins = 100 # ADDED
angle_counts = np.zeros(num_bins)
spawn = 10
delete = spawn + 10
r_max = 0

for i in range(num_p):
    theta = ra.uniform(0, 2 * np.pi)
    px, py = int(center + spawn * np.cos(theta)), int(center + spawn * np.sin(theta))
    stuck = False

    while not stuck:
        px, py = move_2Dt(px, py)
        # if it hits the boundary
        if px < 1 or px >= n-1 or py < 1 or py >= n-1:
            px, py = int(center + spawn * np.cos(theta)), int(center + spawn * np.sin(theta))
            continue

        if sticking_2Dt(px, py, grid, stickiness):
            grid[px, py] = i + 1
            stuck = True
            dist = np.sqrt((px - center)**2 + (py - center)**2)
            if dist + 2 > spawn: spawn = dist + 2
            delete = spawn + 5

# plot
print("")
plt.imshow(grid, cmap='magma')
plt.title(f"Balanced Triangular DLA Fractal (Particles: {num_p})")
plt.show()

# ###################################################################################
# EXTENSION 3 #############################################################
# ###################################################################################

def move_3D(px, py, pz):
    direction = ra.randint(0, 5)
    if direction == 0: px += 1
    elif direction == 1: px -= 1
    elif direction == 2: py += 1
    elif direction == 3: py -= 1
    elif direction == 4: pz += 1
    elif direction == 5: pz -= 1
    return px, py, pz

def balanced_spawn_sphere(counts, num_bins_theta, num_bins_phi): # BETA VERSION NOT FULLY IMPLEMENTED
    weights = 1.0 / (counts + 0.1)
    probabilities = weights.flatten() / np.sum(weights)

    idx = np.random.choice(range(num_bins_theta * num_bins_phi), p=probabilities)
    t_idx = idx // num_bins_phi
    p_idx = idx % num_bins_phi

    theta_width = np.pi / num_bins_theta
    phi_width = (2 * np.pi) / num_bins_phi

    theta = ra.uniform(t_idx * theta_width, (t_idx + 1) * theta_width)
    phi = ra.uniform(p_idx * phi_width, (p_idx + 1) * phi_width)

    return theta, phi, t_idx, p_idx

# ####################################################################################

n = 100
grid = np.zeros((n, n, n))
center = n // 2
grid[center, center, center] = 1
num_p = 1000
spawn, stickiness = 5, 0.9

for i in range(num_p):
    # spawn in a sphere
    phi = ra.uniform(0, 2 * np.pi)
    costheta = ra.uniform(-1, 1)
    theta = np.arccos(costheta)
    px = int(center + spawn * np.sin(theta) * np.cos(phi))
    py = int(center + spawn * np.sin(theta) * np.sin(phi))
    pz = int(center + spawn * np.cos(theta))
    stuck = False

    while not stuck:
        px, py, pz = move_3D(px, py, pz)

        # if it hits the wall
        if px < 1 or px >= n-1 or py < 1 or py >= n-1 or pz < 1 or pz >= n-1:
            break

        if np.any(grid[px-1:px+2, py-1:py+2, pz-1:pz+2] > 0):
            if ra.random() <= stickiness:
                grid[px, py, pz] = 1
                stuck = True
                dist = np.sqrt((px-center)**2 + (py-center)**2 + (pz-center)**2)
                if dist + 2 > spawn: spawn = dist + 2

# 3D Plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
z, y, x = np.nonzero(grid)
ax.scatter(x, y, z, c=z, cmap='magma', s=1)
plt.title(f"Balanced DLA Fractal 3D (Particles: {num_p})")
plt.show()


# ####################################################################################
# END PROJECT 4
# ###################################################################################
