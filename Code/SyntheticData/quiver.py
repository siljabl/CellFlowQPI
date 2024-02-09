import numpy as np
import matplotlib.pyplot as plt

#dir = "../../Data/Synthetic/Field4/"
in_dir = "../../Data/InitialConditions/"
out_dir = ""
u = np.loadtxt(in_dir + "x_velocity.txt").T
v = np.loadtxt(in_dir + "y_velocity.txt").T

step_size = 1

fig, ax = plt.subplots(1,1, figsize=(20,20))
ax.quiver(u[0:-1:step_size, 0:-1:step_size],v[0:-1:step_size, 0:-1:step_size], width=0.0015)
fig.tight_layout()
fig.savefig(out_dir + "velocity_field.png")