import numpy as np
import matplotlib.pyplot as plt

dir = "../Data/Synthetic/Field4/"
u = np.loadtxt(dir + "x_velocity.txt")
v = np.loadtxt(dir + "y_velocity.txt")

step_size = 10

fig, ax = plt.subplots(1,1, figsize=(20,20))
ax.quiver(u[0:-1:step_size, 0:-1:step_size],v[0:-1:step_size, 0:-1:step_size], width=0.0015)
fig.tight_layout()
fig.savefig(dir + "velocity_field.pdf")