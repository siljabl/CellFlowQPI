import os
import numpy as np
import skimage.io as io
import matplotlib.pyplot as plt

folder = "../Data/SyntheticTestData/Basic1"
file   = "Well1-1_resc_reg2_0"
SzW = 5
df = 1
f_min = 1
f_max = 10

u = np.loadtxt(folder + "/tif/x_velocity_tif.txt")
v = np.loadtxt(folder + "/tif/y_velocity_tif.txt")

# normalize
#u = u / np.sqrt(u**2 + v**2).max()
#v = v / np.sqrt(u**2 + v**2).max()

N = np.shape(u)[0] * np.shape(u)[1]

u_err = np.zeros(f_max - f_min)
v_err = np.zeros(f_max - f_min)

for i in range(1, f_max):
    # import data
    x_DIC = io.imread(f"{folder}/dic/x_displacement_SzW_{SzW}/{file}-{str(i)}.tif")
    y_DIC = io.imread(f"{folder}/dic/y_displacement_SzW_{SzW}/{file}-{str(i)}.tif")

    # why does displacement not increase with frame number?
    u_DIC = x_DIC / df
    v_DIC = y_DIC / df

    # normalize
    #u_DIC = u_DIC / np.sqrt(u_DIC**2 + v_DIC**2).max()
    #v_DIC = v_DIC / np.sqrt(u_DIC**2 + v_DIC**2).max()

    # least squares
    u_err[i-1] = np.sum((u - u_DIC)**2 / abs(u + 0.1)) / N
    v_err[i-1] = np.sum((v - v_DIC)**2 / abs(v + 0.1)) / N

frames = np.arange(f_min, f_max)
plt.plot(frames, u_err, 'x-',  label="u_err")
plt.plot(frames, v_err, 'x-', label="v_err")
plt.legend()
plt.savefig(f"Silja/basic1_err_{SzW}.png")

plt.imshow(u_DIC)
plt.savefig("Silja/plot.png")
