import os
import sys
import numpy as np
import skimage.io as io
import matplotlib.pyplot as plt

folder = str(sys.argv[1])
f = 6
path = "../Data/SyntheticTestData/" + folder
file   = "Well1-1_resc_reg2_0"
f_min  = 1
f_max  = 10


SzW = np.array([1, 2, 5, 10, 20])
df = 1
folder = sys.argv[1]
u = np.loadtxt(path + "/tif/x_velocity_tif.txt")
v = np.loadtxt(path + "/tif/y_velocity_tif.txt")
N = np.shape(u)[0] * np.shape(u)[1]

# normalize
#u = u / np.sqrt(u**2 + v**2).max()
#v = v / np.sqrt(u**2 + v**2).max()

# Compute error
frames = np.arange(f_min, f_max)

u_err_SzW = np.zeros_like(SzW)
v_err_SzW = np.zeros_like(SzW)

j = 0
for szw in SzW:
    plot_name = folder + f"_SzW{str(szw)}"

    u_err_frame = np.zeros(f_max - f_min)
    v_err_frame = np.zeros(f_max - f_min)

    for i in range(f_min, f_max):
        # import datafolder = sys.argv[1]
        x_DIC = io.imread(f"{path}/dic/x_displacement_SzW_{szw}/{file}-{str(i)}.tif")
        y_DIC = io.imread(f"{path}/dic/y_displacement_SzW_{szw}/{file}-{str(i)}.tif")

        # why does displacement not increase with frame number?
        u_DIC = x_DIC / df
        v_DIC = y_DIC / df

        # normalize
        #u_DIC = u_DIC / np.sqrt(u_DIC**2 + v_DIC**2).max()
        #v_DIC = v_DIC / np.sqrt(u_DIC**2 + v_DIC**2).max()

        # least squares
        u_err_frame[i-1] = np.sum((u - u_DIC)**2 / abs(u + 0.1)) / N
        v_err_frame[i-1] = np.sum((v - v_DIC)**2 / abs(v + 0.1)) / N

    fig, ax = plt.subplots(1,1)
    ax.plot(frames, u_err_frame, 'x-', label="u_err")
    ax.plot(frames, v_err_frame, 'x-', label="v_err")
    ax.legend()
    fig.savefig(f"Silja/"+plot_name+"_err.png")

    u_err_SzW[j] = u_err_frame[f] #np.mean(u_err_frame)
    v_err_SzW[j] = v_err_frame[f] #np.mean(v_err_frame)
    j += 1


fig, ax = plt.subplots(1,1)
ax.plot(SzW, u_err_SzW, 'x-', label="u_err")
ax.plot(SzW, v_err_SzW, 'x-', label="v_err")
ax.set(xlabel="Window size")
ax.legend()
fig.savefig(f"Silja/" + folder + "_err.png")
