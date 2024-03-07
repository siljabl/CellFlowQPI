import os
import sys
import numpy as np
import skimage.io as io
import matplotlib as mpl
import matplotlib.pyplot as plt

folder = str(sys.argv[1])
frame = 1
eps = 1e-12

# Set files and paths
file   = f"Well1-1_resc_reg2_0-{str(frame)}.tif"
data_path = "../../Data/SyntheticTestData/" + folder
save_path = data_path + "/plots/"
if os.path.isdir(save_path) == 0:
    os.mkdir(save_path)


# take all folders in folder
SzW  = np.array([1, 2, 5, 10, 20, 50])
umax = np.array([1, 2, 4, 6, 8, 12, 16])
df = 1

# empty arrays for computations
u_err = np.zeros([len(umax), len(SzW)])
v_err = np.zeros([len(umax), len(SzW)])

u_mean = np.zeros([len(umax), len(SzW)])
v_mean = np.zeros([len(umax), len(SzW)])

u_true_mean = np.zeros([len(umax)])
v_true_mean = np.zeros([len(umax)])


i = 0
for u in umax:
    # Importing true values
    u_true = np.loadtxt(data_path + f"/umax_{u}/tif/x_velocity_tif.txt")
    v_true = np.loadtxt(data_path + f"/umax_{u}/tif/y_velocity_tif.txt")

    u_true_mean[i] = np.mean(u_true)
    v_true_mean[i] = np.mean(v_true)
    
    j = 0
    for szw in SzW:
        # import datafolders
        x_DIC = io.imread(f"{data_path}/umax_{u}/dic/x_displacement_SzW_{szw}/{file}")
        y_DIC = io.imread(f"{data_path}/umax_{u}/dic/y_displacement_SzW_{szw}/{file}")

        # why does displacement not increase with frame number?
        u_DIC = x_DIC / df
        v_DIC = y_DIC / df

        # mean values
        u_mean[i,j] = np.mean(u_DIC)
        v_mean[i,j] = np.mean(v_DIC)       

        # rms
        u_err[i,j] = np.sqrt(np.mean(((u_true - u_DIC) ** 2)))
        v_err[i,j] = np.sqrt(np.mean(((v_true - v_DIC) ** 2)))

        j += 1
    i += 1

tot_err  = np.sqrt(u_err**2  + v_err**2)
tot_mean = np.sqrt(u_mean**2 + v_mean**2)
tot_true_mean = np.sqrt(u_true_mean**2 + v_true_mean**2)


titles = [r"Total values, as $\sqrt{u_{value}^2 + v_{value}^2}$", \
          "x-component", "y-component"]
means  = [tot_mean, u_mean, v_mean]
true_means = [tot_true_mean, u_true_mean, v_true_mean]
errors = [tot_err,  u_err,  v_err] 
names  = ["total", "x", "y"] 

u_cmap = mpl.cm.plasma
s_cmap = mpl.cm.viridis

u_colors = [u_cmap(i/len(umax)) for i in range(len(umax))];
s_colors = [s_cmap(i/len(SzW))  for i in range(len(SzW))];

u_plot = [0, 2, 3, 4]
s_plot = [0, 2, 3, 4]

for title, mean, true_mean, err, name in zip(titles, means, true_means, errors, names):
    fig, ax = plt.subplots(2, 2, figsize=(8,6))
    fig.suptitle(title)

    for u, s in zip(u_plot, s_plot):
        ax[0,0].plot(umax, mean[:,s], 'o-', color=u_colors[u], label=f"SzW = {SzW[s]}")
        ax[0,1].plot(umax, err[:,s] / umax,  'o-', color=u_colors[u], label=f"SzW = {SzW[s]}")

        ax[1,0].plot(SzW, mean[u], 'o-', color=s_colors[s], label=f"umax = {umax[u]}")
        ax[1,0].plot(SzW,  true_mean[u] * np.ones_like(SzW),  '--', color=s_colors[s])
        ax[1,1].plot(SzW, err[u] / umax[u], 'o-', color=s_colors[s], label=f"umax = {umax[u]}")

    ax[0,0].plot(umax, true_mean, '--', color="gray", label=f"SzW = {SzW[s]}")
    ax[0,0].plot(umax, umax * true_mean[0], 'x', color="gray", label=f"SzW = {SzW[s]}")

    ax[0,0].set(xlabel=r"$u_{max}$", ylabel="Mean value")
    ax[0,1].set(xlabel=r"$u_{max}$", ylabel=r"RMS / $u_{max}$")
    ax[1,0].set(xlabel=r"DIC window size", ylabel="Mean value")
    ax[1,1].set(xlabel=r"DIC window size", ylabel=r"RMS / $u_{max}$")

    for axes in ax.flatten():
        axes.set(xscale="log")
    
    for i in range(2):
        ax[i,1].legend()

    #fig.subplots_adjust(right=0.8)
    fig.tight_layout()#rect=[0, 0, 0.85, 1])
    fig.savefig(save_path + name + "_DIC_performance.png")


# plot imshow...
# fig, ax = plt.subplots(1, 2, figsize=(6,3))
# im_mean = ax[0].imshow(tot_mean)
# im_rms  = ax[1].imshow(tot_err)

# ax[0].set(title="mean", xlabel="Window size", ylabel="U_max")
# ax[1].set(title="rms", xlabel="Window size", ylabel="U_max")

# fig.tight_layout()

# fig.colorbar(im_mean, ax=ax[0])
# fig.colorbar(im_rms,  ax=ax[1])

# fig.savefig(save_path + "DIC_performance_im.png")