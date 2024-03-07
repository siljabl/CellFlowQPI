import sys
import numpy as np
import skimage.io as io
import matplotlib.pyplot as plt

folder = str(sys.argv[1])
f = 1
SzW = [1, 5, 10, 20, 50]

# Set files and paths
file   = f"Well1-1_resc_reg2_0-{f}.tif"
data_path = "../../Data/SyntheticTestData/" + folder
save_path = data_path + "/dic/plots/"
plot_name = folder + f"_SzW{str(SzW)}"


# Import data
u_arr, v_arr = [], []
u_arr.append(np.loadtxt(data_path + "/tif/x_velocity_tif.txt"))
v_arr.append(np.loadtxt(data_path + "/tif/y_velocity_tif.txt"))

# N.B. these are actually displacements and not velocities. For some reason the average displacement does not increase significantly with distance between frames.
for szw in SzW:
    u_DIC = io.imread(f"{data_path}/dic/x_displacement_SzW_{szw}/{file}")
    v_DIC = io.imread(f"{data_path}/dic/y_displacement_SzW_{szw}/{file}")

    u_arr.append(u_DIC)
    v_arr.append(v_DIC)

# Plotting
fig, ax = plt.subplots(2, 6, sharex=True, sharey=True, figsize=(11, 5))
#fig.suptitle(f"DIC window size: {SzW}")

ax[0,0].set(ylabel="x-component", title="True field")
ax[1,0].set(ylabel="y-component")

ax[0,0].imshow(u_arr[0], vmin=np.min(u_arr), vmax=np.max(u_arr))
ax[1,0].imshow(v_arr[0], vmin=np.min(v_arr), vmax=np.max(v_arr))

# Plotting fields
i = 1
for szw in SzW:
    ax[0,i].set(title=f"DIC, SzW: {szw}")

    im_u = ax[0,i].imshow(u_arr[i], vmin=np.min(u_arr), vmax=np.max(u_arr))
    im_v = ax[1,i].imshow(v_arr[i], vmin=np.min(v_arr), vmax=np.max(v_arr))

    i += 1

# Colorbar
fig.subplots_adjust(right=0.9)
cbar_ax_u = fig.add_axes([0.925, 0.54, 0.02, 0.35])
cbar_ax_v = fig.add_axes([0.925, 0.11, 0.02, 0.35])
fig.colorbar(im_u, cax=cbar_ax_u)
fig.colorbar(im_v, cax=cbar_ax_v)

fig.savefig(save_path + f'im_fields_SzW.png')
